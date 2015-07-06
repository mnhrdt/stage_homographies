//gcc-5 -fopenmp -O3 viho.c -I/usr/local/include/libiomp -I/usr/X11/include -I/Users/Hhhh/ENS/Stage_L3_math/homographies/code/jpeg-6b -L/usr/X11/lib -lfftw3 -lX11 -L/usr/local/Cellar/libtiff/4.0.3 -ltiff -ljpeg -lpng


#include "iio.c"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "decomp.h"




/**
  * this program warps an image by an homography
  * the warping is such that the quadrilateral x,y,z,t
  * (supposed convex) is mapped to a rectangle such that
  * height = ratio * width
  * this rectangle is noted
  * A   B
  * D   C
  *
  * the area of the rectangle and the area of the
  * quadrilateral are equal
  */



// compute the inverse homography (inverse of a 3x3 matrix)
double invert_homography(double invH[3][3], double H[3][3]){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = (a[4]*a[8]-a[5]*a[7])/det;
	r[1] = (a[2]*a[7]-a[1]*a[8])/det;
	r[2] = (a[1]*a[5]-a[2]*a[4])/det;
	r[3] = (a[5]*a[6]-a[3]*a[8])/det;
	r[4] = (a[0]*a[8]-a[2]*a[6])/det;
	r[5] = (a[2]*a[3]-a[0]*a[5])/det;
	r[6] = (a[3]*a[7]-a[4]*a[6])/det;
	r[7] = (a[1]*a[6]-a[0]*a[7])/det;
	r[8] = (a[0]*a[4]-a[1]*a[3])/det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
void compose_homographies(double C[3][3], double A[3][3], double B[3][3]){
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}



// Find the homography that changes the canonical projective basis into the given four points (x, y, z, w)
void homography_from_four_points(double H[3][3], double x[2], double y[2], double z[2], double w[2]){
	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;
	if (!isnormal(DET))
		fprintf(stderr, "denormal! DET = %f\n", DET);

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,t) to (a,b,c,d)
void homography_from_eight_points(double H[3][3], double x[2], double y[2], double z[2], double w[2], double a[2], double b[2], double c[2], double d[2]){
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}



//determine if (ab) is a diagonal of the quadrilateral a,b,c,d
bool is_a_diagonal(double a[2], double b[2], double c[2], double d[2]){
    double u,v,w; //(ab) has for equation ux+vy+w=0

    u = b[1]-a[1];
    v = a[0]-b[0];
    w = -u*a[0]-v*a[1];

    bool firstCase = (u*c[0]+v*c[1]+w<0) && (u*d[0]+v*d[1]+w>0),
        secondCase = (u*c[0]+v*c[1]+w>0) && (u*d[0]+v*d[1]+w<0);
    return firstCase || secondCase;
}

//find the two diagonals of the quadrilateral a,b,c,d ; return if it actually found two
bool find_diagonals(double diag1[2][2], double diag2[2][2], double a[2], double b[2], double c[2], double d[2]){
    if(is_a_diagonal(a,b,c,d)){ //diagonals, if there are two, are (ab) and (cd)
        diag1[0][0] = a[0]; diag1[0][1] = a[1];
        diag1[1][0] = b[0]; diag1[1][1] = b[1];
        diag2[0][0] = c[0]; diag2[0][1] = c[1];
        diag2[1][0] = d[0]; diag2[1][1] = d[1];
        return is_a_diagonal(c,d,a,b);
    }else if(is_a_diagonal(a,c,b,d)){ //diagonals, if there are two, are (ac) and (bd)
        diag1[0][0] = a[0]; diag1[0][1] = a[1];
        diag1[1][0] = c[0]; diag1[1][1] = c[1];
        diag2[0][0] = b[0]; diag2[0][1] = b[1];
        diag2[1][0] = d[0]; diag2[1][1] = d[1];
        return is_a_diagonal(b,d,a,c);
    }else{ //diagonals, if there are two, are (ad) and (bc)
        diag1[0][0] = a[0]; diag1[0][1] = a[1];
        diag1[1][0] = d[0]; diag1[1][1] = d[1];
        diag2[0][0] = b[0]; diag2[0][1] = b[1];
        diag2[1][0] = c[0]; diag2[1][1] = c[1];
        return is_a_diagonal(a,d,b,c) && is_a_diagonal(b,c,a,d);
    }
}

//reorder the two diagonals such that the first maps to (AC) and the second maps to (BD)
//A,B,C,D are defined on the head of the file
void reorder_diagonals(double diag1[2][2], double diag2[2][2], double ratio){
    /*
     * to choose which diagonal will maps to (AC), we choose the one "nearer" to (AC),
     * i.e. which forms the smaller angle with AC (modulo pi). The absolute value of
     * the cosine of this angle is an indicator of how much they are near. Thus, the
     * scalar product is also a good indicator, after normalization. Furthermore, it
     * indicates if two vectors have the same orientation or if one must be reverted
     *
     * choosing to compare to (AC) instead of (BD) does not matter, since one is nearer
     * of (AC) and the other is nearer of (BD)
     */
    //copy the coordinates
    double x[2] = {diag1[0][0],diag1[0][1]},
        y[2] = {diag2[0][0],diag2[0][1]},
        z[2] = {diag1[1][0],diag1[1][1]},
        t[2] = {diag2[1][0],diag2[1][1]};

    //compute vectors associated to segments
    double vectAC[2] = {1., ratio},
        vectBD[2] = {-1., ratio},
        vectxz[2] = {z[0]-x[0], z[1]-x[1]},
        vectyt[2] = {t[0]-y[0], t[1]-y[1]}; //vector yt

    //normalize vectors
    double n;
    n = sqrt(pow(vectAC[0],2.) + pow(vectAC[1],2.));
    vectAC[0] /= n; vectAC[1] /= n;
    n = sqrt(pow(vectBD[0],2.) + pow(vectBD[1],2.));
    vectBD[0] /= n; vectBD[1] /= n;
    n = sqrt(pow(vectxz[0],2.) + pow(vectxz[1],2.));
    vectxz[0] /= n; vectxz[1] /= n;
    n = sqrt(pow(vectyt[0],2.) + pow(vectyt[1],2.));
    vectyt[0] /= n; vectyt[1] /= n;

    //compute scalar products with vector AC
    double sp1 = vectAC[0]*vectxz[0]+vectAC[1]*vectxz[1];
    double sp2 = vectAC[0]*vectyt[0]+vectAC[1]*vectyt[1];

    if(fabs(sp1)>fabs(sp2)){ //xz is nearer of AC
        if(sp1>0){ //xz and AC are (quiet) same oriented
            diag1[0][0] = x[0]; diag1[0][1] = x[1];
            diag1[1][0] = z[0]; diag1[1][1] = z[1];
        }else{ //xz and AC are (quiet) opposite oriented
            diag1[0][0] = z[0]; diag1[0][1] = z[1];
            diag1[1][0] = x[0]; diag1[1][1] = x[1];
        }
        if(vectBD[0]*vectyt[0]+vectBD[1]*vectyt[1]>0){ //yt and BD are (quiet) same oriented
            diag2[0][0] = y[0]; diag2[0][1] = y[1];
            diag2[1][0] = t[0]; diag2[1][1] = t[1];
        }else{ //yt and BD are (quiet) opposite oriented
            diag2[0][0] = t[0]; diag2[0][1] = t[1];
            diag2[1][0] = y[0]; diag2[1][1] = y[1];
        }
    }else{ //yt is nearer of AC
        if(sp2>0){ //yt and AC are (quiet) same oriented
            diag1[0][0] = y[0]; diag1[0][1] = y[1];
            diag1[1][0] = t[0]; diag1[1][1] = t[1];
        }else{ //yt and AC are (quiet) opposite oriented
            diag1[0][0] = t[0]; diag1[0][1] = t[1];
            diag1[1][0] = y[0]; diag1[1][1] = y[1];
        }
        if(vectBD[0]*vectxz[0]+vectBD[1]*vectxz[1]>0){ //xz and BD are (quiet) same oriented
            diag2[0][0] = x[0]; diag2[0][1] = x[1];
            diag2[1][0] = z[0]; diag2[1][1] = z[1];
        }else{ //xz and BD are (quiet) opposite oriented
            diag2[0][0] = z[0]; diag2[0][1] = z[1];
            diag2[1][0] = x[0]; diag2[1][1] = x[1];
        }
    }
}



//compute the area of a quadrilateral by separating it in 2 triangles
double area_of_quadrilateral(double a[2], double b[2], double c[2], double d[2]){
    double ab = sqrt(pow(b[0]-a[0],2.)+pow(b[1]-a[1],2.)),
        bc = sqrt(pow(c[0]-b[0],2.)+pow(c[1]-b[1],2.)),
        cd = sqrt(pow(d[0]-c[0],2.)+pow(d[1]-c[1],2.)),
        da = sqrt(pow(a[0]-d[0],2.)+pow(a[1]-d[1],2.)),
        ac = sqrt(pow(c[0]-a[0],2.)+pow(c[1]-a[1],2.));

    double s, area1, area2;
    s = (ab+bc+ac)/2.; //semi sum of sides of the triangle
    area1 = sqrt(s*(s-ab)*(s-bc)*(s-ac));
    s = (cd+da+ac)/2.; //semi sum of sides of the triangle
    area2 = sqrt(s*(s-cd)*(s-da)*(s-ac));

    return area1+area2;
}



int main(int argc,char *argv[]){
	if (argc != 11) {
		printf("usage :\n\t[image.png] x0 x1 y0 y1 z0 z1 t0 t1 ratio\n");
		return 1;
	}



	//image input
	char *filename_in = argv[1];
	float *img;
	int w,h,pd;
	img = iio_read_image_float_vec(filename_in, &w, &h, &pd);



	//coefficients input
	double x[2] = {strtod(argv[2],NULL), strtod(argv[3],NULL)},
        y[2] = {strtod(argv[4],NULL), strtod(argv[5],NULL)},
        z[2] = {strtod(argv[6],NULL), strtod(argv[7],NULL)},
        t[2] = {strtod(argv[8],NULL), strtod(argv[9],NULL)},
        ratio = strtod(argv[10],NULL);

    //reorder x,y,z,t in a,b,c,d to get the simplest homography
    double diag1[2][2], diag2[2][2];
    if(!find_diagonals(diag1, diag2, x, y, z, t)){
        printf("The quadrilateral is not convex. Please select a convex quadrilateral\n");
        exit(1);
    };
    reorder_diagonals(diag1,diag2,ratio);
    double a[2] = {diag1[0][0],diag1[0][1]},
        b[2] = {diag2[0][0],diag2[0][1]},
        c[2] = {diag1[1][0],diag1[1][1]},
        d[2] = {diag2[1][0],diag2[1][1]};

    //get the homography
    double area = area_of_quadrilateral(a,b,c,d);
    int w_out = round(sqrt(area/ratio)),
        h_out = round(sqrt(area*ratio));
    double A[2] = {0,0},
        B[2] = {w_out,0},
        C[2] = {w_out,h_out},
        D[2] = {0,h_out};

	double H[3][3];
	homography_from_eight_points(H,A,B,C,D,a,b,c,d);
        //this H is such that H(A)=a, H(B)=b,...
        //but apply_homo_final resamples by the inverse homography



    //resample (and time)
	float *img_f = malloc(3*w_out*h_out*sizeof(float));
	clock_t debutcpu,fincpu;
	double debutreal,finreal;
	debutcpu = clock();
	debutreal = omp_get_wtime();
	if(pd==3){
        apply_homo_final(img,img_f,w,h,w_out,h_out,H);
	}else{//suppose pd=1
        float *img3 = malloc(3*w*h*sizeof(float));
        for(int i=0;i<w*h;i++){
            for(int l = 0;l<3;l++){
                img3[3*i+l]=img[i];
            }
        }
        apply_homo_final(img3,img_f,w,h,w_out,h_out,H);
	}

	fincpu = clock();
	finreal = omp_get_wtime();
	printf("cputime :%fs\ntime : %fs\n",(double)(fincpu-debutcpu)/CLOCKS_PER_SEC,(double)(finreal-debutreal));



	//output
	iio_save_image_float_vec("img_f.png",img_f,w_out,h_out,3);
	free(img);
	free(img_f);



	return 0;
}
