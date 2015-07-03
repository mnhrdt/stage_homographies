//gcc-5 -fopenmp -O3 viho.c -I/usr/local/include/libiomp -I/usr/X11/include -I/Users/Hhhh/ENS/Stage_L3_math/homographies/code/jpeg-6b -L/usr/X11/lib -lfftw3 -lX11 -L/usr/local/Cellar/libtiff/4.0.3 -ltiff -ljpeg -lpng


#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "decomp.h"



/**
  * this program warps an image by an homography
  * the warping is such that
  *
  * - maps to -
  * x   -->   a
  * y   -->   b
  * z   -->   c
  * t   -->   d
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
		fprintf(stderr, "denormal! DET = %g\n", DET);

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
void homography_from_eight_points(double H[3][3], double x[2], double y[2], double z[2], double w[2], double a[2], double b[2], double c[2], double d[2])
{
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
		return 0;
	if (l < 0) l = 0;
	if (l >= pd) l = pd - 1;
	return x[(w*j+i)*pd+l];
}

// ultra-naive warp
static void warp_homography_nn(float *y, int yw, int yh, double H[3][3],
		float *x, int w, int h, int pd)
{
	//double invH[3][3];
	//invert_homography(invH, H);

	fprintf(stderr, "yw = %d\n", yw);
	fprintf(stderr, "yh = %d\n", yh);
	fprintf(stderr, "w = %d\n", w);
	fprintf(stderr, "h = %d\n", h);

	for (int j = 0; j < yh; j++)
	for (int i = 0; i < yw; i++)
	for (int l = 0; l < pd; l++)
	{
		double ypos[2] = {i, j};
		double xpos[2]; apply_homography(xpos, H, ypos);
		int ii = round(xpos[0]-0.5);
		int jj = round(xpos[1]+0.5);
		int idx = (j * yw + i) * pd + l;
		y[idx] = getsample_0(x, w, h, pd, ii, jj, l);
	}
}

static void warp_homography(float *img, float *img_f, int w, int h, int pd,
		int WOUT, int HOUT, double H[3][3], int method_id)
{
	if (method_id == 1) { // decomposition method
		if (pd == 3) {
			apply_homo_final(img,img_f,w,h,WOUT,HOUT,H);
		} else { //suppose pd=1
			assert(pd == 1);
			float *img3 = malloc(3*w*h*sizeof(float));
			for(int i=0;i<w*h;i++){
				for(int l = 0;l<3;l++){
					img3[3*i+l] = img[i];
				}
			}
			apply_homo_final(img3,img_f,w,h,WOUT,HOUT,H);
		}
	} else if (method_id == 0) { // nearest neighbor interpolation
		warp_homography_nn(img_f, WOUT, HOUT, H, img, w, h, pd);
	} else
		exit(fprintf(stderr,"unrecognized method_id %d\n", method_id));
}


#include "iio.c"
#include "pickopt.c"
int main(int argc,char *argv[])
{
	int method_id = atoi(pick_option(&argc, &argv, "m", "0"));
	if (argc != 21) {
		fprintf(stderr, "usage:\n\t%s i.png o.png w h "
		//                          0 1     2     3 4
			"x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1\n",
		//       5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
			*argv);
		return 1;
	}

	// read input parameters
	char *filename_in = argv[1];
	char *filename_out = argv[2];
	int WOUT = atoi(argv[3]);
	int HOUT = atoi(argv[4]);
	double from[4][2], to[4][2];
	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 2; j++)
	{
		from[i][j] = atof(argv[5  + 2*i + j]);
		to[i][j]   = atof(argv[13 + 2*i + j]);
	}

	// compute the coefficients of the requested homography
	// note: apply_homo_final resamples by the inverse homography
	double H[3][3];
	homography_from_eight_points(H,
			from[0], from[1], from[2], from[3],
			to[0], to[1], to[2], to[3]
			);

	// read input image
	int w,h,pd;
	float *img = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *img_f = malloc(3*WOUT*HOUT*sizeof(float));

	// start timer
	clock_t debutcpu,fincpu;
	double debutreal,finreal;
	debutcpu = clock();
	debutreal = omp_get_wtime();

	// call the actual algorithm
	warp_homography(img, img_f, w, h, pd, WOUT, HOUT, H, method_id);

	// stop timer and print the running time
	fincpu = clock();
	finreal = omp_get_wtime();
	fprintf(stderr, "cputime :%fs\ntime : %fs\n",
			(double)(fincpu-debutcpu)/CLOCKS_PER_SEC,
			(double)(finreal-debutreal));

	// save output image
	iio_save_image_float_vec(filename_out ,img_f, WOUT, HOUT, 3);

	// cleanup and exit
	free(img);
	free(img_f);
	return 0;
}
