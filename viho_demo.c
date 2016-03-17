//gcc -std=c99 -fopenmp -O3 viho_demo.c -lfftw3 -ltiff -ljpeg -lpng


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


// y = H(x)
static void apply_homography_1pt(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

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

// type of the "extrapolator" functions
typedef float (*extrapolator_t)(float*,int,int,int,int,int,int);

// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int n, int p)
{
	int r = n % p;
	r = r < 0 ? r + p : r;
	assert(r >= 0);
	assert(r < p);
	return r;
}

// instance of "extrapolator_t", extrapolate by periodicity
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// instance of "extrapolator_t", extrapolate by a constant value
static float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return value;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}


// type of the "interpolator" functions
typedef float (*interpolator_t)(float*,int,int,int,float,float,int,
		extrapolator_t);

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// instance of "interpolator_t", for bilinear interpolation
static float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}

// instance of "interpolator_t" for nearest neighbor interpolation
static float nearest_neighbor_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = round(p);
	int iq = round(q);
	return pix(x, w, h, pd, ip, iq, l);
}


// one-dimensional cubic interpolation of four data points ("Keys")
static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

// two-dimensional separable cubic interpolation, on a 4x4 grid
static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

// instance of "interpolator_t" for bicubic interpolation
static float bicubic_interpolation_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}

static extrapolator_t obtain_extrapolator(int id)
{
	if (id == -1)
		return getsample_per;
	float top = id;
	getsample_cons(&top, 0, 0, 0, 0, 0, 0);
	return getsample_cons;
}

static interpolator_t obtain_interpolator(int id)
{
	if (id == 0) return nearest_neighbor_at;
	if (id == 2) return bilinear_interpolation_at;
	if (id == 3) return bicubic_interpolation_at;
	return nearest_neighbor_at;
}

static void warp_homography_generic(float *y, int yw, int yh, double H[3][3],
		float *x, int w, int h, int pd,
		int extrapolator, int interpolator)
{
	extrapolator_t OUT  = obtain_extrapolator(extrapolator);
	interpolator_t EVAL = obtain_interpolator(interpolator);

	for (int j = 0; j < yh; j++)
	for (int i = 0; i < yw; i++)
	{
		double p[2] = {i, j};
		apply_homography_1pt(p, H, p);
		//p[0] += 1.0;
		//p[1] += 1.0;
		//p[0] = p[0] * w / (w - 1.0) - 0.5;
		//p[1] = p[1] * h / (h - 1.0) - 0.5;
		for (int l = 0; l < pd; l++)
		{
			int idx = l + pd * (yw * j + i);
			float v = EVAL(x, w, h, pd, p[0], p[1], l, OUT);
			y[idx] = v;
		}
	}
}

static bool insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

// draw the sample positions
static void warp_homography_dots(float *y, int yw, int yh, double H[3][3],
		float *x, int w, int h, int pd)
{
	double invH[3][3];
	invert_homography(invH, H);

	// fill background
	for (int i = 0; i < yw * yh * pd; i++)
		y[i] = 0;

	// traverse samples
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
	{
		double xpos[2] = {i, j};
		double ypos[2]; apply_homography_1pt(ypos, invH, xpos);
		int ii = round(ypos[0]);
		int jj = round(ypos[1]);
		int xidx = (j * w + i) * pd + l;
		int yidx = (jj * yw + ii) * pd + l;
		if (insideP(yw, yh, ii, jj))
			y[yidx] = x[xidx];
	}
}

static void warp_homography(float *img, float *img_f, int w, int h, int pd,
		int ow, int oh, double H[3][3], int method_id)
{
	if (method_id == 1) { // decomposition method
		exit(fprintf(stderr, "decomposintion not here"));
	} else if (method_id == -1) { // positions of the samples
		warp_homography_dots(img_f, ow, oh, H, img, w, h, pd);
	} else if (method_id == 0) { // nearest neighbor interpolation
		warp_homography_generic(img_f, ow,oh, H, img, w,h,pd, 0, 1);
	} else if (method_id == 2) { // bilinear interpolation
		warp_homography_generic(img_f, ow,oh, H, img, w,h,pd, 0, 2);
	} else if (method_id == 3) { // bicubic interpolation
		warp_homography_generic(img_f, ow,oh, H, img, w,h,pd, 0, 3);
	} else
		exit(fprintf(stderr,"unrecognized method_id %d\n", method_id));
}


#include "iio.c"
#include "pickopt.c"
int main(int argc,char *argv[])
{
	int method_id      = atoi(pick_option(&argc, &argv, "m", "0"));
	global_ZOOM_INPUT  = atof(pick_option(&argc, &argv, "s", "0.7"));
	global_ZOOM_OUTPUT = atof(pick_option(&argc, &argv, "S", "0.7"));
	if (argc != 21) {
		fprintf(stderr, "usage:\n%s in.png out.png w h "
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
//	if(pd==3){
//        apply_homography(img,img_f,w,h,WOUT,HOUT,H);
//	}else{//suppose pd=1
//        float *img3 = malloc(3*w*h*sizeof(float));
//        for(int i=0;i<w*h;i++){
//            for(int l = 0;l<3;l++){
//                img3[3*i+l]=img[i];
//            }
//        }
//        apply_homography(img3,img_f,w,h,WOUT,HOUT,H);
//	}

	// call the actual algorithm
	omp_set_num_threads(16);
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
