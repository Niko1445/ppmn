#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"

struct Args{
	double N;
	double Tc;
	double Tr;
};

void A32(
 struct Args *args,
 void f(struct Args *args, double t, gsl_vector* y, gsl_vector *k),
 int i,
 double h,
 double t,
 gsl_vector *y,
 gsl_vector *yn,
 gsl_vector *k,
 gsl_vector *err,
 gsl_matrix *K)
 {
	gsl_vector *k1 = gsl_vector_alloc(K->size1);
	gsl_vector *k2 = gsl_vector_alloc(K->size1);
	gsl_vector *k3 = gsl_vector_alloc(K->size1);

	switch (i) {
		case 0: ;
			f(args, t, y, k);
			break;
		case 1: ;
			gsl_matrix_get_col(k1, K, 0);
			gsl_vector_scale(k1, h*1./2.);
			gsl_vector_add(k1, y);
			f(args, t + 1./2.*h, k1, k);
			break;
		case 2: ;
			gsl_matrix_get_col(k2, K, 1);
			gsl_vector_scale(k2, h*3./4.);
			gsl_vector_add(k2, y);
			f(args, t + 3./4.*h, k2, k);
			break;
		case 3: ;
			gsl_matrix_get_col(k1, K, 0);
			gsl_matrix_get_col(k2, K, 1);
			gsl_matrix_get_col(k3, K, 2);

			gsl_vector_scale(k1, h*2./9.);
			gsl_vector_scale(k2, h*1./3.);
			gsl_vector_scale(k3, h*4./9.);

			gsl_vector_add(k1, k2);
			gsl_vector_add(k1, k3);
			gsl_vector_add(k1, y);

			gsl_vector_memcpy(yn,k1);
			f(args, t+h, yn, k);

			gsl_matrix_get_col(k1, K, 0);
			gsl_matrix_get_col(k2, K, 1);
			gsl_matrix_get_col(k3, K, 2);
			gsl_vector_memcpy(err,k);

			gsl_vector_scale(k1, h*7./24.);
			gsl_vector_scale(k2, h*1./4.);
			gsl_vector_scale(k3, h*1./3.);
			gsl_vector_scale(err, h*1./8.);

			gsl_vector_add(err, k1);
			gsl_vector_add(err, k2);
			gsl_vector_add(err, k3);
			gsl_vector_add(err, y);

			gsl_vector_memcpy(k1,yn); // err stores z_n+1 yn stores y_n+1
			gsl_vector_sub(k1,err);
			gsl_vector_memcpy(err,k1); // err now stores the error vector
			break;
	}
	gsl_vector_free(k1);
	gsl_vector_free(k2);
	gsl_vector_free(k3);

}

void rkstep32(
	struct Args *args,
	void f(struct Args *args, double t, gsl_vector* y, gsl_vector *k),
  double h,
  double t,
  gsl_vector *y,
  gsl_vector *yn,
	gsl_vector *k,
  gsl_vector *err)
	{
	gsl_matrix *K = gsl_matrix_alloc(y->size,4);
	gsl_matrix_set_col(K, 0, k);

	if (gsl_vector_isnull(k)) {
		printf("k is null\n");
		A32(args, f, 0, h, t, y, yn, k, err, K);
		gsl_matrix_set_col(K, 0, k);
	} else {
		/* code */
	}
	for (int i = 1; i < 4; i++) {
		A32(args, f, i, h, t, y, yn, k, err, K);
		gsl_matrix_set_col(K, i, k);
	}
	gsl_matrix_free(K);
}

void driver(
	FILE* data_stream,
	struct Args *args,
	void f(struct Args *args, double t, gsl_vector* y, gsl_vector *k),
	double a,
	gsl_vector *ya,
	double b,
	gsl_vector *yb,
	double acc,
	double eps
) {
	gsl_vector *k = gsl_vector_alloc(ya->size);
	gsl_vector_set_zero(k);
	gsl_vector *err = gsl_vector_alloc(ya->size);
	gsl_vector_set_zero(err);
	gsl_vector *p = gsl_vector_alloc(3);
	double ltol, lerr, h;
	h = 0.025;

	for (double t = a; t < b; t+=h){
		rkstep32(args, f, h, t, ya, yb, k, err);
		ltol = (eps*gsl_vector_length(yb)+acc)*sqrt(h/(b-a));
		lerr = gsl_vector_length(err);
		h = h * pow(ltol/lerr, 1./4.) * 0.95;
		gsl_vector_set(p, 0, t);
		gsl_vector_set(p, 1, ltol);
		gsl_vector_set(p, 2, lerr);

		vector_printEX(data_stream, yb, p);
		//fprintf(data_stream, "%g  %g  %g  %g  %g\n",t, gsl_vector_get(yb, 0), gsl_vector_get(yb, 1), ltol, lerr);
		gsl_vector_memcpy(ya, yb);
	}
	gsl_vector_free(k);
	gsl_vector_free(err);
	gsl_vector_free(p);
}

void harm_osc(struct Args *args, double t, gsl_vector *y, gsl_vector *dxdy){ //u'' + u = 0
	gsl_vector_set(dxdy,0,gsl_vector_get(y,1));
	gsl_vector_set(dxdy,1,-gsl_vector_get(y,0));
}

void SIR(struct Args *args, double t, gsl_vector *y, gsl_vector *dxdy){
	gsl_vector_set(dxdy,0, -(gsl_vector_get(y,0)*gsl_vector_get(y,1))/(args->N*args->Tc));
	gsl_vector_set(dxdy,1, (gsl_vector_get(y,0)*gsl_vector_get(y,1))/(2*args->N*args->Tc) + gsl_vector_get(y,1)*(gsl_vector_get(y,0)/(2*args->N*args->Tc) - 1/args->Tr));
	gsl_vector_set(dxdy,2, gsl_vector_get(y,1)/args->Tr);

}

int main() {
	struct Args args;
	args.N = 5800000.;
	args.Tc = 5.;
	args.Tr = 15.;
	int n = 2;
	double a = 0.;
	double b = 20.;
	double acc = 0.0005;
	double eps = 0.00005;
	gsl_vector *ya = gsl_vector_alloc(n);
	gsl_vector *yb = gsl_vector_alloc(n);

	gsl_vector_set(ya, 0, 0);
	gsl_vector_set(ya, 1, 1);

	FILE* data_stream = fopen("data_A3.txt","w");
	driver(data_stream, &args, &harm_osc, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_free(ya);
	gsl_vector_free(yb);

	n = 3;
	b = 100.;
	acc = 0.0001;
	eps = 0.000005;
	ya = gsl_vector_alloc(n);
	yb = gsl_vector_alloc(n);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));

	data_stream = fopen("data_SIR1.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));
	args.Tc = 4.;

	data_stream = fopen("data_SIR2.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));
	args.Tc = 3.;

	data_stream = fopen("data_SIR3.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));
	args.Tc = 2;

	data_stream = fopen("data_SIR4.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));
	args.Tc = 1;

	data_stream = fopen("data_SIR5.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);

	gsl_vector_set_zero(ya);
	gsl_vector_set(ya, 1, 2000);
	gsl_vector_set(ya, 2, 0);
	gsl_vector_set(ya, 0, args.N-gsl_vector_get(ya, 0)-gsl_vector_get(ya, 1));
	args.Tc = 0.5;

	data_stream = fopen("data_SIR6.txt","w");
	driver(data_stream, &args, &SIR, a, ya, b, yb, acc, eps);
	fclose(data_stream);
  return 0;
}
