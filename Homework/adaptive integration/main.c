#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"
#include <assert.h>

double adapt24 (double f (double ), double a, double b,
 double acc, double eps, double f2, double f3, int nrec){
	assert(nrec<1000000);
	double f1 = f(a+(b-a)/6), f4 = f(a+5*(b-a)/6);
	double Q=(2*f1 + f2 + f3 + 2*f4)/6*(b - a) , q = (f1 + f4 + f2 + f3)/4*(b - a);
	double tolerance = acc + eps*fabs(Q), error=fabs(Q-q) ;
	if (error < tolerance) return Q;
	else {
		double Q1=adapt24(f, a, (a+b)/2, acc/sqrt(2.), eps, f1, f2, nrec+1);
		double Q2=adapt24(f, (a+b)/2, b, acc/sqrt(2.), eps, f3, f4, nrec+1);
		return Q1+Q2; }
}
double adapt (double f(double), double a, double b, double acc, double eps){
	double f2=f(a+2*(b-a)/6), f3=f(a+4*(b-a)/6);
	int nrec=0;
	return adapt24(f ,a ,b ,acc ,eps ,f2 ,f3, nrec);
}



int main() { // uses gcc nested functions
	int calls=0; double a=0, b=1, acc=0.0005, eps = 0.00001;
	double Q;
	double f(double x){
		calls++;
		return 1/sqrt(x);
	}; // nested function

	double fun1(double x){
		calls++;
		return sqrt(x)-2./3.;
	};

	double fun2(double x){
		calls++;
		return 4*sqrt(1-x*x)-M_PI;
	};

	Q = adapt(f, a ,b ,acc ,eps);
	printf("f(x)=1/sqrt(x)  Q = %g, calls = %d\n", Q, calls);
	calls = 0;

	Q = adapt(fun1, a ,b ,acc ,eps);
	printf("f(x)=sqrt(x)-2/3  Q = %g, calls = %d\n", Q, calls);
	calls = 0;

	Q = adapt(fun2, a ,b ,acc ,eps);
	printf("f(x)=4*sqrt(x)-pi  Q = %g, calls = %d\n", Q, calls);
	calls = 0;
	return 0; }
