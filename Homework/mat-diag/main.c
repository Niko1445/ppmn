#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "vec_calc.h"

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_Apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_Aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_Apj);
		gsl_matrix_set(A,q,j,new_Aqj);
		}
}
void jacobi_diag(gsl_matrix *A, gsl_matrix *V){
  int changed;
  do{
  	changed=0;
  	for(int p = 0; p < A->size1-1; p++)
  	for(int q = p+1; q < A->size1; q++){
  		double Apq=gsl_matrix_get(A,p,q);
  		double App=gsl_matrix_get(A,p,p);
  		double Aqq=gsl_matrix_get(A,q,q);
  		double theta=0.5*atan2(2*Apq,Aqq-App);
  		double c=cos(theta),s=sin(theta);
  		double new_App=c*c*App-2*s*c*Apq+s*s*Aqq;
  		double new_Aqq=s*s*App+2*s*c*Apq+c*c*Aqq;
  		if(new_App!=App || new_Aqq!=Aqq) // do rotation
  			{
  			changed=1;
  			timesJ(A,p,q, theta);
  			Jtimes(A,p,q,-theta); // A←J^T*A*J
  			timesJ(V,p,q, theta); // V←V*J
  			}
  	}
  }while(changed!=0);
}

double wave(int j, int i, int n){
    double x;
    x = sin(((double)j+1.0) * M_PI * ((double)((i + 1.0) / (n + 1))));
    return x;
  }

int main() {

  int n = 4;
  gsl_matrix *A = gsl_matrix_alloc(n, n);
  gsl_matrix *A_copy = gsl_matrix_alloc(n, n);

  gsl_matrix *V = gsl_matrix_alloc(n, n);
  gsl_matrix *B = gsl_matrix_alloc(n, n);
  gsl_matrix *D = gsl_matrix_alloc(n, n);

  gsl_matrix_set_identity(V);

  int S[4][4] = {4,-30,60,-35,-30,300,-675,420,60,-675,1620,-1050,-35,420,-1050,700};

  for (size_t x = 0; x < n; x++) {
    for (size_t y = 0; y < n; y++) {
      gsl_matrix_set(A,x,y,S[x][y]);
    }
  }

  printf("Matrix A\n");
  matrix_print(stdout,A);

  gsl_matrix_memcpy(A_copy, A);


  jacobi_diag(A, V);

  printf("Matrix V from jacobi_diag\n");
  matrix_print(stdout,V);


  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, A_copy, 0, B);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, V, 0, D);

  printf("Matrix D from V^TAV\n");
  matrix_print(stdout,D);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, D, 0, B);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, B, V, 0, A);

  printf("Matrix A from V^TDV\n");
  matrix_print(stdout,A);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, B);

  printf("Matrix I from V^TV\n");
  matrix_print(stdout,B);

  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_matrix_free(B);
  gsl_matrix_free(V);
  gsl_matrix_free(D);

  // Part B

  n=60;
  double s=1.0/(n+1);
  gsl_matrix* H = gsl_matrix_alloc(n,n);
  V = gsl_matrix_alloc(n,n);
  gsl_matrix_set_identity(V);

  for(int i=0;i<n-1;i++){
    gsl_matrix_set(H,i,i,-2);
    gsl_matrix_set(H,i,i+1,1);
    gsl_matrix_set(H,i+1,i,1);
    }
  gsl_matrix_set(H,n-1,n-1,-2);
  gsl_matrix_scale(H,-1/s/s);

  jacobi_diag(H, V);

  printf("Energies: calculated, exact\n" );
  for (int k=0; k < n/3; k++){
    double exact = M_PI*M_PI*(k+1)*(k+1);
    double calculated = gsl_matrix_get(H,k,k);
    printf("%i %g %g\n",k,calculated,exact);
  }
  gsl_vector *a = gsl_vector_alloc(n);
  double scale;
  for (int x = 0; x < 3; x++) {
    gsl_matrix_get_col(a,V,x);
    scale = gsl_vector_max(a);
    for (int y = 0; y < n; y++) {
      gsl_matrix_set(V,y,x, gsl_matrix_get(V,y,x)/scale);
    }
  }

  FILE *data_stream = fopen("data_B.txt","w");
  fprintf(data_stream,"%g  %g  %g  %g  %g  %g  %g\n", 0.,0.,0.,0.,0.,0.,0.);
  for(int i=0;i<n;i++){
	   fprintf(data_stream,"%g  %g  %g  %g  %g  %g  %g\n",(i+1.0)/(n+1),
      gsl_matrix_get(V,i,0), gsl_matrix_get(V,i,1), gsl_matrix_get(V,i,2),
      wave(0,i,n), -wave(1,i,n), wave(2,i,n));
  }
  fprintf(data_stream,"%g  %g  %g  %g  %g  %g  %g\n", 1.,0.,0.,0.,0.,0.,0.);
  fclose(data_stream);
  return 0;
}
