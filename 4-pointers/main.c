#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

void set0(double x){
	x=0;
}

void set0p(double *x){ (*x)=0; }

void print_array(int n, double a[]){
	for(int i=0;i<n;i++)printf("print_array: a[%d]=%g\n",i,a[i]);
}

int main(){
	double y=1;
	set0(y);
	printf("y after set0 =%g\n",y);
	set0p(&y);
	printf("y after set0p =%g\n",y);
	int n=5;
	double v[n];
	for(int i=0;i<n;i++){	// i++ : i=i+1
		v[i]=i;
		}
	int i=0; while(i<n) {
		printf("v[%d]=%g\n", i, v[i]);
		i++;
	}
	print_array(n,v);
	int N=7;
	double *a = malloc(N*sizeof(double));
	for(int i=0;i<N;i++) a[i]=i+100;
	print_array(N,a);
	
free(a);
return 0;
}
