#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

/*
struct vector {int n, double a[]};
struct vector my_vector;
typedef struct vector vector;
*/
typedef struct {int n; double* data;} vector;
vector* vector_alloc(int n){
	vector* v=malloc(sizeof(vector));
	(*v).n=n;
	(*v).data=malloc(n*sizeof(double));
	return v;
	}

void vector_free(vector* v){
	free((*v).data);
	free(v);
	}

void vector_set(vector* v,int i, double value){
	assert(i>=0);
	assert(i<(*v).n);
	(*v).data[i]=value;
	}

double vector_get(vector* v,int i){
	assert(i>=0);
	assert(i<(*v).n);
	return (*v).data[i];
}

void vector_print(vector* v){
	for(int i=0;i<(*v).n;i++){
		double vi=vector_get(v,i);
		printf("print_vector: vector[%d]=%g\n",i,vi);
		}
	}


int main(){
	int n=7;
	vector* v=vector_alloc(n);
	for(int i=0;i<n;i++)vector_set(v,i,i+100);
	vector_print(v);

	
vector_free(v);
return 0;
}
