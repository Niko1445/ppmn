#include<stdio.h>
#include<pthread.h>
#include<math.h>

//x->y  =  (*x).y

void* bar(void* arg){
	double* x= (double*) arg;
	for(int i=0;i<1e8;i++) *x=cos(*x);

return NULL;
}

int main(){
	double x=0,y=100,z=-100;
	pthread_t threadx, thready; //threadz
	pthread_attr_t* attributes = NULL;
	pthread_create( &threadx, attributes, bar, (void*)&x );
	pthread_create( &thready, attributes, bar, (void*)&y );
//	pthread_create( &threadz, attributes, bar, (void*)&z );
//	bar((void*)&x);
//	bar((void*)&y);
	bar((void*)&z);
	void* returnvalue = NULL;
	pthread_join(threadx,returnvalue);
	pthread_join(thready,returnvalue);
	printf("x=%g cos(x)=%g\n",x,cos(x));
	printf("y=%g cos(y)=%g\n",y,cos(y));
	printf("z=%g cos(z)=%g\n",z,cos(z));
return 0;
}


