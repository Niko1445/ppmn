#include<stdio.h>
#include<math.h>

int main(){
	double x;
	int items;
	FILE* my_out_stream = fopen("out.file.txt","w");
	do{
		items=fscanf(stdin,"%lg",&x);
		printf("x=%g sin(x)=%g\n",x,sin(x));
		fprintf(stderr,"x=%g sin(x)=%g\n",x,sin(x));
		fprintf(my_out_stream,"x=%g sin(x)=%g\n",x,sin(x));
	}while(items!=EOF);
fclose(my_out_stream);
return 0;
}
