#include<stdio.h>

#define PI 3.1415927

void world(void);
double x;
char hello[] = "file scope";

void bar(void){
	char hello[] = "function scope";
	printf("hello (function scope?) :%s\n",hello);
	{
		char hello[] = "block scope";
		printf("hello (block scope?) :%s\n", hello);
	}
}
//comment line


/*
block comment
*/

int main(){
	printf("hello\n");
	world();
	printf("hello (file scope?) :%s\n",hello);
	bar();
	x=7;
	printf("x=%g\n",x);
return 0;
}

