#include"komplex.h"
#include"stdio.h"
#define TINY 1e-6

int main(){
	komplex a = {1,2}, b = {3,4};
	komplex_print("a =",a);
	komplex_print("b =",b);
	komplex_set(&a,4,6);
	komplex_print("a_set =",a);
	a = komplex_sub(a,b);
  komplex_print("a_sub =",a);
	komplex c = komplex_add(a,b);
	komplex_print("c =",c);
  return 0;
}
