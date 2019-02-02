#include "mylib.h"
#include <stdio.h>

int main(){

	mpz_t n;
	mpz_init(n);
	
	char str[500];
	fgets(str, 500, stdin);

	codifica(n, str);	
	gmp_printf("str");
	strcpy(str,decodifica(n));

	printf("%s \n", str);
	return 0;
}
