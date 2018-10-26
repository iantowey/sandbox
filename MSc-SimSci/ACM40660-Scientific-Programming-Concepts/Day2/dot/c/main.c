#include <dot.h>

/*

maxima check

u : matrix ([1,2,3,4,5]);
v : matrix ([1,2,3,4,5]);
C = u.v;

#check no memory leaks
valgrind --leak-check=full ./dot_app


*/


int main(){
    	dot_main();
	return 0;
}

