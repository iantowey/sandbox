#include <mm.h>

/*
maxima check

A : matrix ([0,1,2],[1,2,3],[2,3,4],[3,4,5],[4,5,6]);
B : matrix([0,-1,-2,-3],[1,0,-1,-2],[2,1,0,-1]);
C = A.B;

[ 5   2  - 1   - 4  ]
[ 8   2  - 4   - 10 ]
[ 11  2  - 7   - 16 ]
[ 14  2  - 10  - 22 ]
[ 17  2  - 13  - 28 ]

#check no memory leaks
valgrind --leak-check=full ./mm_app

*/


int main(){
    	mm_main();
	return 0;
}

