//exe1.c 

#include <stdio.h>
#include<string.h>
int main(void) {
   int account;
   float subtotal;
   float total;

   printf("Enter Account Values :\n");
   scanf("%d", &account);
   printf("Enter Subtotal :\n");
   scanf("%f", &subtotal);
   printf("Enter Total :\n");
   scanf("%f", &total);

   //TODO: Print formatted data to stdout
   printf("----------------------------\n");
   printf("Account %d\tSubtotal: %f\tTotal %f\n",account, subtotal, total);
   printf("----------------------------\n");

  return 0;
}

//run
rm -rf exe1
icc -o exe1 exe1.c 
(echo -en '1\n'; sleep 1; echo -en '1234.56\n'; sleep 1; echo -en '7890.12') | ./exe1

/*
*********************************************************************************************************************
*/

//exe2.c

#include <stdio.h>
#include<string.h>
#include<math.h>

int main(void) {
  	int i, max = 1000;
	float  forward_sum = 0, backward_sum = 0;

	//forward
   	for(i = 0; i < max; i++){
		forward_sum += (float)1/(i+1);
   	}	 

	//forward
   	for(i = max; i > 0; i--){
		backward_sum += (float)1/(i);
   	}	 

   printf("----------------------------\n");
   printf("Forward Sum :: %.30f\n",forward_sum);
   printf("----------------------------\n");

   printf("----------------------------\n");
   printf("Backward Sum :: %.30f\n",backward_sum);
   printf("----------------------------\n");

   printf("Absolute diffDiff :: %.30f\n",fabs(backward_sum-forward_sum));

  return 0;
}

//run
rm -rf exe2 
icc -o exe2 exe2.c -lm
./exe2

/*
*********************************************************************************************************************
*/


//exe3.c

#include <stdio.h>
#include <math.h>

float degree_to_radian(float);

float PI = M_PI;


int main(void) {
   int i, j, dim=13; //Loop index, Counter, Array dimension
   float rad, Tan[dim]; //Return value of the function, Result array

   //TODO: Get table of tan
   j = 0;
   for(i = 0; i < dim; i++){
	Tan[i] = tan(degree_to_radian((float)j));
	//printf("%d\t%f\t%f\n", j, degree_to_radian((float)j), Tan[i]);
	j+=5;
   }

   //TODO: Print the Result array
   j=0;
   for(i = 0; i < dim; i++){
	printf("degrees\t%d => tan\t%f\n",j,Tan[i]);
	j+=5;
   }
   return 0;
}

float degree_to_radian(float degree){
	return PI*degree/180;
}

//run
rm -rf exe3 
icc -o exe3 exe3.c -lm
./exe2

/*
*********************************************************************************************************************
*/


//exe4.c

#include <stdio.h>
#include <math.h>

float degree_to_radian(float);

float PI = M_PI;


int main(void) {


   int i, j, dim=13; //Loop index, Counter, Array dimension
   float rad, Tan[dim], area, coeff; //Return value of the function, Result array, Area, Coefficient

   j = 0;
   for(i = 0; i < dim; i++){
	Tan[i] = tan(degree_to_radian((float)j));
	//printf("%d\t%f\t%f\n", j, degree_to_radian((float)j), Tan[i]);
	j+=5;
   }

   //TODO: Print the Result array
   /*j=0;
   for(i = 0; i < dim; i++){
	printf("degrees\t%d => tan\t%f\n",j,Tan[i]);
	j+=5;
   }*/
   
   //TODO: Calculate the are using Trapezodial rule 
   area = 0;
   printf("i\tcoeff\tTan[i]\tcoeff*Tan[i]\tarea\n");
   for(i = 0; i < dim; i++){
	coeff = (i == 0 || i == (dim-1) ? 1 : 2);
	area += coeff*Tan[i];
	printf("%d\t%f\t%f\t%f\t%f\n",i,coeff, Tan[i], coeff*Tan[i], area);
   }
   area = PI/3*area/(2*dim);

   //TODO: Print the Approximate result and Actual result
   printf("-----------------------------------------------------\n");
   printf("Approximate the value of log(2) using trapezoid rule.\n");
   printf("-----------------------------------------------------\n");
   printf("approx %.20f\t actual %.20f\n", area, log(2));
   printf("-----------------------------------------------------\n");

   return 0;
}

float degree_to_radian(float degree){
	return PI*degree/180;
}

//run
rm -rf exe4 
icc -o exe4 exe4.c -lm
./exe2

/*
*********************************************************************************************************************
*/

