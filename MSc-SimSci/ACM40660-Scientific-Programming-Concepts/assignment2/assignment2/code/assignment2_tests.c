	#include <assignment2_tests.h>

void benchmark_chop_test(){
	printf("\n***************************************************************************************************\n");
	printf("* chopsearch benchmark test\n");
	printf("***************************************************************************************************\n");
	benchmark_chop(2000, 10000, 10, 1000000);
	printf("-----------------------------------------------------------------------------------------------------\n");
	benchmark_chop(2000, 10000, 5000, 1000000);
	printf("-----------------------------------------------------------------------------------------------------\n");
	benchmark_chop(2000, 10000, 9000, 1000000);
	printf("*****************************************************************************************************\n");

}

void benchmark_native_test(){
	printf("\n***************************************************************************************************\n");
	printf("* search benchmark test\n");
	printf("*****************************************************************************************************\n");
	benchmark_native(2000, 10000, 10, 1000000);
	printf("-----------------------------------------------------------------------------------------------------\n");
	benchmark_native(2000, 10000, 5000, 1000000);
	printf("-----------------------------------------------------------------------------------------------------\n");
	benchmark_native(2000, 10000, 9000, 1000000);
	printf("*****************************************************************************************************\n");
}

void chopsearch_tests(){
	int n = 15;
	int arr[] = {8, 101, 1, 101, 9, 3 , 101, 526, 765, 1800,101,101,101,101,101};
	printf("\nUnsorted\n");  //displays string  
	print_array(n, arr);
	bubblesort(n, arr);
	printf("\nSorted\n");  //displays string  
	print_array(n, arr);
	printf("\n***************************Chopsearch Tests**********************\n");
	printf("%d found at idx = %d in sub array [1:3]\n", 9, chopsearch(9, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [1:3]\n", 3, chopsearch(3, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [1:3]\n", 8, chopsearch(8, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [1:3]\n", 526, chopsearch(526, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [0:14]\n", 526, chopsearch(526, n, arr, 0, n-1));
	printf("%d found at idx = %d in sub array [11:14]\n", 526, chopsearch(526, n, arr, 11, n-1));
	printf("%d found at idx = %d in sub array [0:14]\n", 101, chopsearch(101, n, arr, 0, n-1));
	printf("%d found at idx = %d in sub array [0:3]\n", 101, chopsearch(101, n, arr, 0, 3));
	printf("%d found at idx = %d in sub array [7:14]\n", 101, chopsearch(101, n, arr, 7, n-1));
	printf("%d found at idx = %d in sub array [1:3]\n", 102, chopsearch(102, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [1:3]\n", 1801, chopsearch(1801, n, arr, 1, 3));
	printf("%d found at idx = %d in sub array [7:14]\n", 1800, chopsearch(1800, n, arr, 7, n-1));

}

void search_tests(){
	int n = 15;
	int arr[] = {8, 101, 1, 101, 9, 3 , 101, 526, 765, 1800,101,101,101,101,101};
	print_array(n, arr);
	bubblesort(n, arr);
	print_array(n, arr);
	int val = 8;
	printf("\n***************************Search Tests**********************\n");
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
	val = 99;
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
	val = 526;
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
	val = 18001;
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
	val = -1;
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
	val = 1800;
	printf("\n\nSearching for value %d in arr, found at index %d\n", val, search(val, n, arr));
}

void bubblesort_tests(){
	int n, max;
	char buf[1024];

	n = 10, max = 1000;
	int *arr = randomarray(n, max);
	printf("Generate Array\n");
	printf("n = %d\tmax = %d\n", n, max);
	printf("\nUnsorted Array\n");
	print_array(n, arr);
	bubblesort(n, arr);
	printf("\nSorted Array\n");
	print_array(n, arr);
	free(arr);
}

void bubblesort_input_tests(){
	bubblesort_input();
}


void randomarray_tests(){
	int i, n, max;
	int *arr;

	n = 10, max = 100;
	arr = randomarray(n, max);
	printf("\n***************************************************\n");
	printf("Generate Array\n");
	printf("n = %d\tmax = %d\n", n, max);
	print_array(n, arr);
	free(arr);

	n = 20, max = 100;
	arr = randomarray(n, max);
	printf("\n***************************************************\n");
	printf("Generate Array\n");
	printf("n = %d\tmax = %d\n", n, max);
	print_array(n, arr);
	free(arr);

	n = 50, max = 1000;
	arr = randomarray(n, max);
	printf("\n***************************************************\n");
	printf("Generate Array\n");
	printf("n = %d\tmax = %d\n", n, max);
	print_array(n, arr);
	free(arr);
}