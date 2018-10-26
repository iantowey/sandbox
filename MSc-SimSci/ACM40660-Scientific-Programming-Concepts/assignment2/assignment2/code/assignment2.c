#include <assignment2.h>

/*! generic search function  benchmark procedure
	param n 			integer, array
	param max 			integer, max possible random value in array
	param s 			integer, value to find in array
	param mult			integer, number of iterations
	param search_func	function_ptr, function that will be benchmarked
*/
void benchmark(int n, int max, int s, int mult, int (*search_func)(int, int, int*)){
	int i;
	double clk_per_sec;
	int *arr = NULL;
	clock_t start, end;

	printf("\n---------------------------------------------------------------------------------------------------\n");
	printf("\n\tParameters:\n");
	printf("\t\t'n'(length of array) \t\t\t\t:\t%d\n", n);
	printf("\t\t'max'(max random value) \t\t\t:\t%d\n", max);
	printf("\t\t's'(value to search array for) \t\t\t:\t%d\n", s);
	printf("\t\t'mult'(number of iterations of sort to compute) :\t%d\n\n", mult);
	printf("\n---------------------------------------------------------------------------------------------------\n");

	start = clock ();

	for(i = 0; i < mult; i++){
		if(i % CYCLE_FREQUENCY == 0){
			if(arr != NULL) free(arr);
			arr = randomarray(n, max);
			bubblesort(n, arr);
		}
		search_func(s, n, arr);
	}
	free(arr);
	end = clock ();

	clk_per_sec = end - start;
	printf("total time taken \t: %10.3f (CLOCKS_PER_SEC)\n", clk_per_sec);
	printf("avg time \t\t: %10.9f (CLOCKS_PER_SEC)\n", clk_per_sec / mult);
	printf("---------------------------------------------------------------------------------------------------\n");
}


/*! benchmark function 'chopsearch' defaulting to  amin = 0 and amax = n-1.
	param n 	integer, array
	param max 	integer, max possible random value in array
	param s 	integer, value to find in array
	param mult	integer, number of iterations
*/
void benchmark_chop(int n, int max, int s, int mult){
	benchmark(n, max, s, mult, chopsearch_default);
}

/*! benchmark function 'search'.
	param n 	integer, array
	param max 	integer, max possible random value in array
	param s 	integer, value to find in array
	param mult	integer, number of iterations
*/
void benchmark_native(int n, int max, int s, int mult){
	benchmark(n, max, s, mult, search);

}

/*! recursive chopsearch function which returns the index of the value being searched for without offsetting against amin.
	param i 	integer, value to find in array
	param n 	integer, array
	param arr	pointer to integer array
	return index of value i in arr, or -1 if not found  
*/
int chopsearch_inner(int i, int n, int *arr, int amin, int amax){
	int idx = -1;
	int mid = floor((amax+amin)/2);
	if(amin > amax){
		return -1;
	}
	if(i < arr[mid]){
		idx = chopsearch_inner(i, n, arr, amin, mid-1);
	} else if(i > arr[mid]){
		idx = chopsearch_inner(i, n, arr, mid+1, amax);
	} else {
		//i == arr[mid]
		while(i == arr[--mid] && mid >= amin){}
		idx = mid + 1;
	}

	return idx;

}

/*! chopsearch defaulting amin to 0 and amax = n-1, should behave the same as the function search.
	param i 	integer, value to find in array
	param n 	integer, array
	param arr	pointer to integer array
	return index of value i in arr, or -1 if not found  
*/
int chopsearch_default(int i, int n, int *arr){
	return chopsearch(i, n, arr,0, n-1);
}


/*! search subset of Integer array for value i, returning indx from that subset where value found 
	param i 	integer, value to find in array
	param n 	integer, array
	param arr	pointer to integer array
	param amin	inter, subset array min index
	param amax	inter, subset array min index
	return index of value i in arr, or -1 if not found  
*/
int chopsearch(int i, int n, int *arr, int amin, int amax){
	int idx = chopsearch_inner(i, n, arr, amin, amax);
	return idx < 0 ? -1 : idx - amin;
}


/*! search Integer array for value i, returning indx where value found 
	param i 	integer, value to find in array
	param n 	integer, array
	param arr	pointer to integer array
	return index of value i in arr, or -1 if not found  
*/
int search(int i, int n, int *arr){
	int j, idx = -1;
	if (i >= arr[0]){
		for(j = 0; j < n; j++){
			if(arr[j] == i){
				idx = j;
				break;
			}
			if(arr[j] > i){
				break;
			}
		}
	}
	return idx;
}

void print_array(int n, int *arr){
	int i;
	printf("\n***************************************************\n");
	for(i = 0; i < n; i++){
		printf("[%d] = %d\t", i, arr[i]);
	}
	printf("\n***************************************************\n");
}

/*! Return a array of length 'n', populated with random integer values from range [0, max] 
	param n 	integer, array length
	param max 	integer, max possible random value in array
	return pointer to array 
*/
int* randomarray(int n, int max){
	int i = 0;
	int *arr = (int*)malloc(n*sizeof(int));
	MEM_NULL_CHECK(__FILE__, __LINE__, __FUNCTION__, arr, "Error :: allocating memory, aborting!");
	for(i = 0; i < n ; i++){
		arr[i] = rand() % (max+1);
	} 
	return arr;
}

/*! Recursive bubble sort  
	param n 	integer, array length
	param *arr  pointer to integer array
*/
void bubblesort(int n, int *arr){
	int i, tmp;
	if(n == 1) return;
	for(i = 0; i < n-1; i++){
		if(arr[i] > arr[i+1]){
			tmp = arr[i];
			arr[i] = arr[i+1];
			arr[i+1] = tmp;
		}
	}
	bubblesort(n-1, arr);
}

/*! Get user input to create an array of length 'n' (provided by user) and max random value 'max' (provided by user)  
*/
void bubblesort_input(){
	int n, max;
	char buf[1024];
	int *arr;
	
	printf("Generate random array of length n and max random val\n");  //displays string  
	printf("Enter array length (Integer): ");    
	scanf("%s", buf); //reads value of n from user
	n = atoi(buf);
	STR_TO_INT_ERROR_CHK(__FILE__, __LINE__, __FUNCTION__, n, "Input error, 'n' must be a positive integer");
	printf("Enter max random value size (Integer): ");    
	scanf("%s", buf); //reads value of n from user
	max = atoi(buf);
	STR_TO_INT_ERROR_CHK(__FILE__, __LINE__, __FUNCTION__, max, "Input error, 'max' must be a positive integer");
	printf("Generating random array of length n = %d and max = %d random values\n", n, max);  //displays string  
	arr = randomarray(n, max);
	printf("\nBefore sorting\n");  //displays string  
	print_array(n, arr);
	bubblesort(n, arr);
	printf("After sorting\n");  //displays string  
	print_array(n, arr);
	free(arr);
}
