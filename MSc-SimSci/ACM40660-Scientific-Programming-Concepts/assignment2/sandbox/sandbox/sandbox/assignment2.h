#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define MEM_NULL_CHECK(file_name, line_num, function_name, ptr, error_message) if(ptr == NULL) {printf("%s::%d::%s\t%s\n", file_name, line_num, function_name, error_message); exit(-1);}
#define STR_TO_INT_ERROR_CHK(file_name, line_num, function_name, val, error_message) if(val <= 0) {printf("%s::%d::%s\t%s\n", file_name, line_num, function_name, error_message); exit(-1);}

static int CYCLE_FREQUENCY = 1000;

int* randomarray(int n, int max);
void bubblesort(int n, int *arr);
void bubblesort_input();
void print_array(int n, int *arr);
int search(int i, int n, int *arr);
int chopsearch(int i, int n, int *arr, int amin, int amax);
int chopsearch_default(int i, int n, int *arr);
void benchmark_native(int n, int max, int s, int mult);
void benchmark_chop(int n, int max, int s, int mult);