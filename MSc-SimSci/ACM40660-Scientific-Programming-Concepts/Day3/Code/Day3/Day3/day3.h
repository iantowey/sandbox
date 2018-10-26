#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#define SIZE 256

#define local_int(i, val) int i = val
#define FOR(i, i_init, n) for(i = i_init; i < n; i++)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (a > b ? a : b)

typedef struct{
	struct cell *next;
	int val;
	void (* push)();
	int (* pop)();
	int (* peek)();
} cell; //queue

cell* newcell(int cell_val);

typedef struct {
	char* label;
	int size;
	int *arr;
	void (* push)();
	int (* pop)();
	int (* peek)();
} int_stack;


void push(int_stack*, int) ;
int pop(int_stack*) ;
int peek(int_stack*) ;

typedef  int_stack tower;
void move(int size, tower*, tower*, tower*);
void print_towers(int max_size, tower *t1, tower *t2, tower *t3);

int* allocatearray(int size);
void fillwithones(int *arr, int size);
void printarray(int *arr, int size);

void tower_of_hanoi();
double* exp_power_series_n_terms(double  x, int n);
int gcd(int a , int b);
long factorial_nonrecursive(int n);
long factorial_recursive(int n);
void time_eg();

//utils
double sum_array(double* arr, int length);

//Tests
void factorial_test();
void gcd_test();
void exp_test();
void tower_of_hanoi_test();
void mem_alloc_tests();
void cell_tests();