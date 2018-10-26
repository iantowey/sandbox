#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define TRUE 1
#define MAX_PRINT_STRING_LENGTH 40

typedef struct{
    int len;
    int* seq;
} sequence;


int collate(int);

sequence *create_sequence(int);

void increment_sequence(sequence *, int );

void print_sequence(sequence *);

void free_sequence(sequence *);

void collatz_main(int, char**);

