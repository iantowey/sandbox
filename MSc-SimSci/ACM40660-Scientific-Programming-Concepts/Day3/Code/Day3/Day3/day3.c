#include <day3.h>

cell* newcell(int cell_val){
	cell* new_cell = (cell*)malloc(sizeof(cell));
	new_cell->next = NULL;
	new_cell->val = cell_val;
	return new_cell;
}


/*void push(forward_linked_list *cur_link,int nextval) {
	forward_linked_list next_link = (forward_linked_list*)malloc(sizeof(forward_linked_list));
	next_link->next = NULL;
	next_link->val = nextval;
	cur_link->next = next_link;
}

int pop(forward_linked_list *cur_link) {
	cur_link->
	int val= t->arr[t->size-1];
	t->size--;
	if(t->size == 0){
		t->arr = NULL;
	} else {
		t->arr = (int*)realloc(t->arr, t->size * sizeof(int));		
	}
	return val;
}*/

void push(tower *t,int val) {
	t->size++;
	if(t->size == 1){
		t->arr = (int*)malloc(sizeof(int));
	} else {
		t->arr = (int*)realloc(t->arr, t->size*sizeof(int));		
	}
	t->arr[t->size-1] = val;
}

int peek(tower *t) {
	if(t->size == 0 ){
		return INT_MAX;
	}
	return t->arr[t->size-1];
}

int pop(tower *t) {
	if (t->size == 0){
		return INT_MAX;
	}
	int val= t->arr[t->size-1];
	t->size--;
	if(t->size == 0){
		t->arr = NULL;
	} else {
		t->arr = (int*)realloc(t->arr, t->size * sizeof(int));		
	}
	return val;
}

int* allocatearray(int size){
	return (int*)malloc(size*sizeof(int));
}

void fillwithones(int *arr, int size){
	int i = 0;
	for(i = 0; i < size;i++){
		arr[i] = 1;
	}
}
void printarray(int *arr, int size){
	int i = 0;
	for(i = 0; i < size;i++){
		printf("%d\t", arr[i]);
	}
}


void print_towers(int max_size, tower *t1, tower *t2, tower *t3){
	int i = 0, j = 0;
	char buf[1024];
	printf("----------------------------------------------------------------------\n");
	printf("%s\t\t\t\t%s\t\t\t\t%s\n",t1->label,t2->label,t3->label);
	printf("----------------------------------------------------------------------\n");
	for(i = max_size-1; i >= 0; i--){
		if (t1->size-1 >= i){
			for (int j = 0;  j < t1->arr[i];  j++, printf("%c", '*'));
		}
		printf("\t\t\t\t");
		if (t2->size-1 >= i){
			for (int j = 0;  j < t2->arr[i];  j++, printf("%c", '*'));
		}
		printf("\t\t\t\t");
		if (t3->size-1 >= i){
			for (int j = 0;  j < t3->arr[i];  j++, printf("%c", '*'));
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------------\n");
}

void print_tower(tower *t){
	int a = 0, b = 0;
	printf("__________________________\n");
	printf("%s\n", t->label);
	printf("__________________________\n");
	for(a = t->size-1; a >= 0; a--){
		for(b = 0; b < t->arr[a] ; b++){
			printf("*");
		}
		printf("\n");
	}
	printf("__________________________\n");
}

tower *create_tower(char *label,int size){
	int i = 0;
	tower *t = (tower*)malloc(sizeof(tower));
	t->label = label;
	t->size = size;
	if(t->size > 0){
		int *arr = (int*)malloc(t->size * sizeof(int));			
		for(i = 0; i < t->size; i++){
			arr[i] = t->size - i;
		}
		t->arr = arr;
	} else {
		t->arr = NULL;
	}
	t->peek = peek;
	t->push = push;
	t->pop = pop;
	return t;
}

void move(int size, tower* src, tower* tgt, tower* inter){
	if (size > 0){
		move(size-1, src, inter, tgt);
		printf("%s ---%d--->%s\n",src->label, src->peek(src), tgt->label);
		tgt->push(tgt, src->pop(src));
		move(size-1, inter, tgt, src);
	}
}


void tower_of_hanoi(int num_disks){
	//init
	int iter_break = 0;
	tower *A = create_tower("src", num_disks);
	tower *B = create_tower("inter", 0);
	tower *C = create_tower("tgt", 0);

	tower *src = A;
	tower *inter = B;
	tower *tgt = C;

	print_towers(num_disks, A,B,C);

	move(num_disks, src, tgt, inter);

	print_towers(num_disks, A,B,C);

	free(A->arr);
	free(B->arr);
	free(C->arr);

	free(A);
	free(B);
	free(C);

}

int gcd(int a , int b){
	if (b == 0){
		return a;
	} else {
		return gcd(b, a % b);
	}
}

long factorial_nonrecursive(int n){

	long product = 1;
	int k = n;
	while(k > 0){
		product *= (k--);
	}
	return product;
}

long factorial_recursive(int n){

	if(n == 0){
		return 1;
	} else {
		return n * factorial_recursive(n-1);
	}
}

void time_eg(){
  char buffer[SIZE];
  time_t curtime;
  struct tm *loctime;

  /* Get the current time. */
  curtime = time (NULL);

  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  /* Print out the date and time in the standard format. */
  fputs (asctime (loctime), stdout);

  /* Print it out in a nice format. */
  //
  //	http://www.cplusplus.com/reference/ctime/strftime/
  strftime (buffer, SIZE, "Today is %A, %B %d.\n", loctime);
  fputs (buffer, stdout);
  strftime (buffer, SIZE, "The time is %I:%M %p.\n", loctime);
  fputs (buffer, stdout);

}

double* exp_power_series_n_terms(double  x, int n){
	local_int(i, 0);
	double* terms = (double *)malloc(n*sizeof(double));

	FOR(i, 0, n){
		terms[i] = pow(x,i)/factorial_nonrecursive(i);
	}
	return terms;
}



/*
	******* ******* ******* ******* *******
	   *    *		*		   *	*
	   *	******* *******	   *	*******
	   *	*             *    *		  *
	   *	******* *******    *	*******
*/


void factorial_test(){
	printf("n\t\t * factorial_nonrecursive\t * factorial_recursive\n");
	printf("**********************************************************************\n");

	local_int(i, 0);
	FOR(i, 0, 5){
		printf("%d\t\t * %ld\t\t\t\t * %ld\n",i, factorial_nonrecursive(i), factorial_recursive(i));
	}
	printf("**********************************************************************\n");
}

void gcd_test(){
	local_int(i, 0);
	local_int(j, 0);
	printf("i\t\t * j\t * gcd\n");
	printf("**********************************************************************\n");
	FOR(i, 1, 20){
		FOR(j, 1, 20){
			if (i <= j) continue;
			printf("%d\t\t * %d\t\t\t\t * %d\n",i, j, gcd(i, j));
		}
	}
	printf("**********************************************************************\n");
}

double sum_array(double* arr, int length){
	double sum = 0.0;
	local_int(i, 0);
	FOR(i,0,length){
		sum += arr[i];
	}
	return sum;
}

void exp_test(){
	printf("n\t\t * exp estimate\t\t\t\t * estimate error\n");
	printf("**********************************************************************\n");

	double tmp = 0;
	local_int(i, 0);
	local_int(N, 30);
	double *terms = exp_power_series_n_terms(1.0, N);
	FOR(i, 1, N){\
		tmp = sum_array(terms, i);
		printf("%d\t\t * %.20f\t\t * %.30f\n",i, tmp, fabs(tmp-exp(1)));
	}
	printf("**********************************************************************\n");
	free(terms);
}

void tower_of_hanoi_test(){
	tower_of_hanoi(10);
}

void mem_alloc_tests(){
	int n = 10;
	int* arr = allocatearray(n);
	fillwithones(arr, n);
	printarray(arr, n);
	printf("\n");
	free(arr);
	n = 5;
	arr = allocatearray(n);
	fillwithones(arr, n);
	printarray(arr, n);
	free(arr);
}

void print_cells(cell* cell){
	printf("%d\t", cell->val);
	while(cell->next != NULL){
		print_cells(cell->next);
	}
}

void cell_tests(){
	cell *c1  = newcell(15);
	printf("%d\n", c1->val);
//	cell *c2  = newcell(151);
//	c1->next = c2;
//	print_cells(c1);
	free(c1);
//	free(c2);
}