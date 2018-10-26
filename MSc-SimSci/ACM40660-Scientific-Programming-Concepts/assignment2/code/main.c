#include <assignment2_tests.h>


int main(int argc, char** argv){

	int test_idx = 0;

	srand(time(NULL));   // init RNG

	if(argc == 1){
		printf("Usage:\n");
		printf("\t1)randomarray_tests\t\n");
		printf("\t2)bubblesort_tests\t\n");
		printf("\t3)bubblesort_input_tests\t\n");
		printf("\t4)search_tests\t\n");
		printf("\t5)chopsearch_tests\t\n");
		printf("\t6)benchmark_tests\t\n");
		return EXIT_FAILURE;
	}

	test_idx = atoi(argv[1]);

	STR_TO_INT_ERROR_CHK(__FILE__, __LINE__, __FUNCTION__, test_idx, "Input error, 'text_idx' must be a positive integer in range [1-6]");

	switch(test_idx){
		case 1:
			randomarray_tests();
			break;
		case 2:
			bubblesort_tests();
			break;
		case 3:
			bubblesort_input_tests();
			break;
		case 4:
			search_tests();
			break;
		case 5:
			chopsearch_tests();
			break;
		case 6:
			benchmark_native_test();
			benchmark_chop_test();
			break;
		default:
			printf("No test found at index %d, please retry entering an integer value in range [0-6]", test_idx);
		}
	return EXIT_SUCCESS;
}
