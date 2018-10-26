#include <collatz.h>

/**
 *	Sequence Function definition
 */
int collate(int n){
    return (n % 2 == 0 ? n / 2 : 3 * n +1);
}

/**
 *	Create sequence using dynamic memory allocation, return pointer to new memory location
 */
sequence *create_sequence(int size){
    sequence *s_ptr = malloc(sizeof(sequence));
    assert( s_ptr != NULL);
    s_ptr->len = size;
    s_ptr->seq = (int*)malloc(size * sizeof(int));
    assert( s_ptr->seq != NULL);
    return s_ptr;
}

/**
 *	Increment sequence array
 */
void increment_sequence(sequence *s, int val){
    s->len++;
    s->seq = (int*)realloc(s->seq,s->len * sizeof(int));
    s->seq[s->len-1] = val;
}

/**
 *	Free all memory allocated to the sequence
 */
void free_sequence(sequence *s){
    free(s->seq);
    free(s);
}

/**
 *	Print the sequence no allowing more than 40 characters on a line
 */
void print_sequence(sequence *s){

    int i = 0;
    int print_buffer_len, tmp_buffer_len;
    char print_buffer[MAX_PRINT_STRING_LENGTH], tmp_buffer[1000];

    memset(print_buffer,0,sizeof(print_buffer));
    memset(tmp_buffer,0,sizeof(tmp_buffer));

    for(i = 0; i < s->len; i++){
        sprintf(tmp_buffer,"%d, ", s->seq[i]);
        tmp_buffer_len = strlen(tmp_buffer);
        print_buffer_len = strlen(print_buffer);

        if (tmp_buffer_len + print_buffer_len <= MAX_PRINT_STRING_LENGTH){
            strcat(print_buffer, tmp_buffer);
        }else{
            printf("%s\n",print_buffer);
            memset(print_buffer,0,sizeof(print_buffer));
            strcat(print_buffer, tmp_buffer);
        }

        memset(tmp_buffer,0,sizeof(tmp_buffer));
    }
    printf("%s\n",print_buffer);

}


/**
 *	Run collatz sequence procedure
 */
void collatz_main(int argc, char** argv)
{
    int val = 0;
    sequence *s = create_sequence(1);

    if (argc < 2){
        printf("Error: Missing Parameter 1 : Enter sequence seed (postive integer)\n");
        exit(EXIT_FAILURE);
    }

    s->seq[s->len-1] = atoi(argv[1]);
    assert(s->seq[s->len-1] > 0);

    while(TRUE){
        val = collate(s->seq[s->len-1]);
        increment_sequence(s, val);
        if(s->len >= 3){
            if(s->seq[s->len-1] == 1 && s->seq[s->len-2] == 2 && s->seq[s->len-3] == 4)
            break;
        }
    }
    print_sequence(s);
 
    free_sequence(s);
}


