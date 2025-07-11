#include<stdio.h>
#define length 4096
#define RANK_OF_NORM 10
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)
typedef int int_type;

struct data_k{
	int n;
	int iRows_reduced;
	int iCols_reduced;
	int maxRows;
	int_type Lnorm;
	unsigned long long int steps;
	unsigned long long int steps_remainder;
	unsigned long long int copyNum;
	int marginal;
	int_type *mtx_as_vec;
	int *strategy;
};

typedef struct data_k item_calc;

struct data_s{
	int_type **matrix;
	int n_original;
	int iRows;
	int iCols;
	int *original_r;
	int *original_c;
	char fileName[1024];
	char stat;
	int marginal;
	int *original;
	int original_length;
	unsigned long long int total_num_to_calc;
};

typedef struct data_s item;

void matrix_read(item* first){
	int i = 0, j = 0, k = 0;
	int_type value, *row;
	double dValue;
	first->matrix = NULL;
	row = NULL;
	char e='0', g, cNum[256];

	FILE *fp;
	fp = fopen(first->fileName,"r");
	
	do{
        	g = fgetc(fp);
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++; if(g == 'e' || g == 'E') {e = 'e';}}
		else {
            		cNum[i] = '\0'; 
			if(cNum[0] != '\0') {
                		if(e=='e') {sscanf(cNum, "%le", &dValue); value = dValue; e='0'; }
                		else {sscanf(cNum, "%d", &value); }
                		j++;
                		i = 0;
                		row = (int_type*) realloc(row, j * sizeof(int_type));
 	        		if(row != NULL) {row[j-1] = value; }
                		else {perror("Memory allocation failed for row in function matrix_read!\n");}
			}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){
                		first->iCols = j;
                		j = 0;
                		k++;
                		first->matrix = (int_type**) realloc(first->matrix, k * sizeof(int_type*));
                		if(first->matrix != NULL) {first->matrix[k-1] = row; row = NULL;}
                		else {perror("Memory allocation failed for matrix in function matrix_read!\n");}
        		}
		}
	}while(g != EOF);
	first->iRows = k;
	printf("rows: %d, cols: %d\n",first->iRows, first->iCols); 
	fclose(fp);
	if (row != NULL) free(row);
}

void copy_mtx_row(item* first, int* i, int* j){
	int k;
	if(first->original_r[*j] < 0){
		for(k=0; k< first->iCols; k++){
			 first->matrix[*i][k] -= first->matrix[*j][k];
		}
	}
	else {
		for(k=0; k< first->iCols; k++){
			 first->matrix[*i][k] += first->matrix[*j][k];
		}
	}
}

void copy_mtx_col(item* first, int* i, int* j){
	int k;
	if(first->original_c[*j] < 0){
		for(k=0; k< first->iRows; k++){
			 first->matrix[k][*i] -= first->matrix[k][*j];
		}
	}
	else {
		for(k=0; k< first->iRows; k++){
			 first->matrix[k][*i] += first->matrix[k][*j];
		}
	}
}

void eliminate_zero_rows_cols(item* first, item_calc* second){
	short int zero;
	int counter_r, counter_c, i, k;

	first->original_r = (int*) calloc(first->iRows, sizeof(int));
	if(first->original_r == NULL) {perror("Memory allocation failed for original_r in function eliminate_zero_rows_cols!\n");}
	for(i = 0; i < first->iRows; i++){first->original_r[i] = 0;}
	first->original_c = (int*) calloc(first->iCols, sizeof(int));
	if(first->original_c == NULL) {perror("Memory allocation failed for original_c in function eliminate_zero_rows_cols!\n");}
	for(i = 0; i < first->iCols; i++){first->original_c[i] = 0;}
	second->iRows_reduced = first->iRows;
	second->iCols_reduced = first->iCols;
	second->marginal = first->marginal;

	if(first->stat == 'y' || first->stat == 'Y' || first->stat == 'r' || first->stat == 'R' || first->stat == 'c' || first->stat == 'C'){
		counter_r = 0;
		for(i = 0; i < first->iRows; i++){
			if(first->original_r[i] != 0) continue;
			zero = first->matrix[i][0] == 0 ? 1 : 0;
			for(k = 1; (k < first->iCols) && zero; k++){
				if(first->matrix[i][k] != 0) {zero = 0;}
			}
			if(zero) {first->original_r[i] = first->iRows; if(first->marginal == 0 || i > 0) {counter_r++; printf("row %d entries are zeros.\n", i+1);}}
		}
		second->iRows_reduced = first->iRows-counter_r;

		counter_c = 0;
		for(i = 0; i < first->iCols; i++){
			if(first->original_c[i] != 0) continue;
			zero = first->matrix[0][i] == 0 ? 1 : 0;
			for(k = 1; (k < first->iRows) && zero; k++){
				if(first->matrix[k][i] != 0) {zero = 0;}
			}
			if(zero) {first->original_c[i] = first->iCols; if(first->marginal == 0 || i > 0) {counter_c++; printf("col %d entries are zeros.\n", i+1);}}
		}
		second->iCols_reduced = first->iCols-counter_c;
		
		if(first->marginal == 1) {
			if(first->original_r[0] == first->iRows && first->original_c[0] == first->iCols) {second->marginal = 0; second->iRows_reduced--; second->iCols_reduced--; printf("The entries of the first row and column are zeros.\n");}
			else {first->original_r[0] = 0; first->original_c[0] = 0;}
		}
	}
}

void delete_rows(item* first, item_calc* second){
	int *original_r, empty, i, j;
	original_r = (int*) calloc(first->iRows, sizeof(int));
	if(original_r == NULL) {perror("Memory allocation failed for original_r in function delete_rows!\n");}
	for(i=0; i < first->iRows; i++){ original_r[i] = first->original_r[i];}

	empty = 0;
	for(i=0; (i < first->iRows) && (empty < second->iRows_reduced); i++){
		if( (original_r[i] == 0) && (original_r[empty] != 0) ){
			for(j=0; j < first->iCols; j++){
				first->matrix[empty][j] = first->matrix[i][j];
			}				
			original_r[i] = original_r[empty];
			empty++;
		}
		if(original_r[empty] == 0) {empty++;}
	}
	free(original_r);
}

void delete_cols(item* first, item_calc* second){
	int *original_c, empty, i, j;
	original_c = (int*) calloc(first->iCols, sizeof(int));
	if(original_c == NULL) {perror("Memory allocation failed for original_c in function delete_cols!\n");}
	for(i=0; i < first->iCols; i++){original_c[i] = first->original_c[i]; }

	empty = 0;
	for(i=0; (i < first->iCols) && (empty < second->iCols_reduced); i++){
		if( (original_c[i] == 0) && (original_c[empty] != 0) ){
			for(j=0; j < first->iRows; j++){
				first->matrix[j][empty] = first->matrix[j][i];
			}
			original_c[i] = original_c[empty];
			empty++;
		}
		if(original_c[empty] == 0) {empty++;}
	}
	free(original_c);
}

void calc_reduce_matrix_cols_sign(item* first, item_calc* second){
	short int negative, positive;
	int counter_c, i, k, *original, index, mult;
	counter_c = 0;
	original = (int*) calloc(first->iCols, sizeof(int));
	if(original == NULL) {perror("Memory allocation failed for original in function calc_reduce_matrix_cols_sign!\n");}
	for(i = 0; i < first->iCols; i++) {original[i] = 0;}

	counter_c = 0;
	for(i = 0; i < first->iCols; i++){
		if(first->original_c[i] != 0) continue;
		positive = 0;
		negative = 0;
		for(k = 0; (k < first->iRows) && (!positive || !negative); k++){
			if(first->matrix[k][i] > 0) {positive = 1;}
			if(first->matrix[k][i] < 0) {negative = 1;}
		}
		if(positive && negative) {original[i] = 0;}
		else if(!positive && negative) { original[i] = -1; counter_c++; }
		else {original[i] = 1; counter_c++;}
	}

	if(counter_c > 1){
		second->iCols_reduced = second->iCols_reduced-counter_c + 1;
		for(i = 0; original[i] == 0; i++) {}
		index = i;
		for(i = (index + 1); i < first->iCols; i++) {
			if(first->original_c[i] != 0) continue;
			mult = original[index] * original[i];
			if(mult < 0) {
				first->original_c[i] = -(index + 1);
				copy_mtx_col(first, &index, &i);
				printf("col %d was subtracted from %d\n", i+1, index+1);
			}
			if(mult > 0) {
				first->original_c[i] = index + 1;
				copy_mtx_col(first, &index, &i);
				printf("col %d was added to %d\n", i+1, index+1);
			}
		}
	}
	free(original);
}

void calc_reduce_matrix_rows(item* first, item_calc* second){
	short int mult;
	long long int multiple1, multiple2;
	int counter_r, i, j, k, index;
	counter_r = 0;
	
	for(i = first->marginal; i < first->iRows; i++){
		if(first->original_r[i] != 0) continue; 
		for(j = i+1; j < first->iRows; j++){
			if(first->original_r[j] != 0) continue;
			mult = 1;
			for(k = 0; (first->matrix[i][k] == 0) && (first->matrix[j][k] == 0) && (k < first->iCols); k++){}
			if(((first->matrix[i][k] == 0) && (first->matrix[j][k] != 0)) || ((first->matrix[i][k] != 0) && (first->matrix[j][k] == 0))) {mult = 0;}
			
			for(index = k; (k < (first->iCols - 1)) && mult; k++){
				multiple1 = (long long int) first->matrix[i][index] * first->matrix[j][k+1];
				multiple2 = (long long int) first->matrix[j][index] * first->matrix[i][k+1];
				mult = (multiple1 == multiple2) ? 1 : 0;
			}
			if(mult == 0) {first->original_r[j] = 0;}
			else if(mult && (((long long int) first->matrix[i][index] * first->matrix[j][index]) < 0)) {if(first->n_original == 1) {first->original_r[j] = -(i+1); counter_r++; copy_mtx_row(first, &i, &j); printf("row %d was subtracted from %d\n", j+1, i+1);}}
			else {first->original_r[j] = (i+1); counter_r++; copy_mtx_row(first, &i, &j); printf("row %d was added to %d\n", j+1, i+1);}
		}
	}
	second->iRows_reduced = second->iRows_reduced-counter_r;
}

void calc_reduce_matrix_cols(item* first, item_calc* second){
	short int mult;
	long long int multiple1, multiple2;
	int counter_c, i, j, k, index;
	counter_c = 0;
	 
	for(i = first->marginal; i < first->iCols; i++){
		if(first->original_c[i] != 0) continue;
		for(j = i+1; j < first->iCols; j++){
			if(first->original_c[j] != 0) continue;
			mult = 1;
			for(k = 0; (first->matrix[k][i] == 0) && (first->matrix[k][j] == 0) && (k < first->iRows); k++){}
			if(((first->matrix[k][i] == 0) && (first->matrix[k][j] != 0)) || ((first->matrix[k][i] != 0) && (first->matrix[k][j] == 0))) {mult = 0;}
			
			for(index = k; (k < (first->iRows-1)) && mult; k++){
				multiple1 = (long long int) first->matrix[index][i] * first->matrix[k+1][j];
				multiple2 = (long long int) first->matrix[index][j] * first->matrix[k+1][i];
				mult = (multiple1 == multiple2) ? 1 : 0;
			}
			if(mult == 0) {first->original_c[j] = 0;}
			else if(mult && (((long long int) first->matrix[index][i] * first->matrix[index][j]) < 0)) {first->original_c[j] = -(i+1); counter_c++; copy_mtx_col(first, &i, &j); printf("col %d was subtracted from %d\n", j+1, i+1);}
			else {first->original_c[j] = (i+1); counter_c++; copy_mtx_col(first, &i, &j); printf("col %d was added to %d\n", j+1, i+1);}
		}
	}
	second->iCols_reduced = second->iCols_reduced-counter_c;
}

void convert_mtx_to_vec(item* first, item_calc* second){
	int i, j, iShorter; //iShorter is the number of rows or columns, whichever is less
	second->mtx_as_vec = (int_type*) calloc(second->iRows_reduced * second->iCols_reduced, sizeof(int_type));
	if(second->mtx_as_vec == NULL) {perror("Memory allocation failed for second->mtx_as_vec in function convert_mtx_to_vec!\n");}
	if( first->n_original == 1 && second->iRows_reduced > second->iCols_reduced ){ //In this if else sequence the code transposes the matrix if necessary and transforms the matrix into a vector
		for(j = 0; j < second->iCols_reduced; j++){
			for(i = 0; i < second->iRows_reduced; i++){
				second->mtx_as_vec[j * second->iRows_reduced + i] = first->matrix[i][j];
			}
		}
		iShorter = second->iCols_reduced;
		second->iCols_reduced = second->iRows_reduced;
		second->iRows_reduced = iShorter;
		first->original = first->original_c;
		first->original_length = first->iCols;
		printf("Matrix was transposed.\n");
	}
	else{
		for(i = 0; i < second->iRows_reduced; i++){
			for(j = 0; j < second->iCols_reduced; j++){
				second->mtx_as_vec[i * second->iCols_reduced + j] = first->matrix[i][j];
			}
		}
		first->original = first->original_r;
		first->original_length = first->iRows;
	}
	
	if(second->marginal == 1) {
		int_type *row;
	        row = (int_type*) calloc(second->iCols_reduced, sizeof(int_type));
	        if(row == NULL) {perror("Memory allocation failed for row in function convert_mtx_to_vec!\n");}
	        for(i = 0; i < second->iCols_reduced; i++) {
	        	row[i] = second->mtx_as_vec[(second->iRows_reduced - 1) * second->iCols_reduced + i];
	        	second->mtx_as_vec[(second->iRows_reduced - 1) * second->iCols_reduced + i] = second->mtx_as_vec[i];
	        	second->mtx_as_vec[i] = row[i];
	        }
	        free(row);
	}
}

void mtx_free(item* first){
	int i;
	for(i = 0; i < first->iRows; i++){
		free(first->matrix[i]);
	}
	free(first->matrix);
}

void free_first(item* first){
	free(first->original_r);
	free(first->original_c);
	mtx_free(first);
}

void free_second(item_calc* second){
	free(second->strategy);
	free(second->mtx_as_vec);
}

int* calc_Pattern(int* n){
	int i, *iPattern;
	iPattern = (int*) calloc(2 * *n, sizeof(int));
	if(iPattern == NULL) {perror("Memory allocation failed for iPattern in function calc_Pattern!\n");}
	for(i = 0; i < *n; i++){
		iPattern[i] = i;
		iPattern[2 * *n-i-1]=i;
	}
	return iPattern;
}

unsigned long long int* calc_iNumPower(item_calc* second){
	int i;
	unsigned long long int* iNumPower;
	iNumPower = (unsigned long long int*) malloc(second->maxRows * sizeof(unsigned long long int));
	if(iNumPower == NULL) {perror("Memory allocation failed for iNumPower in function calc_iNumPower!\n");}
	iNumPower[0] = 1;
	for(i = 1; i < (second->maxRows-1); i++){iNumPower[i] = second->n * iNumPower[i-1]; }
	return iNumPower;
}

void print_results(item* first, item_calc* second){
	int i, j, val, *strategy;
	char fileOutput[1024]; // The variable 'fileOutput' is the name of the file to which the strategy vector found to be optimal is written
	FILE *fp;
	strategy = (int*) calloc(first->original_length, sizeof(int));
	if(strategy == NULL) {perror("Memory allocation failed for strategy in function print_results!\n");}
	
	if(first->n_original == 1 && first->marginal == 1){
		printf("LM is: %d\n", second->Lnorm); // Write out the value of the L norm to the screen.
		sprintf(fileOutput,"strategy_LM.txt");
	}
	else {
		printf("L%d is: %d\n", first->n_original, second->Lnorm); // Write out the value of the L norm to the screen.
		sprintf(fileOutput,"strategy_L%d.txt", first->n_original);
	}
	fp = fopen(fileOutput, "w");
	j = 0;
	for(i=0; i<first->original_length; i++) {
		if(first->original[i] == 0) {
			if(second->marginal == 1 && second->n == 1) { 
				if(j == 0) {val = 1;}
				else if(j == (second->iRows_reduced - 1)) {val = second->strategy[0];}
				else {val = second->strategy[j];}
			}
			else{
				if(j == (second->iRows_reduced - 1)) {val = second->n == 1 ? 1 : 0; }
				else {val = second->strategy[j]; }
			}
			j++;
		}
		else {
			if((first->original[i] == first->original_length)) {val = second->n == 1 ? 1 : 0; }//When the entry of the first->original[i] vector is equal to first->original_length, it means that the i-th row consists exclusively zeros.
			else if(first->original[i] < 0) {val = (strategy[-first->original[i] - 1] + 1) & 1; } //When the entry of the first->original[i] vector is negative, it means that the strategy is opposite to the (-first->original[i]-1)-th entry of the strategy vector. As the entries of the strategy vectors are zeros and ones, the current entry of the strategy vector should be negated. 
			else {val = strategy[first->original[i] - 1]; } //When the entry of the first->original[i] vector is positive, it means that the strategy is the same as (first->original[i]-1)-th entry of the strategy vector.
		}
		strategy[i] = val;
		if(second->n == 1) {fprintf(fp, "%d\n", 2 * val - 1);}
		else {fprintf(fp, "%d\n", val);}
	}
	free(strategy);
	fclose(fp);
}
