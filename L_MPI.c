/*****************************

WRITTEN BY ISTVÁN MÁRTON

*****************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <mpi.h>
#include "functions.h"

void calc_jmin_jmax(int* index, unsigned long long int* jMin, unsigned long long int* jMax, unsigned long long int* steps, unsigned long long int* steps_remainder){
	*jMax = (*index + 1) * *steps - 1;
	*jMin = *index * *steps;
	if(*index < *steps_remainder) *jMax += *index + 1;
	else *jMax += *steps_remainder;
	if(*index <= *steps_remainder) *jMin += *index;
	else *jMin += *steps_remainder;
}

void arguments_MPI(item* first, int* argc, char** argv){
	FILE *fp;
	int sd, t;
	char msg[] = "Use the following command: mpirun number_of_threads ./L_MPI filename_of_matrix order_of_the_L_norm";
	if(*argc < 3){
		printf("Incorrect number of input arguments. %s\n", msg);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
		
	sprintf(first->fileName,"%s", argv[1]);
	fp = fopen(first->fileName, "r");
	if(fp == NULL){
		printf("Please make sure that the file containig the matrix exists within this directory. %s\n", msg);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	fclose(fp);
	
	sd = sscanf(argv[2], "%d", &t);
	if((sd == 0) || (t < 1)){
		printf("Please make sure that the order of the L norm is a positive integer. %s\n", msg);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	first->n_original = t;
	
	if(*argc < 4){first->stat = 'n';}
	else {first->stat = argv[3][0];}
	
	if(t > RANK_OF_NORM) {printf("The order of the L norm is too large. Please increase the RANK_OF_NORM variable in the code to %d and compile and run it again.\n", t); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
}

void calc_Parameters(item* first, item_calc* second, int* num_of_threads){

	if(first->n_original < second->iRows_reduced || first->n_original == 1) {second->n = first->n_original;}
	else{
		second->n = second->iRows_reduced > 1 ? second->iRows_reduced : 2;
		printf("The preprocessed matrix has a number of rows (%d) less than or equal to the order of the L norm (%d).\n",second->iRows_reduced ,first->n_original);
	}

	if(second->n == 1){
		first->total_num_to_calc = second->iRows_reduced == 0 ? 0 : (unsigned long long int) 1 << (second->iRows_reduced - 1); //The total number of L sums to be calculated is 2^(iRows_reduced-1)
		second->maxRows = NUM_OF_BITS - 1;
	}
	else if(second->n == 2){
		first->total_num_to_calc = second->iRows_reduced == 0 ? 0 : (unsigned long long int) 1 << (second->iRows_reduced - 1);
		second->maxRows = NUM_OF_BITS - 1;
	}
	else if(second->n == 3){
		first->total_num_to_calc = second->iRows_reduced == 0 ? 0 : pow(3, second->iRows_reduced - 1);
		second->maxRows = (int) (floor (NUM_OF_BITS / log2(second->n)) + 1);
	}
	else{
		first->total_num_to_calc = second->iRows_reduced == 0 ? 0 : pow(second->n, second->iRows_reduced - 1);
		second->maxRows = (int) (floor (NUM_OF_BITS / log2(second->n)) + 1);
	}
	second->copyNum = *num_of_threads > first->total_num_to_calc ? first->total_num_to_calc : *num_of_threads; // copyNum is the actual number of threads. It cannot be more than the number of possible L values.
	second->steps = second->copyNum == 0 ? 0 : first->total_num_to_calc/ second->copyNum;
	second->steps_remainder = second->copyNum == 0 ? 0 : first->total_num_to_calc % second->copyNum;

	printf("number of processes: %llu\n", second->copyNum);
	printf("maximum length of strategy vector: %d\n", second->maxRows);
	if(second->maxRows < second->iRows_reduced){printf("Matrix is too large. The matrix has %d rows after reduction.\n", second->iRows_reduced); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
	if(length < second->iCols_reduced) {printf("Matrix is too large. The length variable %d should be bigger or equal than %d.\n", length, second->iCols_reduced); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
}

void L1(int_type* mtx_as_vec, int* index, unsigned long long int steps, unsigned long long int steps_remainder, int_type* L1_vector, int *L1_strategy, int iShorter, int iLonger){ // This function calculates the L1 norm.
	int i, l, vect[NUM_OF_BITS - 1];
	int_type temp[length], product, L1;
	unsigned long long int number, jMax, jMin, iNumofZeros, aux;

	calc_jmin_jmax(index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.

	number = jMin;
	for(l=0; l < iLonger; l++) {temp[l] = mtx_as_vec[(iShorter - 1) * iLonger + l];} // As the code can consider a row of the matrix with a fixed sign, it considers the last row of the matrix with +1.
	product = 0;
	for(i = 0 ; (iShorter - 1) > i; i++){
		iNumofZeros=(unsigned long long int) 1 << i;
		vect[i] = ((number+ iNumofZeros) >> (i+1)) & 1; // floor((j + 2^i)/2^(i+1)) Logical can be 0 and 1. logical is the number-th word and i-th digit of the BRGC.
		if(vect[i] == 1){for(l=0; l < iLonger; l++){temp[l] += mtx_as_vec[i * iLonger + l]; }} // The code determines the vector-matrix multiplication belonging to the number-th word of the BRGC.
		else {for(l=0; l < iLonger; l++){temp[l] -= mtx_as_vec[i * iLonger + l]; }}				
	}
	for(l= 0; l < iLonger; l++) {product += abs(temp[l]);} // The code calculates the L1 Bell value belonging to the number-th word of the BRGC.
	L1 = product; 
	for(l=0; l< (iShorter - 1); l++){L1_strategy[l] = vect[l];} // The program stores the strategy vector belonging to the number-th BRGC word in the L1_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){ //The code determines the BRGC words until number variable reaches jMax.
		product = 0;
		aux = number;
		for(i = 0; (aux & 1) == 0; i++){
			aux = aux >> 1;
		}
		if(vect[i] == 0){vect[i] = 1; for(l=0; l < iLonger; l++){temp[l] += 2 * mtx_as_vec[i * iLonger + l]; product += abs(temp[l]);}} // When the i-th digit is changed, the code changes the result of the vector-matrix multiplication. It only needs to deal with the i-th row of the matrix.
		else {vect[i] = 0; for(l=0; l < iLonger; l++){temp[l] -= 2 * mtx_as_vec[i * iLonger + l]; product += abs(temp[l]);}}
		if(product > L1) {
			L1 = product; // If the current L1 sum, stored in product, is greater than the previous one, it modifies both the value 
			for(l=0; l<(iShorter - 1); l++){L1_strategy[l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	*L1_vector = L1; // Every thread writes the biggest found L1 sum to L1_vector.
}

void L2(int_type* mtx_as_vec, int* index, unsigned long long int steps, unsigned long long int steps_remainder, int_type* L2_vector, int *L2_strategy, int iRows, int iCols){ // This function calculates the L2 norm.
	int i, l, vect[NUM_OF_BITS - 1];
	int_type temp_0[length], temp_1[length], product, L2;
	unsigned long long int number, jMax, jMin, iNumofZeros, aux;

	calc_jmin_jmax(index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) {temp_0[l] = mtx_as_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; } // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
	product = 0;
	for(i = 0 ; (iRows-1) > i; i++){
		iNumofZeros=(unsigned long long int) 1 << i;
		vect[i] = ((number+ iNumofZeros) >> (i+1)) & 1; // floor((j + 2^i)/2^(i+1)) Logical can be 0 and 1. logical is the number-th word and i-th digit of the BRGC.
			if(vect[i] == 1){for(l=0; l < iCols; l++){temp_1[l] += mtx_as_vec[i * iCols + l]; }}
			else {for(l=0; l < iCols; l++){temp_0[l] += mtx_as_vec[i * iCols + l]; }}				
	}

	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]);}
	L2 = product; // The code calculates the first L2 sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){L2_strategy[l] = vect[l];} // The code writes the first possible strategy vector into the L2_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
		product = 0;
		aux = number;
		for(i = 0; (aux & 1) == 0; i++){
			aux = aux >> 1;
		}
		if(vect[i] == 0){vect[i]=1; for(l=0; l < iCols; l++){temp_1[l] += mtx_as_vec[i * iCols + l]; temp_0[l] -= mtx_as_vec[i * iCols + l]; product += abs(temp_0[l]) + abs(temp_1[l]);}} // When the i-th digit is changed, the code changes the result defined by the definition of the L2 norm. It only needs to deal with the i-th row of the matrix.
		else {vect[i]=0; for(l=0; l < iCols; l++){temp_1[l] -= mtx_as_vec[i * iCols + l]; temp_0[l] += mtx_as_vec[i * iCols + l]; product += abs(temp_0[l]) + abs(temp_1[l]);}}	
		if(product > L2) {
			L2 = product; // If the current L2 sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){L2_strategy[l] = vect[l];} // and the corresponding strategy vector as well.
		}
    }
	*L2_vector = L2; // Every thread writes the biggest found L2 sum to the L2_vector.
}

void L3(int_type* mtx_as_vec, int* index, unsigned long long int steps, unsigned long long int steps_remainder, int_type* L3_vector, int *L3_strategy, int iRows, int iCols, unsigned long long int *iNumPower){ // This function calculates the L3 norm.
	int i, l, iPattern[6] = {0, 1, 2, 2, 1, 0}, vect[NUM_OF_BITS - 1]; // iPattern describes the pattern of the ternary reflected Gray code (TRGC).

	int_type temp_0[length], temp_1[length], temp_2[length], product, L3, temporary;
	unsigned long long int number, jMax, jMin, divide, logical;

	calc_jmin_jmax(index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) {temp_0[l] = mtx_as_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; temp_2[l] = 0;} // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/iNumPower[i]) % 6; // Determines the ternary reflected Gray code (TRGC). iNumPower is a vector consisting of the powers of 3.
		vect[i] = iPattern[logical]; // vect is the possible strategy vector. Its elements consist of 0, +1 or +2.
		switch(vect[i]) {
			case 0:
				for(l=0; l < iCols; l++){temp_0[l] += mtx_as_vec[i * iCols + l]; }
				break;
			case 1:
				for(l=0; l < iCols; l++){temp_1[l] += mtx_as_vec[i * iCols + l]; }
				break;
			case 2:
				for(l=0; l < iCols; l++){temp_2[l] += mtx_as_vec[i * iCols + l]; }
				break; // Every time the code finds the change in the TRGC, it stops searching for further changes.
		}				
	}
	
	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);}
	L3 = product; // The code calculates the first L3 sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){L3_strategy[l] = vect[l];} // The code writes the first possible strategy vector into the L3_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
      	product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/iNumPower[i];
			logical = divide % 3; // The code determines if there is a change in the i-th digit in the TRGC.
			if(logical) {
				logical = divide % 6; // If there is a change, the code calculates the value of the TRGC at that position.
				temporary = iPattern[logical];
				if( (vect[i] == 0)  && (temporary == 1) ) {for(l=0; l < iCols; l++){temp_0[l] -= mtx_as_vec[i * iCols + l]; temp_1[l] += mtx_as_vec[i * iCols + l]; }}// When the i-th digit is changed, the code changes the result defined by the definition of the L3 norm. It only needs to deal with the i-th row of the matrix.
				else if((vect[i] == 1)  && (temporary == 2)) {for(l=0; l < iCols; l++){temp_1[l] -= mtx_as_vec[i * iCols + l]; temp_2[l] += mtx_as_vec[i * iCols + l]; }}
				else if((vect[i] == 2)  && (temporary == 1)) {for(l=0; l < iCols; l++){temp_1[l] += mtx_as_vec[i * iCols + l]; temp_2[l] -= mtx_as_vec[i * iCols + l]; }}
				else {for(l=0; l < iCols; l++){temp_0[l] += mtx_as_vec[i * iCols + l]; temp_1[l] -= mtx_as_vec[i * iCols + l]; }}
				vect[i] = temporary;
				break; // Every time the code finds the change in the d-ary Gray code, it stops searching for further changes.
			}				
		}
		
		for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);} // The code calculates the next L3 sum.
		if(product > L3) {
			L3 = product; // If the current L3 sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){L3_strategy[l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	*L3_vector = L3; // Every thread writes the biggest found L3 sum to the L3_vector.
}

void Ln(int_type* mtx_as_vec, int* index, int* iPattern, unsigned long long int steps, unsigned long long int steps_remainder, int_type* Ln_vector, int *Ln_strategy, int iRows, int iCols, int n, unsigned long long int *iNumPower){ // This function calculates the Ld norm.
	int i, l, vect[NUM_OF_BITS - 1], temporary;
	int_type temp[RANK_OF_NORM][length], product, Ln;
	unsigned long long int number, jMax, jMin, divide, logical;

	calc_jmin_jmax(index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) { // Initializes the temp variable.
		temp[0][l] = mtx_as_vec[(iRows-1) * iCols + l]; // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
		for(i=1; i<n; i++){
			temp[i][l] = 0;
		}
	}
	
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/iNumPower[i]) % (2*n); 
		vect[i] = iPattern[logical]; //iPattern is a vector consisting of the powers of n. It helps determine the words of the n-ary Gray code.
		for(l=0; l < iCols; l++){temp[vect[i]][l] += mtx_as_vec[i * iCols + l]; }				
	}

	for(l= 0; l < iCols; l++) {
		for(i=0; i < n; i++){
			product += abs(temp[i][l]);
		}
	}
	Ln = product;  // The code calculates the first Ln sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){Ln_strategy[l] = vect[l];} // The code writes the first possible strategy vector into the Ln_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
		product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/iNumPower[i];
			logical = divide % n; // The code calculates if there is a change in the i-th digit in the n-ary Gray code.
			if(logical) {
				logical = divide % (2*n); // If there is a change, the code determines the value of the n-ary Gray code at that position.
				temporary = iPattern[logical];
				for(l=0; l < iCols; l++) {temp[vect[i]][l] -= mtx_as_vec[i * iCols + l]; temp[temporary][l] += mtx_as_vec[i * iCols + l]; } // When the i-th digit is changed, the code changes the result defined by the definition of the Ln norm. It only needs to deal with the i-th row of the matrix.
				vect[i] = temporary;
				break;
			}			
		}
		for(l= 0; l < iCols; l++) { // The code calculates the next Ln sum.
			for(i=0; i < n; i++){product += abs(temp[i][l]);}
		}
		if(product > Ln) {
			Ln = product; // If the current Ln sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){Ln_strategy[l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	*Ln_vector = Ln; // Every thread writes the biggest found Ln sum to the Ln_vector.
}

void calc_Lnorm(item_calc* second, int* index) {
	int i, iMax, Ln_vector;

	if(second->n == 1) { // If the order of the L norm is 1 then this part of the code will be executed.
		L1(second->mtx_as_vec, index, second->steps, second->steps_remainder, &Ln_vector, second->strategy, second->iRows_reduced, second->iCols_reduced);
	}
	else if(second->n == 2){ // if the order of the L norm is 2 then this part of the code will be executed.
		L2(second->mtx_as_vec, index, second->steps, second->steps_remainder, &Ln_vector, second->strategy, second->iRows_reduced, second->iCols_reduced);
	}
	else if(second->n == 3){ // if the order of the L norm is 3 then this part of the code will be executed.
		unsigned long long int *iNumPower;
		iNumPower = calc_iNumPower(second);
		L3(second->mtx_as_vec, index, second->steps, second->steps_remainder, &Ln_vector, second->strategy, second->iRows_reduced, second->iCols_reduced, iNumPower);
		free(iNumPower);
	}
	else{ // if the order of the L norm is greater than 3, then this part of the code will be executed.
		unsigned long long int *iNumPower;
		int *iPattern; // iPattern describes the d-ary Gray code.
		iPattern = calc_Pattern(&(second->n));
		iNumPower = calc_iNumPower(second);
		Ln(second->mtx_as_vec, index, iPattern, second->steps, second->steps_remainder, &Ln_vector, second->strategy, second->iRows_reduced, second->iCols_reduced, second->n, iNumPower);
		free(iNumPower);
		free(iPattern);
	}

	if(*index != 0) {
		MPI_Send(&Ln_vector, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	} else {
		second->Lnorm = Ln_vector;
		iMax = 0;
		for(i = 1; i < second->copyNum; i++) {
			MPI_Recv(&Ln_vector, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(second->Lnorm < Ln_vector){
				second->Lnorm = Ln_vector;
				iMax = i;
			}
		}
	}
	MPI_Bcast(&iMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(*index == iMax) {MPI_Send(second->strategy, second->iRows_reduced - 1, MPI_INT, 0, 0, MPI_COMM_WORLD);}
	if(*index == 0) {
		MPI_Recv(second->strategy, second->iRows_reduced - 1, MPI_INT, iMax, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} //Write out the strategy belonging to the L norm to file.
}

void load_parameters(item* first, item_calc* second, int* argc, char** argv, int* num_of_threads){
	arguments_MPI(first, argc, argv);
	matrix_read(first);
	eliminate_zero_rows_cols(first, second);
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R') {calc_reduce_matrix_rows(first, second); }
	if(first->stat == 'y' || first->stat == 'c' || first->stat == 'Y' || first->stat == 'C') {calc_reduce_matrix_cols(first, second); if(first->n_original > 1) {calc_reduce_matrix_cols_sign(first, second);} }
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R' || first->stat == 'c' || first->stat == 'C') { delete_rows(first, second); delete_cols(first, second);}
	convert_mtx_to_vec(first, second);
	calc_Parameters(first, second, num_of_threads); // Calculates the necessary parameters for the calculation.
}

void create_mpi_datatype(item_calc* second, MPI_Datatype *mpi_item_calc){
	int lengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Aint displacements[8];
	MPI_Datatype types[8] = { MPI_INT, MPI_INT, MPI_INT,  MPI_INT, MPI_INT, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG };
	MPI_Aint base_address;
	MPI_Get_address(second, &base_address);
	MPI_Get_address(&second->n, &displacements[0]);
	MPI_Get_address(&second->iRows_reduced, &displacements[1]);
	MPI_Get_address(&second->iCols_reduced, &displacements[2]);
	MPI_Get_address(&second->maxRows, &displacements[3]);
	MPI_Get_address(&second->Lnorm, &displacements[4]);
	MPI_Get_address(&second->steps, &displacements[5]);
	MPI_Get_address(&second->steps_remainder, &displacements[6]);
	MPI_Get_address(&second->copyNum, &displacements[7]);
	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);
	displacements[4] = MPI_Aint_diff(displacements[4], base_address);
	displacements[5] = MPI_Aint_diff(displacements[5], base_address);
	displacements[6] = MPI_Aint_diff(displacements[6], base_address);
	displacements[7] = MPI_Aint_diff(displacements[7], base_address);
	MPI_Type_create_struct(8, lengths, displacements, types, mpi_item_calc);
	MPI_Type_commit(mpi_item_calc);
}

int main(int argc, char *argv[]) {
	item first;
	item_calc second;
	int index, num_of_threads;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &index);
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_threads);
	if (index == 0) {
		load_parameters(&first, &second, &argc, argv, &num_of_threads);
	}

	MPI_Datatype mpi_item_calc;
	create_mpi_datatype(&second, &mpi_item_calc);
	MPI_Bcast(&second, 1, mpi_item_calc, 0, MPI_COMM_WORLD);
	
	if(index != 0){
		second.mtx_as_vec = (int*) malloc(second.iRows_reduced * second.iCols_reduced * sizeof(int));
	}
	second.strategy = (int*)calloc(second.iRows_reduced - 1, sizeof(int));
	MPI_Bcast(second.mtx_as_vec, second.iRows_reduced * second.iCols_reduced, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(second.iRows_reduced == 0 ) {
		if(index == 0) {printf("This is a zero matrix.\n"); second.Lnorm = 0; }
	}
	else{
		if(index < second.copyNum) {calc_Lnorm(&second, &index);} // The function 'calc_Lnorm' calculates the L norm of order n of the input matrix.
	}
	if(index == 0){
		print_results(&first, &second);
		free_first(&first);
	}
	if(second.iRows_reduced > 0) free_second(&second);
	MPI_Type_free(&mpi_item_calc);
	MPI_Finalize();
	return 0;
}
