/*****************************

WRITTEN BY ISTVÁN MÁRTON

*****************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "functions.h"

__device__ void calc_jmin_jmax(int* index, unsigned long long int* jMin, unsigned long long int* jMax, unsigned long long int* steps, unsigned long long int* steps_remainder){
	*jMax = (*index + 1) * *steps - 1;
	*jMin = *index * *steps;
	if(*index < *steps_remainder) *jMax += *index + 1;
	else *jMax += *steps_remainder;
	if(*index <= *steps_remainder) *jMin += *index;
	else *jMin += *steps_remainder;
}

__global__ void L1(int_type* d_mtx_as_vec, unsigned long long int steps, unsigned long long int steps_remainder, int_type *d_L1_vector, int *d_L1_strategy, int iShorter, int iLonger){ // This function calculates the L1 norm.
	int i, l, index, vect[NUM_OF_BITS - 1];
	int_type temp[length], product, L1;
	unsigned long long int number, jMax, jMin, iNumofZeros, aux;

	index = blockIdx.x * blockDim.x + threadIdx.x; // Index of threads.
	calc_jmin_jmax(&index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.

	number = jMin;
	for(l=0; l < iLonger; l++) {temp[l] = d_mtx_as_vec[(iShorter - 1) * iLonger + l];} // As the code can consider a row of the matrix with a fixed sign, it considers the last row of the matrix with +1.
	product = 0;
	for(i = 0 ; (iShorter - 1) > i; i++){
		iNumofZeros=(unsigned long long int) 1 << i;
		vect[i] = ((number+ iNumofZeros) >> (i+1)) & 1; // floor((j + 2^i)/2^(i+1)) Logical can be 0 and 1. logical is the number-th word and i-th digit of the BRGC.
		if(vect[i] == 1){for(l=0; l < iLonger; l++){temp[l] += d_mtx_as_vec[i * iLonger + l]; }} // The code determines the vector-matrix multiplication belonging to the number-th word of the BRGC.
		else {for(l=0; l < iLonger; l++){temp[l] -= d_mtx_as_vec[i * iLonger + l]; }}				
	}
	for(l= 0; l < iLonger; l++) {product += abs(temp[l]);} // The code calculates the L1 Bell value belonging to the number-th word of the BRGC.
	L1 = product; 
	for(l=0; l< (iShorter - 1); l++){d_L1_strategy[index * (iShorter - 1) + l] = vect[l];} // The program stores the strategy vector belonging to the number-th BRGC word in the d_L1_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){ //The code determines the BRGC words until number variable reaches jMax.
		product = 0;
		aux = number;
		for(i = 0; (aux & 1) == 0; i++){
			aux = aux >> 1;
		}
		if(vect[i] == 0){vect[i] = 1; for(l=0; l < iLonger; l++){temp[l] += 2 * d_mtx_as_vec[i * iLonger + l]; product += abs(temp[l]);}} // When the i-th digit is changed, the code changes the result of the vector-matrix multiplication. It only needs to deal with the i-th row of the matrix.
		else {vect[i] = 0; for(l=0; l < iLonger; l++){temp[l] -= 2 * d_mtx_as_vec[i * iLonger + l]; product += abs(temp[l]);}}
		if(product > L1) {
			L1 = product; // If the current L1 sum, stored in product, is greater than the previous one, it modifies both the value 
			for(l=0; l<(iShorter - 1); l++){d_L1_strategy[index * (iShorter - 1) + l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	d_L1_vector[index] = L1; // Every thread writes the biggest found L1 sum to d_L1_vector.
}

__global__ void L2(int_type* d_mtx_as_vec, unsigned long long int steps, unsigned long long int steps_remainder, int_type *d_L2_vector, int *d_L2_strategy, int iRows, int iCols){ // This function calculates the L2 norm.
	int i, l, index, vect[NUM_OF_BITS - 1];
	int_type temp_0[length], temp_1[length], product, L2;
	unsigned long long int number, jMax, jMin, iNumofZeros, aux;

	index = blockIdx.x * blockDim.x + threadIdx.x; // Index of threads.

	calc_jmin_jmax(&index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) {temp_0[l] = d_mtx_as_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; } // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
	product = 0;
	for(i = 0 ; (iRows-1) > i; i++){
		iNumofZeros=(unsigned long long int) 1 << i;
		vect[i] = ((number+ iNumofZeros) >> (i+1)) & 1; // floor((j + 2^i)/2^(i+1)) Logical can be 0 and 1. logical is the number-th word and i-th digit of the BRGC.
			if(vect[i] == 1){for(l=0; l < iCols; l++){temp_1[l] += d_mtx_as_vec[i * iCols + l]; }}
			else {for(l=0; l < iCols; l++){temp_0[l] += d_mtx_as_vec[i * iCols + l]; }}				
	}

	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]);}
	L2 = product; // The code calculates the first L2 sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){d_L2_strategy[index * (iRows - 1) + l] = vect[l];} // The code writes the first possible strategy vector into the d_L2_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
		product = 0;
		aux = number;
		for(i = 0; (aux & 1) == 0; i++){
			aux = aux >> 1;
		}
		if(vect[i] == 0){vect[i]=1; for(l=0; l < iCols; l++){temp_1[l] += d_mtx_as_vec[i * iCols + l]; temp_0[l] -= d_mtx_as_vec[i * iCols + l]; product += abs(temp_0[l]) + abs(temp_1[l]);}} // When the i-th digit is changed, the code changes the result defined by the definition of the L2 norm. It only needs to deal with the i-th row of the matrix.
		else {vect[i]=0; for(l=0; l < iCols; l++){temp_1[l] -= d_mtx_as_vec[i * iCols + l]; temp_0[l] += d_mtx_as_vec[i * iCols + l]; product += abs(temp_0[l]) + abs(temp_1[l]);}}	
		if(product > L2) {
			L2 = product; // If the current L2 sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){d_L2_strategy[index * (iRows - 1) + l] = vect[l];} // and the corresponding strategy vector as well.
		}
    }
	d_L2_vector[index] = L2; // Every thread writes the biggest found L2 sum to the d_L2_vector.
}

__global__ void L3(int_type* d_mtx_as_vec, unsigned long long int steps, unsigned long long int steps_remainder, int_type *d_L3_vector, int *d_L3_strategy, int iRows, int iCols, unsigned long long int *d_iNumPower){ // This function calculates the L3 norm.
	int i, l, index, vect[NUM_OF_BITS - 1], temporary, d_iPattern[6] = {0, 1, 2, 2, 1, 0}; // d_iPattern describes the pattern of the ternary reflected Gray code (TRGC).

	int_type temp_0[length], temp_1[length], temp_2[length], product, L3;
	unsigned long long int number, jMax, jMin, divide, logical;

	index = blockIdx.x * blockDim.x + threadIdx.x; // Index of threads.

	calc_jmin_jmax(&index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) {temp_0[l] = d_mtx_as_vec[(iRows-1) * iCols + l]; temp_1[l] = 0; temp_2[l] = 0;} // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/d_iNumPower[i]) % 6; // Determines the ternary reflected Gray code (TRGC). d_iNumPower is a vector consisting of the powers of 3.
		vect[i] = d_iPattern[logical]; // vect is the possible strategy vector. Its elements consist of 0, +1 or +2.
		switch(vect[i]) {
			case 0:
				for(l=0; l < iCols; l++){temp_0[l] += d_mtx_as_vec[i * iCols + l]; }
				break;
			case 1:
				for(l=0; l < iCols; l++){temp_1[l] += d_mtx_as_vec[i * iCols + l]; }
				break;
			case 2:
				for(l=0; l < iCols; l++){temp_2[l] += d_mtx_as_vec[i * iCols + l]; }
				break; // Every time the code finds the change in the TRGC, it stops searching for further changes.
		}				
	}
	
	for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);}
	L3 = product; // The code calculates the first L3 sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){d_L3_strategy[index * (iRows - 1) + l] = vect[l];} // The code writes the first possible strategy vector into the d_L3_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
      	product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/d_iNumPower[i];
			logical = divide % 3; // The code determines if there is a change in the i-th digit in the TRGC.
			if(logical) {
				logical = divide % 6; // If there is a change, the code calculates the value of the TRGC at that position.
				temporary = d_iPattern[logical];
				if( (vect[i] == 0)  && (temporary == 1) ) {for(l=0; l < iCols; l++){temp_0[l] -= d_mtx_as_vec[i * iCols + l]; temp_1[l] += d_mtx_as_vec[i * iCols + l]; }}// When the i-th digit is changed, the code changes the result defined by the definition of the L3 norm. It only needs to deal with the i-th row of the matrix.
				else if((vect[i] == 1)  && (temporary == 2)) {for(l=0; l < iCols; l++){temp_1[l] -= d_mtx_as_vec[i * iCols + l]; temp_2[l] += d_mtx_as_vec[i * iCols + l]; }}
				else if((vect[i] == 2)  && (temporary == 1)) {for(l=0; l < iCols; l++){temp_1[l] += d_mtx_as_vec[i * iCols + l]; temp_2[l] -= d_mtx_as_vec[i * iCols + l]; }}
				else {for(l=0; l < iCols; l++){temp_0[l] += d_mtx_as_vec[i * iCols + l]; temp_1[l] -= d_mtx_as_vec[i * iCols + l]; }}
				vect[i] = temporary;
				break; // Every time the code finds the change in the d-ary Gray code, it stops searching for further changes.
			}				
		}
		
		for(l= 0; l < iCols; l++) {product += abs(temp_0[l]) + abs(temp_1[l]) + abs(temp_2[l]);} // The code calculates the next L3 sum.
		if(product > L3) {
			L3 = product; // If the current L3 sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){d_L3_strategy[index * (iRows - 1) + l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	d_L3_vector[index] = L3; // Every thread writes the biggest found L3 sum to the d_L3_vector.
}

__global__ void Ln(int_type* d_mtx_as_vec, int* d_iPattern, unsigned long long int steps, unsigned long long int steps_remainder, int_type *d_Ln_vector, int *d_Ln_strategy, int iRows, int iCols, int n, unsigned long long int *d_iNumPower){ // This function calculates the Ld norm.
	int i, l, index, vect[NUM_OF_BITS - 1], temporary;
	int_type temp[RANK_OF_NORM][length], product, Ln;
	unsigned long long int number, jMax, jMin, divide, logical;

	index = blockIdx.x * blockDim.x + threadIdx.x; // Index of threads.

	calc_jmin_jmax(&index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	
	number = jMin;
	for(l=0; l < iCols; l++) { // Initializes the temp variable.
		temp[0][l] = d_mtx_as_vec[(iRows-1) * iCols + l]; // As the code can consider a row of the matrix with a fixed label, it considers the last row of the matrix with 0.
		for(i=1; i<n; i++){
			temp[i][l] = 0;
		}
	}
	
	product = 0;
	for(i = 0 ; (iRows - 1) > i; i++){
		logical = (number/d_iNumPower[i]) % (2*n); 
		vect[i] = d_iPattern[logical]; //d_iPattern is a vector consisting of the powers of n. It helps determine the words of the n-ary Gray code.
		for(l=0; l < iCols; l++){temp[vect[i]][l] += d_mtx_as_vec[i * iCols + l]; }				
	}

	for(l= 0; l < iCols; l++) {
		for(i=0; i < n; i++){
			product += abs(temp[i][l]);
		}
	}
	Ln = product;  // The code calculates the first Ln sum associated with the thread.
	for(l=0; l<(iRows - 1); l++){d_Ln_strategy[index * (iRows - 1) + l] = vect[l];} // The code writes the first possible strategy vector into the d_Ln_strategy vector.

	for(number=jMin + 1; number <= jMax; number++){
		product = 0;
		for(i = 0 ; (iRows - 1) > i; i++){
			divide = number/d_iNumPower[i];
			logical = divide % n; // The code calculates if there is a change in the i-th digit in the n-ary Gray code.
			if(logical) {
				logical = divide % (2*n); // If there is a change, the code determines the value of the n-ary Gray code at that position.
				temporary = d_iPattern[logical];
				for(l=0; l < iCols; l++) {temp[vect[i]][l] -= d_mtx_as_vec[i * iCols + l]; temp[temporary][l] += d_mtx_as_vec[i * iCols + l]; } // When the i-th digit is changed, the code changes the result defined by the definition of the Ln norm. It only needs to deal with the i-th row of the matrix.
				vect[i] = temporary;
				break;
			}			
		}
		for(l= 0; l < iCols; l++) { // The code calculates the next Ln sum.
			for(i=0; i < n; i++){product += abs(temp[i][l]);}
		}
		if(product > Ln) {
			Ln = product; // If the current Ln sum is greater than the previous one, it modifies both the value
			for(l=0; l<(iRows - 1); l++){d_Ln_strategy[index * (iRows - 1) + l] = vect[l];} // and the corresponding strategy vector as well.
		}
	}
	d_Ln_vector[index] = Ln; // Every thread writes the biggest found Ln sum to the d_Ln_vector.
}

void calc_Parameters(item* first, item_calc* second, int* num_of_blocks, int* num_of_threads_per_block){

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

	*num_of_threads_per_block = *num_of_threads_per_block > first->total_num_to_calc ? first->total_num_to_calc : *num_of_threads_per_block;
	if(*num_of_threads_per_block == 0) {*num_of_blocks = 0;}
	else {*num_of_blocks = first->total_num_to_calc/ *num_of_threads_per_block > *num_of_blocks ? *num_of_blocks : first->total_num_to_calc/ *num_of_threads_per_block;}
	second->copyNum = *num_of_blocks * *num_of_threads_per_block;

	second->steps = second->copyNum == 0 ? 0 : first->total_num_to_calc/ second->copyNum;
	second->steps_remainder = second->copyNum == 0 ? 0 : first->total_num_to_calc % second->copyNum;

	printf("num_ofBlock: %d, num_ofThread in block: %d\n", *num_of_blocks, *num_of_threads_per_block);
	printf("maximum length of strategy vector: %d\n", second->maxRows);
	if(second->maxRows < second->iRows_reduced){printf("Matrix is too large. The matrix has %d rows after reduction. The maximum number of rows cannot be more than %d.\n", second->iRows_reduced, second->maxRows); exit(-1);}
	if(second->iCols_reduced > length) {printf("Matrix is too large. The length variable %d should be greater or equal than %d.\n", length, second->iCols_reduced); exit(-1);}
}

void calc_Lnorm(item_calc* second, int* num_of_blocks, int* num_of_threads_per_block){
	int i, iMax, *Ln_strategy, *d_Ln_strategy; //iMax: the index of the strategy vector found to be optimal; mtx_as_vec: the input matrix is converted to a vector in the host; d_mtx_as_vec: the converted matrix in the device; maxRows: them maximal number of rows (in the case of L1, the maximal number of rows or columns) of the matrix the program can deal with (this number is determined by the order of the L norm that should be calculated); Ln_vector and d_Ln_vector are the two vectors containing the possible L norms belonging to a given thread in the host and device, respectively; Ln_strategy and d_Ln_strategy are the two vector containing all of the possible strategy vectors belonging to different threads; num_ofBlock is the number of blocks the program uses; num_ofThread is the number of threads in a block
	int_type *d_mtx_as_vec, *Ln_vector, *d_Ln_vector;
	Ln_vector = (int_type*)calloc(second->copyNum, sizeof(int_type)); // The code allocates memory in the host for the possible L norms.
	Ln_strategy = (int*) malloc(second->copyNum * (second->iRows_reduced - 1) * sizeof(int)); // The code allocates memory for the strategies belonging to the possible L norms in the host.
	second->strategy = (int*) calloc(second->iRows_reduced - 1, sizeof(int));
	cudaMalloc((void**)&d_mtx_as_vec, second->iRows_reduced * second->iCols_reduced * sizeof(int_type)); // Allocating memory for the matrix in the device.
	cudaMalloc((void**)&d_Ln_vector, second->copyNum * sizeof(int_type)); // The code allocates memory in the device for the possible L norms.
	cudaMalloc((void**)&d_Ln_strategy, second->copyNum * (second->iRows_reduced - 1) * sizeof(int)); // The code allocates memory for the strategies belonging to the possible L norms in the device.
	cudaMemcpy(d_mtx_as_vec, second->mtx_as_vec, second->iRows_reduced * second->iCols_reduced * sizeof(int_type), cudaMemcpyHostToDevice); // The matrix is copied from RAM to GPU memory.
	if(second->n == 1){ // If the order of the L norm is 1 then this part of the code will be executed.
		L1<<<*num_of_blocks, *num_of_threads_per_block>>>(d_mtx_as_vec, second->steps, second->steps_remainder, d_Ln_vector, d_Ln_strategy, second->iRows_reduced, second->iCols_reduced); // The calculation of the L1 norm with GPU.
	}
	else if(second->n == 2){ // If the order of the L norm is 2 then this part of the code will be executed.
		L2<<<*num_of_blocks, *num_of_threads_per_block>>>(d_mtx_as_vec, second->steps, second->steps_remainder, d_Ln_vector, d_Ln_strategy, second->iRows_reduced, second->iCols_reduced);
	}
	else if(second->n == 3){ // If the order of the L norm is 3 then this part of the code will be executed.
		unsigned long long int *iNumPower, *d_iNumPower;// iNumPower is copied to the device memory to speed up the calculation of the d-ary Gray code
		iNumPower = calc_iNumPower(second);
		cudaMalloc((void**)&d_iNumPower, (second->maxRows-1) * sizeof(unsigned long long int));
		cudaMemcpy(d_iNumPower, iNumPower, (second->maxRows-1) * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		L3<<<*num_of_blocks, *num_of_threads_per_block>>>(d_mtx_as_vec, second->steps, second->steps_remainder, d_Ln_vector, d_Ln_strategy, second->iRows_reduced, second->iCols_reduced, d_iNumPower);
		free(iNumPower);
		cudaFree(d_iNumPower);
	}
	else{ // If the order of the L norm is greater than 3, then this part of the code will be executed.
		unsigned long long int *iNumPower, *d_iNumPower;// iNumPower is copied to the device memory to speed up the calculation of the d-ary Gray code
		int *iPattern, *d_iPattern; // iPattern describes the d-ary Gray code.
		iPattern = calc_Pattern(&(second->n));
		iNumPower = calc_iNumPower(second);
		cudaMalloc((void**)&d_iNumPower, (second->maxRows-1) * sizeof(unsigned long long int));
		cudaMalloc((void**)&d_iPattern, 2 * second->n * sizeof(int));
		cudaMemcpy(d_iNumPower, iNumPower, (second->maxRows-1) * sizeof(unsigned long long int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_iPattern, iPattern, 2 * second->n * sizeof(int), cudaMemcpyHostToDevice);
		Ln<<<*num_of_blocks, *num_of_threads_per_block>>>(d_mtx_as_vec, d_iPattern, second->steps, second->steps_remainder, d_Ln_vector, d_Ln_strategy, second->iRows_reduced, second->iCols_reduced, second->n, d_iNumPower);
		free(iNumPower);
		free(iPattern);
		cudaFree(d_iNumPower);
		cudaFree(d_iPattern);
	}
	cudaMemcpy(Ln_vector, d_Ln_vector, second->copyNum * sizeof(int_type), cudaMemcpyDeviceToHost); // Copy the possible L norm values from device to host.
	cudaMemcpy(Ln_strategy, d_Ln_strategy, second->copyNum * (second->iRows_reduced - 1) * sizeof(int), cudaMemcpyDeviceToHost); // Copy the possible strategies belonging to possible L norm values from device to host.
	second->Lnorm = Ln_vector[0];
	iMax = 0; // Determining the maximal element of Ln_vector, which is the L norm, and the index of the corresponding strategy vector as well.
	for(i = 1; i < second->copyNum; i++){ if(second->Lnorm < Ln_vector[i]) {second->Lnorm = Ln_vector[i]; iMax = i;}}
	for(i = 0; i < (second->iRows_reduced - 1); i++){second->strategy[i] = Ln_strategy[iMax * (second->iRows_reduced - 1) + i];}
	
	free(Ln_vector); // Deallocates the vectors in the host memory.
	free(Ln_strategy);
	
	cudaFree(d_Ln_vector); // Deallocates the vectors in the device memory.
	cudaFree(d_Ln_strategy);
	cudaFree(d_mtx_as_vec);
}

void arguments_CUDA(item* first, int* num_of_blocks, int* num_of_threads_per_block, int* argc, char** argv, cudaDeviceProp* devProp){
	FILE *fp;
	int sd, t;
	char msg[] = "Use the following command: ./L_CUDA number_of_blocks number_of_threads_per_block filename_of_matrix order_of_the_L_norm";
	if(*argc < 5){
		printf("Incorrect number of input arguments. %s\n", msg);
		exit(-1);
	}
	
	sd = sscanf(argv[1], "%d", &t);
	if((sd == 0) || (t < 1)){
		printf("Please make sure that the number of blocks is a positive integer. %s\n", msg);
		exit(-1);
	}
	*num_of_blocks = t;
	
	sd = sscanf(argv[2], "%d", &t);
	if((sd == 0) || (t < 1)){
		printf("Please make sure that the number of threads is a positive integer. %s\n", msg);
		exit(-1);
	}
	*num_of_threads_per_block = t;
	
	sprintf(first->fileName,"%s", argv[3]);
	fp = fopen(first->fileName, "r");
	if(fp == NULL){
		printf("Please make sure that the file containig the matrix exists within this directory. %s\n", msg);
	}
	fclose(fp);
	
	sd = sscanf(argv[4], "%d", &t);
	if((sd == 0) || (t < 1)){
		printf("Please make sure that the order of the L norm is a positive integer. %s\n", msg);
		exit(-1);
	}
	first->n_original = t;
	
	if(*argc < 6){first->stat = 'n';}
	else {first->stat = argv[5][0];}
	
	if(t > RANK_OF_NORM) {printf("The order of the L norm is too large. Please increase the RANK_OF_NORM variable in the code to %d and compile and run it again.\n", t); exit(-1);}
	if(*num_of_threads_per_block > devProp->maxThreadsPerBlock) {printf("The maximum number of threads per block cannot be more than %d.\n", devProp->maxThreadsPerBlock); exit(-1);}
}

void load_parameters(item* first, item_calc* second, int* argc, char** argv, cudaDeviceProp* devProp, int* num_of_blocks, int* num_of_threads_per_block){
	arguments_CUDA(first, num_of_blocks, num_of_threads_per_block, argc, argv, devProp);
	matrix_read(first);
	eliminate_zero_rows_cols(first, second);
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R') {calc_reduce_matrix_rows(first, second); }
	if(first->stat == 'y' || first->stat == 'c' || first->stat == 'Y' || first->stat == 'C') {calc_reduce_matrix_cols(first, second); if(first->n_original > 1) {calc_reduce_matrix_cols_sign(first, second);} }
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R' || first->stat == 'c' || first->stat == 'C') { delete_rows(first, second); delete_cols(first, second);}
	if(first->n_original == 1) {convert_mtx_to_vec_Transpose(first, second); }
	else {convert_mtx_to_vec_noTranspose(first, second);}
	calc_Parameters(first, second, num_of_blocks, num_of_threads_per_block); // Calculates the necessary parameters for the calculation.
}

int main(int argc, char *argv[]){
	int num_of_threads_per_block, num_of_blocks;
	item first;
	item_calc second;
	cudaDeviceProp devProp; // devProp contains the number of cores a warp contain.
	cudaGetDeviceProperties(&devProp, 0);
	load_parameters(&first, &second, &argc, argv, &devProp, &num_of_blocks, &num_of_threads_per_block);
	if(second.iRows_reduced == 0) {printf("This is a zero matrix.\n"); second.Lnorm = 0;}
	else {calc_Lnorm(&second, &num_of_blocks, &num_of_threads_per_block);} // The function 'calc_Lnorm' calculates the L norm of order n of the input matrix.
	print_results(&first, &second);
	
	free_first(&first);
	if(second.iRows_reduced > 0) free_second(&second);
	return 0;
}
