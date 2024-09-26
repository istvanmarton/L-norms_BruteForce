#Brute force calculation of the L norms of a matrix

This repository accompanies the article "Classical bounds on correlation-type Bell expressions and linear prepare-and-measure witnesses: efficient computation in parallel environments such as graphics processing units" by I. Márton, E. Bene, G. Drótos. In this article, we present the calculation of the $L_d$ norms, for any d, of an n×m matrix having real entries. The $L_1$ norm is the local bound of the bipartite correlation-type Bell expression represented by the given matrix. The $L_d$ norm for d ≥ 2 is the classical d-dimensional bound of the linear prepare-and-measure witness represented by the given matrix, i.e., the bound of the given witness associated with the communication of d-level classical messages in the prepare-and-measure scenario.

A C implementation utilizing CUDA, OpenMP, and MPI has been made accessible in this repository under the source files L_CUDA.cu, L_OpenMP.c, and L_MPI.c, respectively; for the input, the matrix entries need to be converted to integers by an appropriate rescaling of the matrix (as the $L_d$ norms scale linearly with a multiplicative factor).

##Compilation and execution

The "functions.h" file is needed to compile any of the above codes successfully.

###CUDA

To compile and run the CUDA code, one needs access to a CUDA compatible GPU and to install the CUDA compiler and a GPU driver. The CUDA code can be compiled with the following command:

	nvcc -o L_CUDA L_CUDA.cu

Then the program can be invoked as

	./L_CUDA number_of_blocks number_of_threads_per_block filename_of_matrix d opt_flag_for_preproc

where
 + 'number_of_blocks' is the number of blocks (grid size) the GPU is allowed to utilize,
 + 'number_of_threads_per_block' is the number of threads in one block (the block size),
 + 'filename_of_matrix' is the name of the file containing the input matrix,
 + 'd' is a positive integer indicating the order of the L norm one intends to calculate.

Furthermore, 'opt_flag_for_preproc' stands for an optional character that controls preprocessing of the matrix before the actual computation; please refer to the article for a discussion. When the value of this variable is not specified, the code does not perform any preprocessing. If the value is
 + 'y' or 'Y', the code will search for rows and columns that can be united or have only zero entries,
 + 'c' or 'C', it will search for columns that can be united and rows and columns with only zero entries,
 + 'r' or 'R' it will search for rows that can be united and rows and columns with only zero entries.

If the value of the optional character is different from the previously listed ones, it is equivalent to not specifying any character.

*Example*: The command

	./L_CUDA 512 32 matrix.txt 1

invokes the GPU with 512 blocks and 32 threads in one block. The name of the file containing the matrix is matrix.txt, the order of the L norm is 1, and as the optional character is not specified, the program will search for neither rows nor columns that can be united or have only zero entries. In case one specifies the optional flag as in

	./L_CUDA 512 32 matrix.txt 1 y

the program will preprocess the matrix for both rows or columns to reduce the matrix.

###OpenMP

The OpenMP code can be compiled, e.g., with the following command:

	gcc L_OpenMP.c -fopenmp -o L_OpenMP -lm

Please note that we have tested GCC only. For GCC, the default is no optimization; to specify the optimization level, please use the -O option. After compilation, the program can be invoked as

	./L_OpenMP number_of_threads filename_of_matrix d opt_flag_for_preproc

The argument 'number_of_threads' specifies the number of threads. For the rest of the arguments, please refer to the CUDA version.

###MPI

In order to compile and execute the MPI version, an MPI library should be installed. This code was tested with Open MPI. The MPI version can be compiled, e.g., with the following command:

	mpicc -o L_MPI L_MPI.c -lm

Regarding the optimization level, please refer to the documentation of the compiler to be used. After compilation, the program can be invoked as

	mpirun -n number_of_processes ./L_MPI filename_of_matrix d opt_flag_for_preproc

or

	mpiexec -n number_of_processes ./L_MPI filename_of_matrix d opt_flag_for_preproc

The argument 'number_of_processes' specifies the number of processes. For the rest of the arguments, please refer to the CUDA version.

##Options for tuning in 'functions.h'

&nbsp;
 + The variable 'length' defines the maximal number of columns allowed for the matrix for d > 1. For d = 1, it limits the greater of the number of rows and columns.
 + The maximal order of the L norm, i.e., d, that the code can calculate is defined in the variable 'RANK_OF_NORM'.

Please note: If the product of the variables 'length' and 'RANK_OF_NORM' is too large for a given GPU, the CUDA version may fail silently and return a zero L norm value. For instance, when NVIDIA RTX 6000 Ada was used along with 'length = 4096' and 'RANK_OF_NORM = 20', and the matrix entries were represented as int, the code returned with zero for the $L_4$ norm of the matrix, and the CUDA kernel was not even invoked.

The code is prepared for representing the matrix entries as an arbitrarily chosen integer type. The default is 'int'. Suppose one intends to use 'long long int'. In this case, the following modifications are needed in 'functions.h':
 + line 11: 'typedef int int_type;' -> 'typedef long long int int_type;'
 + line 59 (in function 'matrix_read'): 'sscanf(cNum, "%d", &value);' -> 'sscanf(cNum, "%lld", &value);'
 + line 412 (in function 'print_results'): 'printf("L%d is: %d\n", second->n, second->Lnorm);' -> 'printf("L%d is: %lld\n", second->n, second->Lnorm);'

The abs() functions in *all of the source files* may also need to be modified to llabs(). We need to mention that calculating with 'long long int' can cause significantly longer execution times compared to the case when 'int' is used. Matrix preprocessing can be used when multiplying any two entries of the matrix does not cause overflow, namely the the resulting number can be represented with a long long int.

**Input**

In the file containing the matrix:
 + The separator between entries of the same row can be any character apart from a numeral or any of the following characters: 'e', 'E', '.', '+', '-'.
 + Consecutive rows must be placed in new lines.
 + The entries of the matrix must be integers.

**Output**

Output is shared between the standard output stream (where the actual L norm value is printed among other pieces of information) and a file dedicated to the optimal strategy vector.

*Standard output*

The program prints here
 + the number of rows and columns of the input matrix as 'rows: r, cols: c', where r and c stand for the number of rows and the number of columns of the input matrix, respectively,
 + the number of blocks and the number of threads for the CUDA version, the number of threads for the OpenMP version, and number of processes for the MPI version,
 + the maximal length of the strategy vector the computer can deal with, and
 + the value of the L norm.

Furthermore, in case the program was allowed to preprocess the matrix, the program writes out the performed reduction steps. For instance:
 + 'row 0 entries are zeros.' means the entries of the zeroth row of the matrix are zeros, and that row was deleted.
 + 'col 3 was subtracted from 1' means column 3 was subtracted from column 1, and column 3 was deleted.
 + 'row 5 was added to 1.' means row 5 was added to row 1, and row 5 was deleted.

*Output file*

The program also writes the optimal strategy, defining the L norm of a given order, into a file. The filename has the following format: strategy_Ld.txt, where 'd' means the order of the L norm. For example, the file 'strategy_L2.txt' contains the optimal strategy vector associated with the L_2 norm of the matrix. These files consist of entries ±1 if d = 1. If d ≥ 2, the entries of the strategy vectors in the files are from {0,1,...,d-1}.

##Testing the code for performance

In this section, we describe how to run calculations to figure out the optimal grid and block sizes for the CUDA code. These auxiliary codes are written in bash, Octave, and gnuplot.

*gen_mtx.m*

To test the code for performance, one needs to generate a matrix for which the calculations are performed. The Octave script 'gen_mtx.m' generates a random matrix. It has two important parameters, namely the variables 'n' and 'r'. The script produces a square matrix, the size (number of rows or columns) of which is given by 'n'. The variable 'r' defines the interval, as \[-r,r\], from which the entries of the matrix are taken as random integers. The output of the script will be a file containing a matrix. The name of the matrix describes its size. For instance, if the name of the file is 'mtx17.txt', the size of the matrix is 17×17.

*time_order_1.sh*

The file 'time_order_1.sh' is a template for running the CUDA program and measuring the runtimes. To run this script, 'L_CUDA.cu', 'functions.h', the shell script itself, and the file containing the matrix should be copied to the working directory. The variable 'order' in the shell script stands for the order of the L norm the program will compute. The variable 'iteration' determines the number of times the program will calculate the L norm of the matrix (in order to sample variability in the runtime). This sample code assumes an input matrix file named 'mtx38.txt'; this can be adjusted in line 25. The variables 'threads' and 'blocks' determine the block (number of threads in a block) and grid (number of blocks) sizes the program will use. By default, the 'threads' vector is 1024, 512, etc., and the 'blocks' vector is 8, 16, etc. This means that the CUDA program will first be invoked 'iterations' times with block size 1024 and grid size 8; after that, it will be invoked with block size 512 grid size 16; and so on. As one can see, the total number of threads (the grid size × the block size) is the same in the subsequent steps.

To run the script one may need to make it executable as, e.g., 'chmod +x time_order_1.sh'. Then it can simply be invoked as './time_order_1.sh', without any arguments.

The script generates two kinds of output files. The names of the files indicate parameters.
	1. The first kind uses names of the following form: 'meas_thread_1024_block_16_order_1.dat'. The block size is 1024, the grid size is 16, and the order of the L norm is 1 in this example. Such a file contains the standard output of the corresponding CUDA program along with that of the built-in 'time' command.
	2. For the second kind, names take the following form: 'proc_thread_1024_block_16_order_1.dat'. The parameters are the same. Such a file contains the runtimes of the corresponding CUDA program in seconds.

*calc_Statistics.m*
	
The Octave script 'calc_Statistics.m' can postprocess the 'proc_thread_...' files to generate the statistics of the execution times. The following variables must be defined in this script:
 1. 'subdir' should be defined as empty if the Octave script can be found in the same directory as the 'proc_thread_...' files. In case the Octave script is in one directory and the 'proc_thread_...' files can be found in a subfolder of that directory, 'subdir' should be the name of this subfolder.
 2. 'order' refers to the order of the L norm for which the calculations were performed.
 3. 'threads' must be the same as in the shell script 'time_order_1.sh'.
 4. Similarly for 'blocks'.

The script writes output to a file named like 'grid_8192_order_1.txt', where 8192 stands for the total number of threads and 1 for the order of the L norm. The output files have 8 columns. The columns stand for
 1. the grid size (number of blocks),
 2. the block size (number of threads in one block),
and
 3. the average,
 4. the minimum,
 5. the first quartile,
 6. the median,
 7. the third quartile,
 8. and the maximum
of the runtimes corresponding to the given grid and block size.

*runTime.gnu*

The output of 'calc_Statistics.m' can be plotted by the 'runTime.gnu' gnuplot script. It will generate a file 'runTime.eps' which depicts the runtime as a function of the block size. The candlesticks depict the minimal runtime, the first and third quartile, and the maximum of the runtimes. The total number of threads is the same for candlesticks plotted with the same color.
