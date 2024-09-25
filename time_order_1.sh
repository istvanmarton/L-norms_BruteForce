#!/bin/bash

# Compile the CUDA program
nvcc -o L_CUDA L_CUDA.cu

# Define variables
order=1
iterations=50  # Number of times to run the command
base_filename='meas'

threads=(1024 512 256 128 64 32 16 8) # Different thread values to use
blocks=(8 16 32 64 128 256 512 1024) # Different block values to use
# Check if both arrays have the same length
if [ ${#threads[@]} -ne ${#blocks[@]} ]; then
    echo "Error: The length of threads and blocks arrays must be the same."
    exit 1
fi
for i in "${!threads[@]}"; do
    t=${threads[$i]}
    b=${blocks[$i]}
    printf -v var "${base_filename}_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    printf -v var1 "proc_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    for j in $(seq 1 $iterations); do
        echo "Thread count: $t, Iteration: $j" >> "$var"
        { time ./L_CUDA "$b" "$t" mtx38.txt "$order"; } >> "$var" 2>&1
    done
    grep "real" "$var" | awk '{print $2}' >> "$var1"
done

threads=(1024 512 256 128 64 32 16 8) # Different thread values to use
blocks=(16 32 64 128 256 512 1024 2048) # Different block values to use
# Check if both arrays have the same length
if [ ${#threads[@]} -ne ${#blocks[@]} ]; then
    echo "Error: The length of threads and blocks arrays must be the same."
    exit 1
fi
# Create a base filename
# Run the executable multiple times in different threads and redirect all output to the file
order=1
for i in "${!threads[@]}"; do
    t=${threads[$i]}
    b=${blocks[$i]}
    printf -v var "${base_filename}_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    printf -v var1 "proc_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    for j in $(seq 1 $iterations); do
        echo "Thread count: $t, Iteration: $j" >> "$var"
        { time ./L_CUDA "$b" "$t" mtx38.txt "$order"; } >> "$var" 2>&1
    done
    grep "real" "$var" | awk '{print $2}' >> "$var1"
done

order=1
threads=(512 256 128 64 32 16 8) # Different thread values to use
blocks=(21 42 84 168 336 672 1344) # Different block values to use
# Check if both arrays have the same length
if [ ${#threads[@]} -ne ${#blocks[@]} ]; then
    echo "Error: The length of threads and blocks arrays must be the same."
    exit 1
fi

for i in "${!threads[@]}"; do
    t=${threads[$i]}
    b=${blocks[$i]}
    printf -v var "${base_filename}_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    printf -v var1 "proc_thread_%d_block_%d_order_%d.dat" "$t" "$b" "$order"
    for j in $(seq 1 $iterations); do
        echo "Thread count: $t, Iteration: $j" >> "$var"
        { time ./L_CUDA "$b" "$t" mtx38.txt "$order" n; } >> "$var" 2>&1
    done
    grep "real" "$var" | awk '{print $2}' >> "$var1"
done