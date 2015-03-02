#!/bin/bash
# Account name
#SBATCH -A TG-CCR100030
# Job name,  will be displayed on the showq command
#SBATCH -J MG_SER
# Filename for standard output 
# At end of job, it is in directory from which sbatch was invoked
#SBATCH -o MG_SER.o%j
#SBATCH -n 2 
#SBATCH -p development
#  The requested wall clock job time limit in HH:MM:SS
#  Your job will end when it exceeds this time limit
#SBATCH -t 00:10:00 

export OMP_NUM_THREADS=8

# Run
echo ">>> Job Starts"

ibrun -np 2 ./mg

echo ">>> Job Ends"

date

# Print out the processor information
cat /proc/cpuinfo


# Print out the environment
printenv
