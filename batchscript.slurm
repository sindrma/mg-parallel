#!/bin/bash
# Account name
#SBATCH -A TG-CCR100030
# Job name,  will be displayed on the showq command
#SBATCH -J MG_SER
# Filename for standard output 
# At end of job, it is in directory from which sbatch was invoked
#SBATCH -o MG_SER.o%j
#SBATCH -n 2 
<<<<<<< Updated upstream
#SBATCH -p development
#  The requested wall clock job time limit in HH:MM:SS
#  Your job will end when it exceeds this time limit
#SBATCH -t 00:10:00 
=======
#SBATCH -p gpu
#  The requested wall clock job time limit in HH:MM:SS
#  Your job will end when it exceeds this time limit
#SBATCH -t 00:5:00 
>>>>>>> Stashed changes

export OMP_NUM_THREADS=8

# Run
echo ">>> Job Starts"

<<<<<<< Updated upstream
ibrun -np 2 ./mg
=======
ibrun -np 1 ./mg -n 256
>>>>>>> Stashed changes

echo ">>> Job Ends"

date

# Print out the processor information
<<<<<<< Updated upstream
cat /proc/cpuinfo


# Print out the environment
printenv
=======
# sbatch batchscript.slurm
# cat /proc/cpuinfo
>>>>>>> Stashed changes
