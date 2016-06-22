#!/bin/sh
#SBATCH --job-name conduction
#SBATCH -N 1
#SBATCH -t 00:05:00

module load OpenMPI mpi4py numpy matplotlib
mpirun -np 8 python /home/tab10/code/conduction/conduction/mpi_run.py