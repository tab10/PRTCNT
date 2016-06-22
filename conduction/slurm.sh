#!/bin/sh
#SBATCH --job-name conduction
#SBATCH -N 1
#SBATCH -t 00:00:30

module load OpenMPI
mpiexec -n 8 python mpi_run.py