#!/bin/sh
#SBATCH --job-name conduction
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH --exclusive

mpirun -np 8 python /home/tab10/code/conduction/conduction/mpi_run.py --dim 3 --orientation random --num_tubes --tube_length