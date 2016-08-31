#!/bin/sh
#SBATCH --job-name conduction
#SBATCH -N 5
#SBATCH -t 24:00:00

for i in {1..25}; do
mpirun -np 40 python /home/tab10/code/conduction/conduction/mpi_run.py --dim 3 --orientation horizontal --num_tubes X --tube_length 15 --gen_plots True
done
