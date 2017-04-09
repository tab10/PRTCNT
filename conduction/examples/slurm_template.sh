#!/bin/sh
#SBATCH --job-name conduction
#SBATCH -N 10
#SBATCH -t 24:00:00

for i in {1..1}; do
mpirun -np 80 python /home/tab10/code/conduction/conduction/mpi_run.py --dim 3 --orientation horizontal --num_tubes X --tube_length 15 --gen_plots True --rules_test True
done


mpirun -np 2 python mpi_run.py --dim 3 --orientation horizontal --num_tubes 5000 --tube_length 15 --gen_plots True --rules_test True --timesteps 35000 --num_walkers 100