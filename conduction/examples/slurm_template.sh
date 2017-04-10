#!/bin/sh -login
#SBATCH --job-name conduction
#SBATCH --partition=normal
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=20
#SBATCH -t 06:00:00

module load mpi4py/1.3.1-intel-2015b-Python-2.7.10
module load matplotlib/1.5.0-intel-2016a-Python-2.7.10
module load numpy/1.10.4-intel-2016a-Python-2.7.10

for i in {1..1}; do
mpirun python /home/tab10/code/conduction/conduction/mpi_run.py --dim 3 --orientation horizontal --num_tubes X --tube_length 15 --gen_plots True --rules_test True
done


mpirun -np 2 python mpi_run.py --dim 3 --orientation horizontal --num_tubes 5000 --tube_length 15 --gen_plots True --rules_test True --timesteps 35000 --num_walkers 100