#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name conduction
#SBATCH --ntasks=200
#SBATCH --ntasks-per-node=20
#SBATCH -t 00:10:00

module load mpi4py/1.3.1-intel-2015b-Python-2.7.10
module load numpy/1.10.4-intel-2016b-Python-2.7.10
module load matplotlib/1.5.0-intel-2016a-Python-2.7.10
module load scipy/0.17.0-intel-2016a-Python-2.7.12

for i in {1..1}; do
mpirun python /home/tab10/code/conduction/conduction/mpi_run.py --dim 3 --orientation random --prob_m_cn 0.5 --num_tubes X --tube_length 15 --gen_plots True --num_walkers 210000 --timesteps 35000 --model kapitza
done
