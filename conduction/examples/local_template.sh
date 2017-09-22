#!/usr/bin/env bash

# for local mac/linux
mpirun -np 2 python /Users/tim/Documents/PycharmProjects/conduction/conduction/mpi_run.py --dim 2 --orientation random --num_tubes 10 --tube_length 10 --gen_plots True --rules_test False --timesteps 15000 --num_walkers 500 --model tunneling_w_vol --grid_size 30

# for local windows
#mpiexec -n 2 python C:\Users\Tim\Documents\PycharmProjects\conduction\conduction\mpi_run.py --dim 3 --orientation random --num_tubes 1 --tube_length 10 --gen_plots True --rules_test True --timesteps 15000 --num_walkers 250 --model kapitza --grid_size 20 --prob_m_cn 0.5 --disable_func False
