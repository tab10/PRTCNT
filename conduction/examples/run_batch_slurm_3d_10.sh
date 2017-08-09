#!/usr/bin/env bash
rand_vals=( 1794
3609
5403
7212
8992
10801
12604
14429
16228
18060
19845 )
hv_vals=( 1813
3640
5453
7266
9080
10906
12719
14533
16376
18190 )
echo "SCRIPT DESIGNED FOR TUBE LENGTH 10 3D"
if [ "$1" = "random" ]
then
        for i in ${rand_vals[@]}; do
        cp slurm_template.sh temp_$i.sh
        sed -i 's/X/'$i'/g' temp_$i.sh
        sbatch temp_$i.sh
        rm temp_$i.sh
        done
elif [ "$1" = "horizontal" ] || [ "$1" = "vertical" ]
then
        for i in ${hv_vals[@]}; do
        cp slurm_template.sh temp_$i.sh
        sed -i 's/X/'$i'/g' temp_$i.sh
        sbatch temp_$i.sh
        rm temp_$i.sh
        done
fi