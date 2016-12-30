#!/usr/bin/env bash
rand_vals=( 951
1904
2863
3822
4757
5728
6660
7627
8638
9562
10519 )
hv_vals=( 950
1904
2854
3808
4758
5712
6662
7616
8575
9529 )
echo "SCRIPT DESIGNED FOR TUBE LENGTH 20 3D"
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