#!/usr/bin/env bash
rand_vals=( 14
30
45
58
75
88
102
115
131
143
156
172
191
201
220 )
hv_vals=( 15
29
44
58
73
87
102
116
131
145
160
175
189
204 )
echo "SCRIPT DESIGNED FOR TUBE LENGTH 20 2D"
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