#!/usr/bin/env bash
rand_vals=( 8
15
23
30
39
46
53
60
70
75
82
89
99
108
116 )
hv_vals=( 8
15
23
30
38
46
53
61
69
76
84
91
99
107 )
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