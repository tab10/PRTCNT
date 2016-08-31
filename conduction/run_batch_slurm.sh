for i in `seq 50 50 500`; do
cp slurm_template.sh temp_$i.sh
sed -i 's/X/'$i'/g' temp_$i.sh
sbatch temp_$i.sh
rm temp_$i.sh
done
