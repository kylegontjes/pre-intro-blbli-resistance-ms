for i in $(ls *sbat | grep -v model)
do
echo $i
sbatch $i
done