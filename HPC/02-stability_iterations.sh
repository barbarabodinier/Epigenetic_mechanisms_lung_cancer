ttd_list=(1 2 3 5)

for k in ${ttd_list[*]}
do
echo $k
sed "s/{ttd_input}/${k}/g" template_stability.sh > run1.sh
for j in $(seq 1 1 2)
do
echo $j
sed "s/{model_id_input}/${j}/g" run1.sh > run.sh
qsub run.sh
done
done

