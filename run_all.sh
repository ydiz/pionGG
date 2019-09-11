cur_dir=$PWD
# ensembles=(Pion_24ID Pion_32ID Pion_32IDF Pion_48I Pion_24ID_disc)
ensembles=(Pion_48I)
targets=(imag_analytic real form_factor)
declare -A file_p3s
file_p3s=([Pion_24ID]="$cur_dir/integral/24ID/data.txt"
[Pion_24ID_disc]="$cur_dir/integral/24ID/data.txt"
[Pion_32ID]="$cur_dir/integral/32ID/data.txt"
[Pion_32IDF]="$cur_dir/integral/32IDF/data.txt"
[Pion_48I]="$cur_dir/integral/48I/data.txt"
)

declare -A time_cutoff_ends
time_cutoff_ends=([Pion_24ID]=16
[Pion_24ID_disc]=16
[Pion_32ID]=16
[Pion_32IDF]=16
[Pion_48I]=24
)


program=/gpfs/mira-home/yidizhao/cooley/pionGG/jackknife


mkdir -p submit_scripts

for ensemble in "${ensembles[@]}"; do
  for target in "${targets[@]}"; do
    echo $ensemble $target
    mkdir -p ./rst/$ensemble

    outfile=$cur_dir/rst/$ensemble/$target.txt

    submit_script=./submit_scripts/"$ensemble"_"$target".sh
    sed -e "s|@file_p3@|${file_p3s[$ensemble]}|" -e "s|@ensemble@|$ensemble|" -e "s|@target@|$target|" -e "s|@outfile@|$outfile|" submit_template.sh -e "s|@time_cutoff_end@|${time_cutoff_ends[$ensemble]}|" > $submit_script
    chmod a+x $submit_script
    qsub -n 1 -t 00:10:00 -A CSC249ADSE03 $submit_script

    sleep 0.5
  done
done
