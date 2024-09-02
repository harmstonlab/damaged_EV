#!/usr/bin/bash

for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    cutadapt -g XGAACGACATGGCTACGATCCGACTT -a AGTCGGAGGCCAAGCGGTCTTAGGAAGACAAX -o /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_trimmed.fq.gz /home/qianhui/damagedmuscle_project/data/01_data_raw/${i}.fq.gz         --info-file /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_info_file.tsv >> /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_report.txt

done