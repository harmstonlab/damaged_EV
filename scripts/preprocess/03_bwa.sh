# Alignment 
for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    bwa aln -t 10 \
            /home/shared/genomes/mm10/BWAIndex/mm10.fa \
            /home/qianhui/damagedmuscle_project/data/01_data_raw/${i}.fq.gz > \
            /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/${i}.sai

done

# Output to bam file 
for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    bwa samse /home/shared/genomes/mm10/BWAIndex/mm10.fa \
             /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/${i}.sai \
             /home/qianhui/damagedmuscle_project/data/01_data_raw/${i}.fq.gz | \
             samtools sort -o /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/${i}_sorted.bam 
done