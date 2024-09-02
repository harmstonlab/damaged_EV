# Run featureCounts on STAR output, save all to one giant results file

featureCounts -T 10 \
              -F GTF \
              -t miRNA \
              -g ID \
              -a /home/shared/resources/mirna/mmu.gff3 \
              -o /home/qianhui/damagedmuscle_project/data/02_data_processed/04_featurecounts/bwa_counts.txt \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_1_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_2_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_1_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_2_sorted.bam \


