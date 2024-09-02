
# Run featureCounts on STAR output, save all to one giant results file

featureCounts -T 10 \
              -F GTF \
              -t miRNA \
              -g ID \
              -a /home/shared/resources/mirna/mmu.gff3 \
              -o /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_02_featurecounts/star_counts.txt \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_01_STAR/damaged_muscle_exosome_1Aligned.sortedByCoord.out.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_01_STAR/damaged_muscle_exosome_2Aligned.sortedByCoord.out.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_01_STAR/muscle_exosome_1Aligned.sortedByCoord.out.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_01_STAR/muscle_exosome_2Aligned.sortedByCoord.out.bam

