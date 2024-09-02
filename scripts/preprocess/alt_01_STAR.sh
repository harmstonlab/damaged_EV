for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
  STAR  --genomeDir /home/shared/genomes/mm10/StarIndex/ensembl102 \
        --sjdbGTFfile /home/shared/genomes/mm10/StarIndex/ensembl102/Mus_musculus.GRCm38.102.chr.gtf \
        --readFilesIn /home/qianhui/damagedmuscle_project/data/01_data_raw/${i}.fq.gz \
        --runThreadN 20 \
        --twopassMode Basic \
        --outWigType bedGraph \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat \
        --runDirPerm All_RWX \
        --outFileNamePrefix /home/qianhui/damagedmuscle_project/data/02_data_processed/alt_01_STAR/${i}
done