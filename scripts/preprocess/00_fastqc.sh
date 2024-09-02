# FastQC script

fastqc  *.fq.gz \
        --adapters adapters.txt \
        --threads 4 \
        --outdir /home/qianhui/damagedmuscle_project/data/02_data_processed 