# Differences in miRNA cargoes of EVs secreted from normal and damaged muscle

Sarcopenic muscles exhibit reduced EV secretion and altered EV cargoes, promoting cancer progression. miRNA sequencing of C2C12 myotubes with and without hydrogen peroxide treatment identified miR-7a-5p, a tumor-suppressor microRNA highly enriched in healthy muscle EVs but diminished in damaged muscle EVs. 

# Quick navigation
* [notebooks](notebooks): Main files containing analysis code. Start here. 
* [output](output): Intermediate output files from notebooks (RDS, csv, etc). 
* [scripts](scripts): Scripts. Contains R function scripts and bash scripts for FastQC, cutadapt, bwa and featureCount. 
* [data](data): Raw, aligned and processed data (RNA-seq raw counts, bwa and featureCount outputs). 

# Experimental design overview

- C2C12 mouse myoblast cell lines were differentiated into myocytes. Two replicates were treated with hydrogen peroxide. 
- Exosomes were collected and sent for sequencing. 

#  Pipeline overview

## Summary
1. Quality control with `FastQC`
2. Adapter trimming with `cutadapt`
3. Quality control of trimmed reads with `FastQC`
4. Align reads with `bwa`
5. Count reads with `featureCount`

## Details

### 1. Quality control with `FastQC`

```
fastqc  *.fq.gz \
        --adapters adapters.txt \
        --threads 4 \
        --outdir data/02_data_processed 
```
**Link:** [scripts/preprocess/00_fastqc.sh](scripts/preprocess/00_fastqc.sh)

### 2. Adapter trimming with `cutadapt`

```
# Trim 5' and 3' adapters with non-internal matching

for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    cutadapt -g XGAACGACATGGCTACGATCCGACTT -a AGTCGGAGGCCAAGCGGTCTTAGGAAGACAAX -o damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_trimmed.fq.gz damagedmuscle_project/data/01_data_raw/${i}.fq.gz         --info-file damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_info_file.tsv >> damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_report.txt

done

```
Link: [scripts/preprocess/01_cutadapt.sh](scripts/preprocess/01_cutadapt.sh)

### 3. Quality control of trimmed reads with `FastQC`

We use `FastQC` to check read length after trimming: 

```bash
fastqc  *.fq.gz \
        --threads 4 \
        --outdir data/02_data_processed/02_fastqc_trimmed 
```
Link: [scripts/preprocess/02_fastqc_trimmed.sh](scripts/preprocess/02_fastqc_trimmed.sh)

### 4. Align reads with `bwa`

```
# Alignment 
for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    bwa aln -t 10 \
    /home/shared/genomes/mm10/BWAIndex/mm10.fa \
    damagedmuscle_project/data/01_data_raw/${i}.fq.gz > \
    damagedmuscle_project/data/02_data_processed/03_bwa/${i}.sai

done

# Output to bam file 
for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    bwa samse /home/shared/genomes/mm10/BWAIndex/mm10.fa \
             damagedmuscle_project/data/02_data_processed/03_bwa/${i}.sai \
             damagedmuscle_project/data/01_data_raw/${i}.fq.gz | \
             samtools sort -o damagedmuscle_project/data/02_data_processed/03_bwa/${i}_sorted.bam 
done
```
Link: [scripts/preprocess/03_bwa.sh](scripts/preprocess/03_bwa.sh)

### 5. Count reads with `featureCount`

```bash
# Run featureCounts on bwa output, save all to one giant results file

featureCounts -T 10 \
              -F GTF \
              -t miRNA \
              -g ID \
              -a /home/shared/resources/mirna/mmu.gff3 \
              -o damagedmuscle_project/data/02_data_processed/04_featurecounts/bwa_counts.txt \
              damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_1_sorted.bam \
              damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_2_sorted.bam \
              damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_1_sorted.bam \
              damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_2_sorted.bam \
```
Link: [scripts/preprocess/04_featurecounts.sh](scripts/preprocess/04_featurecounts.sh)


# Creating .bai files for viewing on IGV

We perform a visual check using IGV to ensure that genomes are all aligned properly. To do so, we first generate the `.bai` files for all the `.bam` files using `samtools`:

```
# Run this in the 03_bwa folder where the .bam files are:

samtools index *.bam

```

# Additional information 

## Adapter sequences

```
# List of adapters used 

3prime	AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
5prime  GAACGACATGGCTACGATCCGACTT
```
Link: [data/01_data_raw/adapters.txt](data/01_data_raw/adapters.txt)
## Version info

**FastQC**: v0.11.8

**Cutadapt**: v2.4

**BWA**: 0.7.17-r1198-dirty

**samtools**: v1.9

**STAR**: v2.7.1a

**featureCounts**: v2.0.1

**Mouse annotations**: mm10.fa
