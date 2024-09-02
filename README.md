# damagedmuscle
Analysis code for the damaged muscle exosomes project, in collaboration with Hongwen's lab.

# Quick navigation
* [notebooks](notebooks): Main files containing analysis code. Start here. 
* [output](output): Intermediate output files from notebooks (RDS, csv, etc). Large RDS files are stored in server under `home/qianhui/damagedmuscle_project/output`.
* [scripts](scripts): Scripts. Contains R function scripts and bash scripts for FastQC, cutadapt, bwa and featureCount. 
* [data](data): Raw, aligned and processed data (RNA-seq raw counts, bwa and featureCount outputs). Stored in server under `home/qianhui/damagedmuscle_project/data`.

# A running list of issues

1. Initial FastQC shows read counts of varying lengths (14 - 41 bp). Pre-trimming done by Illumina? 
2. `Muscle_exosome_1` shows unusually low read counts (half that of other samples). Why?

# Experimental design overview

- 4 wells of C2C12 mouse myoblast cell lines were differentiated into myocytes using the (insert protocol here) protocol. 2 wells were treated with hydrogen peroxide. Exosomes were collected and sent for sequencing. 

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
        --outdir /home/qianhui/damagedmuscle_project/data/02_data_processed 
```
**Link:** [scripts/preprocess/00_fastqc.sh](scripts/preprocess/00_fastqc.sh)

**Discussion:**

As expected, sequences are all between 15-45 bp, with a peak at 22bp. See [data/02_data_processed/00_fastqc/](data/02_data_processed/00_fastqc/). This aligns with our expectations of miRNA lengths. 

However, this also implies that some trimming has been done by Illumina - not all reads are at 45bp. 

We check again that FastQC doesn't automatically trim adapters (it shouldn't) by running it without specifying the adapter sequences: 

```
fastqc  *.fq.gz \
        --threads 4 \
        --outdir /home/qianhui/damagedmuscle_project/data/02_data_processed/00_fastqc_noadapters

```

Yup they look the same - there's definitely pre-trimming done, as fastqc doesn't even pick up any adapter contamination. See [data/02_data_processed/00_fastqc_noadapters/](data/02_data_processed/00_fastqc/noadapters)

### 2. Adapter trimming with `cutadapt`

For tools like `STAR` - soft-clipping happens and so often having small portions of adapter sequences at the 3' end of the read doesn't affect alignments. However, when using something like `BWA` this can lead to poor mapping rates. So lets use `cutadapt` to trim adapter sequences from the 3' end of the read:

```
# Trim 5' and 3' adapters with non-internal matching

for i in damaged_muscle_exosome_1 damaged_muscle_exosome_2 \
         muscle_exosome_1 muscle_exosome_2 \

do
    cutadapt -g XGAACGACATGGCTACGATCCGACTT -a AGTCGGAGGCCAAGCGGTCTTAGGAAGACAAX -o /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_trimmed.fq.gz /home/qianhui/damagedmuscle_project/data/01_data_raw/${i}.fq.gz         --info-file /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_info_file.tsv >> /home/qianhui/damagedmuscle_project/data/02_data_processed/01_cutadapt/${i}_report.txt

done

```
Link: [scripts/preprocess/01_cutadapt.sh](scripts/preprocess/01_cutadapt.sh)

**Discussion:**

Cutadapt shows little to no adapter contamination - see [data/02_data_processed/01_cutadapt/](data/02_data_processed/01_cutadapt/).

Reads with adapters constitute 2.5 - 5.3% of total reads, and most trimmed reads are between 3 to 5 base pairs, which are just partial matches on the 3' end. We do not see any large 31bp adapters being trimmed. This means that our RNA sequences are of a good quality, and that there is little to no adapter contamination (as expected). 

However, there is a possibility of pre-trimming done by the company, which we need to clarify in the next meeting. This is further compounded by the fact that read counts are not of a uniform length even prior to adapter trimming - see fastqc reports in `preprocess/00_fastqc`.

All this is good. Since we have no adapter contamination, we can and proceed to `bwa` with the raw `fq.gz` files without any `cutadapt` trimming. 


### 3. Quality control of trimmed reads with `FastQC`

We use `FastQC` to check read length after trimming: 

```bash
fastqc  *.fq.gz \
        --threads 4 \
        --outdir /home/qianhui/damagedmuscle_project/data/02_data_processed/02_fastqc_trimmed 
```
Link: [scripts/preprocess/02_fastqc_trimmed.sh](scripts/preprocess/02_fastqc_trimmed.sh)


**Discussion**: 

Before trimming, reads had read lengths between 15 - 45 bps. See [data/02_data_processed/00_fastqc/](data/02_data_processed/00_fastqc/).

After cutadapt trimming, reads and lengths between 0 - 45 bps, with a peak around 23 bps. However, the general pattern remains the same, which means that some 31 bp fragments were adapters that got cut? See [data/02_data_processed/00_fastqc/](data/02_data_processed/02_fastqc_trimmed/). 

This further supports our conclusion above that cutadapt trimming is unnecessary for this dataset. We proceed to align the raw reads with `bwa`: 

### 4. Align reads with `bwa`

```
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
              -o /home/qianhui/damagedmuscle_project/data/02_data_processed/04_featurecounts/bwa_counts.txt \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_1_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_2_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_1_sorted.bam \
              /home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_2_sorted.bam \
```
Link: [scripts/preprocess/04_featurecounts.sh](scripts/preprocess/04_featurecounts.sh)


# Creating .bai files for viewing on IGV

We perform a visual check using IGV to ensure that genomes are all aligned properly. To do so, we first generate the `.bai` files for all the `.bam` files using `samtools`:

```
# Run this in the 03_bwa folder where the .bam files are:

samtools index *.bam

```

# Output

```
Status	/home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_1_sorted.bam	/home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/damaged_muscle_exosome_2_sorted.bam	/home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_1_sorted.bam	/home/qianhui/damagedmuscle_project/data/02_data_processed/03_bwa/muscle_exosome_2_sorted.bam
Assigned	6901394	7350345	4577668	7638663
Unassigned_Unmapped	1673410	1564508	1018240	1272528
Unassigned_Read_Type	0	0	0	0
Unassigned_Singleton	0	0	0	0
Unassigned_MappingQuality	0	0	0	0
Unassigned_Chimera	0	0	0	0
Unassigned_FragmentLength	0	0	0	0
Unassigned_Duplicate	0	0	0	0
Unassigned_MultiMapping	0	0	0	0
Unassigned_Secondary	0	0	0	0
Unassigned_NonSplit	0	0	0	0
Unassigned_NoFeatures	19011267	18905412	23753141	19211198
Unassigned_Overlapping_Length	0	0	0	0
Unassigned_Ambiguity	307420	221581	149658	180944
```
Link: [data/02_data_processed/04_featurecounts/bwa_counts.txt.summary](data/02_data_processed/04_featurecounts/bwa_counts.txt.summary)



# Additional information 

## Alternative pipeline with STAR 

As miRNAs are not spliced, we should not need to use a splice-aware aligner like STAR. However, just to be sure, we ran a STAR alignment and compared the results with BWA. 

We found that STAR was not ideal for miRNAs, as it is unable to handle multimapping reads (which are typical for miRNAs since they are very short sequences).
To find the STAR pipeline and results, see the following: 

- STAR alignment: [scripts/preprocess/alt_01_STAR.sh](scripts/preprocess/alt_01_STAR.sh)
- featureCount for STAR: [scripts/preprocess/alt_02_featurecounts.sh](scripts/preprocess/alt_02_featureCounts.sh)
- featureCount output for STAR: [data/02_data_processed/alt_02_featurecounts/star_counts.txt.summary](data/02_data_processed/alt_02_featurecounts/star_counts.txt.summary)

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

**Mouse annotations**: Mus_musculus.GRCm38.102.chr.gtf, mm10.fa
