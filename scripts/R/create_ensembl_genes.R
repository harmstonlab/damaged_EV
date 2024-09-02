## Title: Create ensembl_genes.RDS
## Author: Qian Hui TAN
## Date modified: 27-Oct-2022

## Description: 
## Creates a dataframe of gene annotations. For every gene 
## in the mouse genome, we get the ensembl ID, gene biotype, chromosome, 
## ranges, strand, gene name, etc. 



## -- Setup -- ## 
# Load libraries
library(biomaRt)
library(GenomicFeatures)

## -- Creating transcript annotations -- ## 

# Create transcript annotations from the mouse genome
mm.gtf.db <- makeTxDbFromGFF(
  "../../data/refseq/Mus_musculus.GRCm38.102.chr.gtf",
  format = "gtf")

# Extract genes from the transcript annotation.
ensembl.genes = genes(mm.gtf.db)



## -- Retrieving biomaRt annotations -- ##

# For each ensembl gene in the drosophila genome, we retrieve additional 
# information such as the gene biotype, gene name, entrez id. These will be
# useful for downstream analysis steps. 

# Set fly 
mouse = useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "https://asia.ensembl.org",
  dataset = "mmusculus_gene_ensembl",
  version = "102"
)

# For each gene, retrieve the following annotations from biomaRt
bm.annotations = getBM(
  attributes = c("ensembl_gene_id",# FBgn, highly specific identifier"
                 "external_gene_name", # normal gene name
                 "gene_biotype",
                 "go_id", # go term accession
                 "name_1006", # go term name,
                 "entrezgene_id",# entrezid for pathway
                 "description"
  ),
  mart = mouse,
  filters = "ensembl_gene_id",
  values = ensembl.genes$gene_id,
  uniqueRows = TRUE
)

# Checking bm.annotations
head(bm.annotations)

## -- Creating the object -- ##
# Create ensembl.genes object
ensembl.genes$ensembl_gene_id = bm.annotations$ensembl_gene_id[
  match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]

ensembl.genes$external_gene_name = bm.annotations$external_gene_name[
  match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]

ensembl.genes$gene_biotype = bm.annotations$gene_biotype[
  match(ensembl.genes$gene_id,bm.annotations$ensembl_gene_id)]
ensembl.genes$go_id = bm.annotations$go_id[
  match(ensembl.genes$gene_id,bm.annotations$ensembl_gene_id)]
ensembl.genes$name_1006 = bm.annotations$name_1006[
  match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]
ensembl.genes$entrezgene_id = bm.annotations$entrezgene_id[
  match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]

head(ensembl.genes)

# Save the output
saveRDS(ensembl.genes, file = "../../output/00_ensembl_genes.RDS")
