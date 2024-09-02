Main folder for analysis code used in this project. 

## File types
Each analysis file comes in three flavours, and are used for various purposes. 

* To view on GitHub: `.md`
* To reproduce analysis: `.ipynb`, `.Rmd`. Both files are paired via [Jupytext](https://jupytext.readthedocs.io/en/latest/).

I use [Quarto](https://quarto.org) to render `.ipynb` files, generating the github-flavored markdown (gfm) `md` file. 

## Useful links 
* Preprocessing scripts for `STAR` and `RSEM`: [scripts/preprocess](../scripts/preprocess/)
* Outputs: [output/](../output/). All `RDS` files are on the server. 
* Aligned count matrix: [data/02_data_processed/04_featurecounts/bwa_counts.txt](../data/02_data_processed/04_featurecounts/bwa_counts.txt)

        
