# Name of pipeline
create-weight-dbs

## Description
Generates PredictDB sqlite dbs using weights calculated by generate-prediction-models. These dbs are used for gene expression prediction by PrediXcan and for gene-level association by MetaXcan.

This takes all files in the input folder and creates dbs
## Input files

- `tissue-name.allBetas.txt.gz`: has the gene-rsid pairs with the corresponding weights (called beta here)
> gene    rsid    ref     alt     beta    alpha

  - `gene` is a unique identifier such as ensid or other
  - `rsid` is rs number of SNP. This will be changed to chr_position_A1_A2_build later
  - `alt` is the effect allele
  - `ref` is the non effect allele
  - `beta` has the weight for the gene-rsid pair
  - `alpha` describes model (for now numeric = mixing parameter in Elastic Net)
  - additional columns are ignored

- `tissue-name.allResults.txt`: has descriptiones of the model runs, CV performance (R2)
> gene    alpha   cvm     lambda.iteration        lambda.min      n.snps  R2      pval    genename

  - `gene` is a unique identifier such as ensid or other
  - `alpha` describes model (for now numeric = mixing parameter in Elastic Net)
  - `genename` is the HUGO name of the gene
  - `R2` is the cross validated correlation squared between predicted and observed expression

## Output files

- `tissue-name_alpha.db`: as many files as alpha values in the input betas/results files will be generated.

## Usage
> generate_sqlite_dbs.py

## Working Example
- Download [link]()
- Untar the file ()
- cd to the untarred folder
- Run the following line
> generate_sqlite_dbs.py

## Options
No options.
