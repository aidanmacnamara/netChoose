# netChoose
An assessment of approaches for choosing proxies of genetic targets

This is the code and data used to generate the results for:

'Network and pathway expansion of genetic disease associations identifies successful drug targets'
https://www.biorxiv.org/content/10.1101/2020.04.22.051144v1

## Index

### scripts/

The R scripts for data collection, collation, and statistical analysis. The most relevant are listed below:

* pipeline_all_diseases.R: Reads in the data sources and produces a long data frame (genes * GWAS studies) that is used as input for the statistical analysis.

* pipeline_all_diseases_analysis.R: Runs the statistical analysis on the data frame produced by 'pipeline_all_diseases.R' (Stratified Fisher's Test).

* pipeline_all_diseases_final_report.R: Produces the figures for the main paper and the supplementary data.

### GWAS/

The data sources used to define High Confidences Genetic Hits, along with annotation files.

### final_report/

An R Markdown report of the results, together with the R objects used to generate the embedded figures.

### interactions/

The interaction data used to assess network approaches for target proxy identification.

### genesets/

Gene sets from MSigDB used to define 'pathway buckets' for Pascal enrichment.

### gwas_*

R lists showing the results of Pascal pathway enrichment.

