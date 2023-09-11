# Acute Myeloid Leukemia Heatmap

Author: CCDL for ALSF - Adapted for this repository by Candace Savonen
Date: October 2021

_This analysis has been adapted from this [refine.bio-examples notebook](https://alexslemonade.github.io/refinebio-examples/03-rnaseq/clustering_rnaseq_01_heatmap.html)_

## Purpose of the Analysis

In this analysis, we will use the [acute myeloid leukemia sample dataset](https://www.refine.bio/experiments/SRP070849) from [Shih et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28193779/) and pre-processed by [refinebio](https://www.refine.bio/).

The dataset contains 19 samples obtained from 19 acute myeloid leukemia (AML) model mice, and it includes RNA-sequencing results for different types of AML under controlled treatment conditions.

## Setup Analysis Folders

To set up the analysis, follow these steps:

1. Create the "data" folder if it doesn't exist.
2. Define the file path to the "plots" directory and create it if it doesn't exist.
3. Define the file path to the "results" directory and create it if it doesn't exist.

## Clustering Heatmap - RNA-seq

### Install Libraries

Install the "pheatmap" library by Slowikowski et al., 2017 for clustering and creating a heatmap.

R if (!("pheatmap" %in% installed.packages())) { install.packages("pheatmap", update = FALSE) }


### Import and Set Up Data

R
Read in metadata TSV file

metadata <- readr::read_tsv(metadata_file)
Read in data TSV file

expression_df <- readr::read_tsv(data_file) %>% tibble::column_to_rownames("Gene")


### Choose Genes of Interest

R
Calculate the variance for each gene

variances <- apply(expression_df, 1, var)
Determine the upper quartile variance cutoff value

upper_var <- quantile(variances, 0.75)
Filter the data choosing only genes whose variances are in the upper quartile

df_by_var <- data.frame(expression_df) %>% dplyr::filter(variances > upper_var)


### Prepare Metadata for Annotation

R annotation_df <- metadata %>% dplyr::mutate( mutation = dplyr::case_when( startsWith(refinebio_title, "TET2") ~ "TET2", startsWith(refinebio_title, "IDH2") ~ "IDH2", startsWith(refinebio_title, "WT") ~ "WT", TRUE ~ "unknown" ) ) %>% dplyr::select( refinebio_accession_code, mutation, refinebio_treatment ) %>% tibble::column_to_rownames("refinebio_accession_code")


### Create Annotated Heatmap

R heatmap_annotated <- pheatmap( df_by_var, cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = annotation_df, show_rownames = TRUE, show_colnames = TRUE )


This README provides an overview of the Acute Myeloid Leukemia Heatmap analysis and steps to set up
