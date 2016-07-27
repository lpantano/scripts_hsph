Rscript -e 'library(rmarkdown);library(myRfunctions);ORIGINAL=TRUE;render(meta_proteomic_original_norm.Rmd)'
Rscript -e 'library(rmarkdown);library(myRfunctions);ORIGINAL=FALSE;render(meta_proteomic.Rmd)'
