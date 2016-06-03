ROOT="/orch/scratch/vishal_mirna_kidney/publish"
mkdir -p $ROOT/results/FA-model
mkdir -p $ROOT/results/UUO-model

# Run FA rnaseq
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/FA-model/rnaseq/final/2016-02-05_mrna_bio";library(rmarkdown);render(FA-model/mrna-summary.Rmd)'
mkdir -p $ROOT/results/FA-model/rnaseq
cp -r $ROOT/FA-model/rnaseq/final/2016-02-05_mrna_bio/files_publish/* $ROOT/results/FA-model/rnaseq/.

# Run FA srnaseq
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/FA-model/srnaseq/final";library(rmarkdown);render(FA-model/srna-summary.Rmd)'

# Run seqcluster
if [ ! -e /orch/scratch/vishal_mirna_kidney/publish/FA-model/srnaseq/final/files_publish/pairs ] ; then
    /usr/local/conda/bin/seqcluster target --input /orch/scratch/vishal_mirna_kidney/publish/FA-model/srnaseq/final/files_publish/mirna_sing.txt \
--sps mmu -o /orch/scratch/vishal_mirna_kidney/publish/FA-model/srnaseq/final/files_publish/pairs \
--annotation /orch/scratch/vishal_mirna_kidney/publish/srna-data
fi

# Run FA miRNA-mRNA interactome
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/FA-model";library(rmarkdown);render("FA-model/interactome.Rmd")'
mkdir -p $ROOT/results/FA-model/srnaseq
cp -r $ROOT/FA-model/srnaseq/final/files_publish/* $ROOT/results/FA-model/srnaseq/.

# Run FA 3D omics
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/FA-model";library(rmarkdown);render("FA-model/omics.Rmd")'
mkdir -p $ROOT/results/FA-model/protein
cp -r $ROOT/FA-model/protein/files_publish/* $ROOT/results/FA-model/protein/.
mkdir -p $ROOT/results/FA-model/omics
cp -r $ROOT/FA-model/omics/files_publish/* $ROOT/results/FA-model/omics/.

# Run UUO rnaseq
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/UUO-model/rnaseq/final/2016-02-03_mrna";library(rmarkdown);render("UUO-model/mrna-summary.Rmd")'
mkdir -p $ROOT/results/UUO-model/rnaseq
cp -r $ROOT/UUO-model/rnaseq/final/2016-02-03_mrna/files_publish/* $ROOT/results/UUO-model/rnaseq/.

# Run UUO srnaseq
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/UUO-model/srnaseq/final";library(rmarkdown);render("UUO-model/srna-summary.Rmd")'

# Run seqcluster
if [ ! -e /orch/scratch/vishal_mirna_kidney/publish/FA-model/srnaseq/final/files_publish/pairs ] ; then
    /usr/local/conda/bin/seqcluster target --input /orch/scratch/vishal_mirna_kidney/publish/UUO-model/srnaseq/final/files_publish/mirna_sing.txt \
--sps mmu -o /orch/scratch/vishal_mirna_kidney/publish/UUO-model/srnaseq/final/files_publish/pairs \
--annotation /orch/scratch/vishal_mirna_kidney/publish/srna-data
fi

# Run UUO miRNA-mRNA interactome
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/UUO-model";library(rmarkdown);render("UUO-model/interactome.Rmd")'
mkdir -p $ROOT/results/UUO-model/srnaseq
cp -r $ROOT/UUO-model/srnaseq/final/files_publish/* $ROOT/results/UUO-model/srnaseq/.

# Run UUO 3D omics
Rscript -e 'root_path="/orch/scratch/vishal_mirna_kidney/publish/UUO-model";library(rmarkdown);render("UUO-model/omics.Rmd")'
mkdir -p $ROOT/results/UUO-model/protein
cp -r $ROOT/UUO-model/protein/files_publish/* $ROOT/results/UUO-model/protein/.
mkdir -p $ROOT/results/UUO-model/omics
cp -r $ROOT/UUO-model/omics/files_publish/* $ROOT/results/UUO-model/omics/.

mkdir -p  $ROOT/results/code_and_html/.
cp -r * $ROOT/results/code_and_html/.