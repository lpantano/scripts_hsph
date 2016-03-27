# installing package. I am using dev version of everything, so be aware
devtools::install_github(c("Bioconductor-mirror/GenomicAlignments",
                           "Bioconductor-mirror/rtracklayer"))
BiocInstaller::biocLite(c("BSgenome","BSgenome.Hsapiens.UCSC.hg19"))
devtools::install_github("raerose01/deconstructSigs")

# load some data
root = file.path(system.file(package = "deconstructSigs", "data/"))

library(deconstructSigs)
library(useful)
data("sample.mut.ref")
head(sample.mut.ref)
tail(sample.mut.ref)

# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")
dim(sigs.input)
corner(sigs.input)


# how to compare to signatures
data(signatures.nature2013)
dim(signatures.nature2013)
corner(signatures.nature2013)

data(randomly.generated.tumors)
corner(randomly.generated.tumors)
dim(randomly.generated.tumors)
# Determine the signatures contributing an already normalized sample
test = whichSignatures(tumor.ref = randomly.generated.tumors, 
                       signatures.ref = signatures.nature2013, 
                       sample.id = 2)
test$weights


# test sample inside the signature
noise = data.frame(t(t(signatures.nature2013) + rnorm(96,0,0.001)),check.names = F)

test = whichSignatures(tumor.ref = noise, 
                       signatures.ref = signatures.nature2013, 
                       sample.id = "Signature.2")
test$weights
