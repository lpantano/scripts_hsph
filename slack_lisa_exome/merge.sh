bcftools  merge -f PASS vardict-java/bN*effects.vcf.gz -o vardict-java/samples.vcf.gz -O z
bcftools query -i 'FILTER==PASS' -f '[\t%AF][\t%DP][\t%SAMPLE]\n' vardict-java/samples.vcf.gz > ../final/2015-12-07_tumor_mm10/samples_af_depth.tsv
