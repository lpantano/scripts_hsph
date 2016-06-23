cd /n/data1/cores/bcbio/walsh_asd/batch_1/batch1/final

for VCF in `ls */*vcf.gz` ; do
    echo $VCF
    NAME=`basename $VCF | sed 's/.vcf.gz//'`
    DB="2016-04-12_batch1/"$NAME.db
    bsub -R "rusage[mem=3000]" -W 5:00 -n 2 -q mcore /groups/bcbio/bcbio/anaconda/bin/gemini  annotate -f $VCF -a extract  -c af -t float -e AF -o first $DB
done

cd 2016-04-12_batch1

for DB in `ls *db`;  do
    /groups/bcbio/bcbio/anaconda/bin/gemini query --header -q "select *,(gts).(*) from variants where impact_severity='HIGH' and af<0.05 and in_exac=0" $DB > filtered/$DB.aaf_5.no_exac.high.tsv
done
