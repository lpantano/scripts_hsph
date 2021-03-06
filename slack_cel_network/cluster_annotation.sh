#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

PATH=~/soft/bcbiome/anaconda/bin:~/soft/bcbiotools/bin:$PATH
GENOME="/groups/bcbio/bcbio/genomes/Celegans/WBcel235/bowtie/WBcel235"
DB="/n/data1/cores/bcbio/celegans_network/celegans_network_raw/data"
DIR=/n/data1/cores/bcbio/celegans_network/celegans_network_raw/celegans_smallrna
OUT=cluster
cd $DIR
rm -rf cluster/*
#if [ ! -e $OUT ] ; then
     mkdir -p  $OUT/pre
     mkdir -p  $OUT/res
     seqcluster prepare -c $1 -o $OUT/pre --minc 2 --minl 17 --maxl 42 
     bowtie -a --best --strata -m 5000 -f $GENOME $OUT/pre/seqs.fa -S  $OUT/pre/seqs.sam
     samtools view -Sbh $OUT/pre/seqs.sam | samtools sort -o /dev/stdin  tmp > $OUT/pre/seqs.sort.bam
     seqcluster cluster -a $OUT/pre/seqs.sort.bam -m $OUT/pre/seqs.ma -o $OUT/res  -b $DB/tRNA.bed,$DB/miRNA.bed,$DB/repeat.bed,$DB/rRNA.bed,$DB/snoRNA.bed,$DB/snRNA.bed,$DB/ncRNA.bed 
#fi

