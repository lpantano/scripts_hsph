#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

DIR=/n/data1/cores/bcbio/tina_exosomes/work
PATH=/home/lp113/soft/bcbiome/anaconda/bin:$PATH

cd $DIR 

while read file cell ; do
 NAME=$cell
 fastq=$file
 seqcluster prepare -c seqcluster.files -o res -l 17 --minc 1 --maxl 50
 STAR 
 seqcluster cluster -m res/seqs.ma -a res/seqs.bam -b $DB/tRNA.bed,$DB/mirna.bed
done
