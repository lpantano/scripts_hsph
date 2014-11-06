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
DB=/n/data1/cores/bcbio/tina_exosomes/data

cd $DIR 

seqcluster prepare -c ../seqcluster.files -o res -l 17 --minc 1 --maxl 50
STAR --genomeDir  /groups/bcbio/bcbio/genomes/Hsapiens/hg19/star --readFilesIn res/seqs.fastq --outFilterMultimapNmax 5000 --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate
 
seqcluster cluster -m res/seqs.ma -a Aligned.sortedByCoord.out.bam -b $DB/tRNA.bed,$DB/miRNA.bed,$DB/refGene.bed,$DB/snoRNA.bed,$DB/rmsk.bed,$DB/wgRNA.bed -o res/cluster

samtools index Aligned.sortedByCoord.out.bam

seqcluster stats -j res/cluster/seqcluster.json -m res/seqs.ma -a Aligned.sortedByCoord.out.bam -o res/stats 
