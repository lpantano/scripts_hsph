#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

SB_HOME=~/soft/seqbuster/
SB_DB="./"
DIR=/n/data1/cores/bcbio/tina_exosomes/work
ADAPTER=AGATCGGAAGAGCACACGTCT
PATH=/home/lp113/soft/bcbiome/anaconda/bin:$PATH

cd $DIR 

if [ ! -e DB ] ; then
    mkdir -p DB
    cd DB
    echo "download miRBase files"
    wget -q ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
    wget -q ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz
    gunzip hairpin.fa.gz miRNA.str.gz
    awk '{if ($0~/>/){print name"\t"seq;name=$1;seq=""}else{seq=seq$0}}' hairpin.fa | sed "s/U/T/g" | sed 's/\t/\n/' > hairpin.hsa.fa
    cd ..
fi

rm -rf ../mirna.files ../seqcluster.files


while read file cell ; do
 NAME=$cell
 fastq=$file
 if [ ! -e  $NAME.fastq ] ; then
     cutadapt --adapter=$ADAPTER --minimum-length=8 --untrimmed-output=$NAME.notfound.fastq -o $NAME.fastq -m 17 --overlap=8 ../raw/$file
 fi
 if [ ! -e $NAME.collapse.fastq ] ; then
     seqcluster collapse -f $NAME.fastq -o collapse
 fi
 #awk '{i=i+1; print ">seq_"i"_x"$2"\n"$1}' $NAME.ad > $NAME.fa
 if [ ! -e $NAME-ann.mirna ] ; then
    java -jar "$SB_HOME"/modules/miraligner/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i collapse/$NAME.fastq -db DB -o $NAME-ann -pre -freq -minl 17
 fi
 
 #echo $NAME-ann.mirna $NAME>>../mirna.files
 echo $NAME.fasq $NAME>>../seqcluster.files
done < $1 #smalrna.csv


