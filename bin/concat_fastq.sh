#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

DIR=/n/data1/cores/bcbio/celegans_network/celegans_network_raw/celegans_mrna

cd $DIR

for FOLDER in `ls`; do
 cd $FOLDER
 echo $FOLDER
 if [ ! -e $FOLDER.fastq.gz ] ; then
  zcat */*fastq.gz | gzip > $FOLDER.fastq.gz
 fi
 cd ..
done


