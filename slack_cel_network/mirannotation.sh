#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

SB_HOME=~/repos/seqbuster/
SB_DB="$SB_HOME/modules/miraligner/DB"
RNAB_DB="~/soft/srnabench-db/sRNAbenchDB"
RNAB="~/soft/srnabench"

echo "download miRBase files"
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz

gunzip hairpin.fa.gz miRNA.str.gz

awk '{if ($0~/>/){print name"\t"seq;name=$1;seq=""}else{seq=seq$0}}' hairpin.fa | grep "^>hsa" | sed "s/U/T/g" | sed 's/\t/\n/' > hairpin.hsa.fa


java -jar  "$SB_HOME"/modules/miraligner/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i $NAME -db $SB_DB -o $NAME.out -pre



