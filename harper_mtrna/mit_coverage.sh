set -o pipefail  # trace ERR through pipes
set -o errtrace  # trace ERR through 'time command' and other functions
set -o nounset   ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

for NAME in D1 D2 D3 T1 T2 T3 ; do
    FILE=$NAME.bam
    echo $FILE
    if  [ ! -e $FILE.bai ] ; then
        sambamba index -p $FILE 
    fi
    if [ ! -e $NAME"_mt-tm.bam" ] ; then
        sambamba slice $FILE chrM:3932-4744 >| $NAME"_mt-tm.bam"
        sambamba index -p $NAME"_mt-tm.bam"
    fi
    if [ ! -e $NAME"_mt-tl.bam" ] ; then
        sambamba slice $FILE chrM:8230-8450 >| $NAME"_mt-tl.bam"
        sambamba index -p $NAME"_mt-tl.bam"
    fi
    if [ ! -e $NAME"_coverage.dat" ] ; then
        bamtools coverage -in $NAME"_mt-tm.bam" -out $NAME"_coverage.dat"
    fi
    if [ ! -e $NAME"_coverage_L.dat" ] ; then
        bamtools coverage -in $NAME"_mt-tl.bam" -out $NAME"_coverage_L.dat"
    fi
done


