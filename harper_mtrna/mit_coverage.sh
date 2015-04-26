set -o pipefail  # trace ERR through pipes
set -o errtrace  # trace ERR through 'time command' and other functions
set -o nounset   ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

while read chr source type start end name; do 
    start=$(( start - 100  ))
    end=$(( end + 100  ))
    echo $chr $start $end $name
    for NAME in D1 D2 D3 T1 T2 T3 ; do
        FILE=$NAME.bam
        echo $FILE
        if  [ ! -e $FILE.bai ] ; then
            sambamba index -p $FILE 
        fi
        if [ ! -e slice/$NAME"_mt-"$name".bam" ] ; then
            sambamba slice $FILE chrM:$start-$end >| slice/$NAME"_mt-"$name".bam"
            sambamba index -p slice/$NAME"_mt-"$name".bam"
        fi
        if [ ! -e coverage_region/$NAME"_coverage_"$name".dat" ] ; then
            bamtools coverage -in slice/$NAME"_mt-"$name".bam" -out coverage_region/$NAME"_coverage_"$name".dat"
        fi
    done

done < region.tab
