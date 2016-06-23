#!/bin/sh
#BSUB -q priority
#BSUB -J walsh
#BSUB -oo main.out
#BSUB -n 1
#BSUB -R "rusage[mem=12024]"
#BSUB -W 336:00

date

bcbio_nextgen.py batch1.yaml -n 16 -t ipython -s lsf -q mcore '-rW=336:00' -r mincores=2 -rminconcores=2 --retries 3 --timeout 180 --tag walsh  
# -r 'R=\"select[model!=XeonE52680]\"'

date
