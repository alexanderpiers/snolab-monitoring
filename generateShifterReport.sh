#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N snolab_monitoring
#$ -l h=compute-3-*.local
#$ -pe mpi 24 

source /home/apiers/env/damicm/bin/activate


# Run monitoring code
python /home/apiers/snolab-monitor/shifterReport.py
