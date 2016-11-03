#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=24G
#$ -l arch=linux-x64
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
#$ -q lab.q
#$ -l hostname="!opt*"			    #-- anything but opt*
#$ -pe ompi 8
#$ -t 5
#$ -N cluster5
#########################################
  

source ~/.bash_profile
module load imp/last_ok_build

if [ -z $1 ]; then
    PREFILTER="1000.0"
else
    PREFILTER="$1"
fi
echo "PREFILTER = $PREFILTER"

if [ -z $2 ]; then
    NMODS="1000"
else
    NMODS="$2"
fi
echo "NMODS = $NMODS"


# write hostname and starting time 
hostname
date

# run
mpirun -np $NSLOTS python ./clustering.py -mpi True -preload False -nmods $NMODS -nclusters $SGE_TASK_ID -prefilter $PREFILTER

# done
hostname
date

