#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=24G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
#$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
#$ -pe ompi 16
##$ -t 5
#$ -t 1-10                           #-- specify the number of tasks
#$ -N rmsf_precision
#########################################
  
source ~/.bash_profile
module load python/scikit/0.15.2
module load python/matplotlib/1.2.0
module load openmpi-x86_64
module load relion/1.4

export OMP_NUM_THREADS=$NSLOTS


echo "NSLOTS = $NSLOTS"
echo "JOB_ID = $JOB_ID"
echo "SGE_TASK_ID = $SGE_TASK_ID"

if [ -z $1 ]; then
    NMODS="500"
else
    NMODS="$1"
fi
echo "NMODS = $NMODS"


# write hostname and starting time 
hostname
date

echo "$IMP python ./precision_rmsf.py -test False -dir kmeans_${NMODS}_${SGE_TASK_ID}"
mpirun -np $NSLOTS python ./precision_rmsf.py -test False -dir kmeans_${NMODS}_${SGE_TASK_ID}

# done
hostname
date

