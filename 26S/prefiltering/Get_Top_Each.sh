#!/bin/bash

while read p; 
do
    python prefilter.py -nmods 5 -preload False -mpi False -nclusters 1 -prefilter 1000.0 -dir $p
done < dir.txt

mkdir Top_1000
for i in `seq 0 499`;
do
    rm kmeans_5_${i}/all_models.4/*.pdb
    cp kmeans_5_${i}/all_models.4/0.0.rmf3 Top_1000/${i}_0.rmf3
    cp kmeans_5_${i}/all_models.4/1.0.rmf3 Top_1000/${i}_1.rmf3
done



