#!/usr/bin/bash

for i in `seq 1 500`; 
do 
    for j in `seq 0 1`; 
    do 
	cp Top_10_Each/kmeans_10_modeling${i}/all_models.9/${j}.rmf3 Top_2_Combined/${i}_${j}.rmf3;  
    done
done

