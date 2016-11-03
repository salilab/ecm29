#!/usr/bin/bash

rm Random_Models_H1.dat Random_Models_H2.dat
rm Cluster_*.txt

norm1=0
norm2=0

for i in `seq 0 249`; 
do
    ls -1v ./all_models.999/${i}_0.rmf3 | awk 'BEGIN{FS="/"}{print $3}' >> Random_Models_H1.dat
    if [ -f ./all_models.999/${i}_0.rmf3 ]
    then
	(( norm1++ ))
    fi
    ls -1v ./all_models.999/${i}_1.rmf3 | awk 'BEGIN{FS="/"}{print $3}' >> Random_Models_H1.dat
    if [ -f ./all_models.999/${i}_1.rmf3 ]
    then
        (( norm1++ ))
    fi           
done

for i in `seq 250 499`;
do
    ls -1v ./all_models.999/${i}_0.rmf3 | awk 'BEGIN{FS="/"}{print $3}' >> Random_Models_H2.dat
    if [ -f ./all_models.999/${i}_0.rmf3 ]
    then
	(( norm2++ ))
    fi
    ls -1v ./all_models.999/${i}_1.rmf3 | awk 'BEGIN{FS="/"}{print $3}' >> Random_Models_H2.dat
    if [ -f ./all_models.999/${i}_1.rmf3 ]
    then
	(( norm2++ ))
    fi
done

echo $norm1, $norm2

for i in `seq 1 7`; 
do 
    for j in `seq 0 $(($i-1))`; 
    do 
	echo \
	    `while read p; do find ./kmeans_1000_${i}/cluster.${j} -name *-${p}; done < Random_Models_H1.dat | awk 'BEGIN{FS="/"}{print $3}' | awk 'BEGIN{FS="."}{print $2}' | sort | uniq -c | awk '{print $1 / 492}'` \
	    `while read p; do find ./kmeans_1000_${i}/cluster.${j} -name *-${p}; done < Random_Models_H2.dat | awk 'BEGIN{FS="/"}{print $3}' | awk 'BEGIN{FS="."}{print $2}' | sort | uniq -c | awk '{print $1 / 492}'` \
	    >> Cluster_${i}.txt
    done; 
done
