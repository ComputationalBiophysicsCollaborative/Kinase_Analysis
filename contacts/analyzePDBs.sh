#!/usr/bin/env bash

#for dataset in CA CB NHA SCNHA SCC
for dataset in NHA
do

for weight in 0.0 0.1 0.4
do
    mkdir -p contacts$dataset/weight_$weight
    
    ./contactFreq.py ../structure/pdblist_filt ../structure/pdbseqs_ID $weight regionsLong ./aligned$dataset contacts$dataset/weight_$weight 175
done

done
