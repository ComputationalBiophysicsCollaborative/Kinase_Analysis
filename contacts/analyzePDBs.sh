#!/usr/bin/env bash
SCRIPTDIR=$(dirname "$0")
struct_dir=$1

for dataset in CA CB NHA SCNHA SCC
do

for weight in 0.0 0.1 0.4
do
    mkdir -p contacts$dataset/weight_$weight
    
    $SCRIPTDIR/contactFreq.py $struct_dir/pdblist_filt $struct_dir/pdbseqs_ID $weight $SCRIPTDIR/regionsLong ./aligned$dataset contacts$dataset/weight_$weight 175
done

done
