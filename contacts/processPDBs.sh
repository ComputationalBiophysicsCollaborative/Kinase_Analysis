#!/usr/bin/env bash

SCRIPTDIR=$(dirname "$0")
struct_dir=$1

for dataset in CA CB NHA SCNHA SCC
do
    $SCRIPTDIR/alignContacts.py ./aligned$dataset $struct_dir/distances$dataset/ $struct_dir/pdbseqIndices $struct_dir/pdbseqs_full_ID $SCRIPTDIR/regionsLong &
done

wait

#(
#    cd alignedPDB
#    mkdir combined
#    ls aln/* | parallel convert {} pdb/{/} +append combined/{/}
#)
