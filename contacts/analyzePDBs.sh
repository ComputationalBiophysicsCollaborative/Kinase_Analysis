#!/usr/bin/env bash

outdir=alignedPDB_peng_analysis
for weight in none unique m1 0.1 0.4
do
    mkdir -p $outdir/weight_$weight
    
    ../scripts/distFuncs.py ../structure/PDB_ABC/inList_filt ../structure/PDB_ABC/outList_filt 70 30 ../structure/PDB_ABC/pdb_counts_peng/weights_${weight}_ID ../alignment/regionsLong alignedPDB $outdir/weight_$weight 175
done

outdir=alignedPDB_vij_analysis
for weight in none unique m1 0.1 0.4
do
    mkdir -p $outdir/weight_$weight
    
    ../scripts/distFuncs.py ../SeqClasses/vijayan/DFG_IN_Ids ../SeqClasses/vijayan/DFG_OUT_Ids 70 30 ../structure/PDB_ABC/pdb_counts_vijayan/weights_${weight}_ID ../alignment/regionsLong alignedPDB $outdir/weight_$weight 175
done
