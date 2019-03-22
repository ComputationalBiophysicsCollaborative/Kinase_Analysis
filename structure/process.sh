#!/usr/bin/env bash

datadir=$PWD
hmmfile=alignment.hmm
unaligned='0-4,180-240'
HMMER=~/hmmer/binaries/hmmalign

echo "$hmmfile  '$unaligned'" >datasource

./getPDBseqs.py PDB 100 >seq_log.txt
$HMMER $hmmfile pdbseqs.fa >pdbseqs_aligned.stk
./getAlnSeqs.py pdbseqs_aligned.stk orig --showid --unalign $unaligned >rawseqs
awk 'gsub(/-/,"-",$2)<32 {print}' rawseqs >pdbseqs_full_ID
rm pdbseqs_aligned.stk rawseqs

sort -u pdbseqs_full_ID -o pdbseqs_full_ID
awk '{gsub(/[a-z.]/,"",$2); print }' pdbseqs_full_ID >pdbseqs_ID
cut -d ' ' -f 1 pdbseqs_ID >IDs
cut -d ' ' -f 2 pdbseqs_ID >pdbseqs

#now taken care of in pdb_counts
echo "" >>Neffs
for weight in none unique m1 0.1 0.4
do
    echo $weight " " >>Neffs
    phyloWeights.py $weight pdbseqs weights_$weight >>Neffs
done

#../../scripts/alphaMap.py pdbseqs ../../alignment/map8 >pdbseqs8
#paste IDs pdbseqs8 >pdbseqs8_ID

#../getDists.py PDB

#Note: 2LGC.pdb removed because I can't process "MODEL" in pdbs
#Note: Corrected PHE to GLY at residue 550 in 3V5Q chain B
