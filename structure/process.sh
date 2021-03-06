#!/usr/bin/env bash

datadir=$PWD
hmmfile=alignment.hmm
unaligned='0-4,180-240'
HMMALIGN=hmmalign
SCRIPTDIR=$(dirname "$0")

echo "$hmmfile  '$unaligned'" >datasource

$SCRIPTDIR/getPDBseqs.py PDB 100 >seq_log.txt
$HMMALIGN $SCRIPTDIR/$hmmfile pdbseqs.fa >pdbseqs_aligned.stk
$SCRIPTDIR/getAlnSeqs.py pdbseqs_aligned.stk orig --showid --unalign $unaligned >rawseqs
awk 'gsub(/-/,"-",$2)<32 {print}' rawseqs >pdbseqs_full_ID
rm pdbseqs_aligned.stk rawseqs

sort -u pdbseqs_full_ID -o pdbseqs_full_ID
awk '{gsub(/[a-z.]/,"",$2); print }' pdbseqs_full_ID >pdbseqs_ID
cut -d ' ' -f 1 pdbseqs_ID >IDs
cut -d ' ' -f 2 pdbseqs_ID >pdbseqs

$SCRIPTDIR/alnResinds.py pdbseqResinds pdbseqs_full_ID >alnMap

#now taken care of in pdb_counts
echo "" >Neffs
for weight in none unique m1 0.1 0.4
do
    echo $weight " " >>Neffs
    phyloWeights.py $weight pdbseqs weights_$weight >>Neffs
done

$SCRIPTDIR/filterPDBs.py IDs pdbseqResinds pdbseqs_full_ID >pdbfilt_log,txt

$SCRIPTDIR/getDists.py PDB 100 | tee distlog.txt

#Note: 2LGC.pdb removed because I can't process "MODEL" in pdbs
#Note: Corrected PHE to GLY at residue 550 in 3V5Q chain B
