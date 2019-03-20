#!/usr/bin/env bash


../scripts/alignContacts.py ./alignedPDB ../structure/PDB_ABC/distances/ ../structure/PDB_ABC/pdbseqIndices ../structure/PDB_ABC/pdbseqs_full_ID ../alignment/regionsLong

(
    cd alignedPDB
    mkdir combined
    ls aln/* | parallel convert {} pdb/{/} +append combined/{/}
)

#./DFGinout_Thresh.py ../SeqClasses/inList ../SeqClasses/outList 70 30 8 ../structure/weights_m1_ID ./aligned ./8A_70-30
