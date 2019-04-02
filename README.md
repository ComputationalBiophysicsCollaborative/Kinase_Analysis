Scripts to process the PDB protein-kinase dataset.
Requires hmmer to be installed, and IvoGPU with the utils folder on the PATH.

Rough outline of steps to do the analysis:
```bash
cd structure
mkdir my_dataset
cd my_dataset
mv ~/my_pdb_file_collection PDB
../process.sh

cd ../../contacts
mkdir my_dataset
cd my_dataset
../processPDBs.sh ./../structure/my_dataset
../analyzePDBs.sh
```
