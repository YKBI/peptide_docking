#!/bin/bash


pdb=$1
seq=$2
Ht=$3
outdir="$pdb"_"$seq"_"$Ht"
time python peptide_docking/peptide_docking.py $pdb $seq 500 peptide_docking/DockingDB-0508.txt
rm -rf "$outdir"/PDB_*/*.out
time python peptide_docking/crebfa.py "$pdb"_"$seq"_"$Ht"

time tar cf - $outdir |pbzip2 > "$outdir".tar.bz2
ssh user1@10.1.5.9 "mkdir -p /arc10-2-03/200529_NeoInven/output/"
ssh user1@10.1.5.9 "mkdir -p /arc10-2-03/200529_NeoInven/output_comp/"
time scp "$outdir".ba2 user1@10.1.5.9:/arc10-2-03/200529_NeoInven/output_comp/
time scp "$outdir"/*matrix/*total.txt user1@10.1.5.9:/arc10-2-03/200529_NeoInven/output/


