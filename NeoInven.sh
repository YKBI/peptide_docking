#!/bin/bash


pdb=$1
seq=$2
Ht=$3
outdir="$pdb"_"$seq"_"$Ht"
time python peptide_docking/peptide_docking.py $pdb $seq 100 peptide_docking/clucab.txt > "$outdir".logA

time python peptide_docking/crebfa.py "$pdb"_"$seq"_"$Ht" > "$outdir".logB

rm -rf "$outdir"/PDB_*/*.out "$outdir"/pdbs "$outdir"/PDB_*/*env
time tar cf - $outdir |pbzip2 > "$outdir".tar.bz2
ssh user1@10.1.5.9 "mkdir -p /arc10-2-03/200610_clucab_test/output/"
ssh user1@10.1.5.9 "mkdir -p /arc10-2-05/200610_clucab_test/output_comp/logs/"
time scp "$outdir".tar.bz2 user1@10.1.5.9:/arc10-2-05/200610_clucab_test/output_comp/
time scp "$outdir".log* user1@10.1.5.9:/arc10-2-05/200610_clucab_test/output_comp/logs/
time scp "$outdir"/*matrix/*total.txt user1@10.1.5.9:/arc10-2-03/200610_clucab_test/output/

rm -rf "$outdir" "$outdir".tar.bz2 "$outdir".log*

