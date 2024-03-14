#!/bin/bash
#SBATCH -J OLIGOTYPING
#SBATCH -pshort -N1 -c4

# include blastn and makeblastdb into path
PATH="$PATH:$HOME/opt/ncbi/blast+-2.11.0/bin"

positions=$(cat filtered_positions)
aln="mothur2oligo.fasta"

entropy=$aln"-ENTROPY"
out_dir=$aln".position_oligotype."$(echo $positions | sed 's/,/_/g')

. $HOME/.local/env/python-3.10.10-venv-generic/bin/activate

rm -rf $out_dir # clean up old results

oligotype -M 0 -s 3 -C $positions -N $SLURM_CPUS_PER_TASK -o $out_dir \
	$aln $aln"-ENTROPY"

ln -sfT $out_dir mothur2oligo.fasta.oligo_final

script/plot.oligo_size_histogram.py \
	-p $out_dir.oligo_size_histogram.png \
	$out_dir

deactivate
