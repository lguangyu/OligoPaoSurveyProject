#!/bin/bash
#SBATCH -J MOTHUR2OLIGO
#SBATCH -pshort -N1 -c8 --time 12:00:00

# USAGE: sh mothur2oligo.sh
# This is a shell script for transforming mothur output to appropriate format for 
# A. murat Eren's oligotyping pipeline 

## Set variables

# Adjust the file names to your own study - these are the files from the mothur SOP
mothur="$HOME/opt/mothur/1.45.3/bin/mothur"
in_prefix="mothur.output.seqs"
processors=$SLURM_CPUS_PER_TASK

# Set the taxon you want to select, separate taxonomic levels with ";" 
# Do not touch inner and outer quotes
taxon="'$(cat extract_taxon)'"

# Call mothur commands for getting taxon-specific seqs
$mothur "#set.current(processors=$processors); get.lineage(taxonomy=${in_prefix}.taxonomy, taxon=$taxon, count=${in_prefix}.count_table);
	list.seqs(count=current);
	get.seqs(accnos=current, fasta=${in_prefix}.fasta)"

# re-align using mafft
$HOME/opt/mafft/7.505/bin/mafft \
	--thread $processors \
	${in_prefix}.pick.fasta > ${in_prefix}.pick.mafft.fasta

# Call mothur commands for generating deuniqued sequences
$mothur "#set.current(processors=$processors); deunique.seqs(fasta=${in_prefix}.pick.mafft.fasta, count=${in_prefix}.pick.count_table);"

# Replace all "_" in fasta header with a ":"
sed 's/_/:/g' ${in_prefix}.pick.redundant.groups > intermediate1
# Make a file which maps sample names to sequence headers
paste ${in_prefix}.pick.redundant.groups intermediate1 | awk 'BEGIN{FS="\t"}{print $1"\t"$2"_"$3}' > intermediate2

# Perl script to rename the headers of the fasta to include the sample name at the beginning followed by a "_"
perl ./script/renamer.pl ${in_prefix}.pick.mafft.redundant.fasta intermediate2 
ln -sfT ${in_prefix}.pick.mafft.redundant.fasta_headers-replaced.fasta final.fasta

# clean up
rm -f intermediate1 intermediate2
