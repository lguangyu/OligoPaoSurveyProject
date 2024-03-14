#!/bin/bash
#SBATCH -J OLIGO_MSA
#SBATCH -o /dev/null
#SBATCH -pshort -N1 -c12 --time 8:00:00

# include blastn and makeblastdb into path
PATH="$PATH:$HOME/opt/ncbi/blast+-2.11.0/bin"

positions=$(cat filtered_positions)
aln="mothur2oligo.fasta"

# remove old results
msa_dir=$aln".msa_test"
rm -rf $msa_dir
mkdir $msa_dir
res_file=$aln".msa_test.tsv"
rm -f $res_file
rm -f $res_file.tmp

# the parameter value to test
msa_list="msa.list"

. $HOME/.local/env/python-3.10.10-venv-generic/bin/activate
# iterative test -M
entropy=$aln"-ENTROPY"
for msa in $(cat $msa_list); do
    # output directory for each run
    out_dir=$msa_dir"/"$msa
    # run test
    oligotype -M $msa -s 3 -C $positions -N $SLURM_CPUS_PER_TASK -o $out_dir \
        $aln $entropy
    # extract results
	if [[ -s $out_dir"/HTML-OUTPUT/index.html" ]]; then
		res=$(grep -Ee '-M.*elimination' $out_dir"/HTML-OUTPUT/index.html" |\
    	    grep -oP '\d+\D*$' | grep -oP '\d+')
		echo -e "$msa\t$res" >> $res_file.tmp
	fi
done
deactivate
