#!/bin/bash

# get the compositional data transform for dirname 
# i.e. if called in `analysis.log-clr`, the transform is `log-clr`
basename $(pwd) | grep -oE '[^.]*$' > compositional_transform.txt

# transform data into abundance/compositional
trans_meth=$(cat compositional_transform.txt)
pseudo_count=0.10
abund_thres=0.10
## acinetobacter
#script/transform.oligo_abund.py \
#	--pseudo-count $pseudo_count \
#	--transform $trans_meth \
#	--output data/acinetobacter.compositional.tsv \
#	--output-abund-table data/acinetobacter.abund.tsv \
#	--output-alpha-diversity data/acinetobacter.alpha_diversity.tsv \
#	--oligo-rename-prefix ACIN \
#	--oligo-rename-result data/acinetobacter.rename.txt \
#	--abund-oligo-thres $abund_thres \
#	--abund-oligo-output data/acinetobacter.abund_oligo.list \
#	../oligo.acinetobacter/oligotyping/mothur2oligo.fasta.oligo_final/MATRIX-COUNT.txt
## candidatus_accumulibacter
#script/transform.oligo_abund.py \
#	--pseudo-count $pseudo_count \
#	--transform $trans_meth \
#	--output data/candidatus_accumulibacter.compositional.tsv \
#	--output-abund-table data/candidatus_accumulibacter.abund.tsv \
#	--output-alpha-diversity data/candidatus_accumulibacter.alpha_diversity.tsv \
#	--oligo-rename-prefix ACC \
#	--oligo-rename-result data/candidatus_accumulibacter.rename.txt \
#	--abund-oligo-thres $abund_thres \
#	--abund-oligo-output data/candidatus_accumulibacter.abund_oligo.list \
#	../oligo.candidatus_accumulibacter/oligotyping/mothur2oligo.fasta.oligo_final/MATRIX-COUNT.txt
## dechloromonas
#script/transform.oligo_abund.py \
#	--pseudo-count $pseudo_count \
#	--transform $trans_meth \
#	--output data/dechloromonas.compositional.tsv \
#	--output-abund-table data/dechloromonas.abund.tsv \
#	--output-alpha-diversity data/dechloromonas.alpha_diversity.tsv \
#	--oligo-rename-prefix DECH \
#	--oligo-rename-result data/dechloromonas.rename.txt \
#	--abund-oligo-thres $abund_thres \
#	--abund-oligo-output data/dechloromonas.abund_oligo.list \
#	../oligo.dechloromonas/oligotyping/mothur2oligo.fasta.oligo_final/MATRIX-COUNT.txt
## tetrasphaera
#script/transform.oligo_abund.py \
#	--pseudo-count $pseudo_count \
#	--transform $trans_meth \
#	--output data/tetrasphaera.compositional.tsv \
#	--output-abund-table data/tetrasphaera.abund.tsv \
#	--output-alpha-diversity data/tetrasphaera.alpha_diversity.tsv \
#	--oligo-rename-prefix TETR \
#	--oligo-rename-result data/tetrasphaera.rename.txt \
#	--abund-oligo-thres $abund_thres \
#	--abund-oligo-output data/tetrasphaera.abund_oligo.list \
#	../oligo.tetrasphaera/oligotyping/mothur2oligo.fasta.oligo_final/MATRIX-COUNT.txt


# plot
for orgn in {acinetobacter,candidatus_accumulibacter,dechloromonas,tetrasphaera}; do
	## plot alpha diversity
	#script/plot.alpha_diversity.py \
	#	-p $orgn.alpha_diversity.png \
	#	data/$orgn.alpha_diversity.tsv

	## plot abundance hca
	#script/plot.oligo_abund_with_hca.py \
	#	-p $orgn.oligo_abund_with_hca.png \
	#	-l data/$orgn.abund_oligo.list \
	#	--colorbar-label $trans_meth \
	#	data/$orgn.compositional.tsv

	## plot abundance decomposition
	#script/plot.oligo_abund_decomp.py \
	#	-m t-sne \
	#	-p $orgn.oligo_abund_decomp.png \
	#	data/$orgn.compositional.tsv

	# calculate anosim
	script/calc.wwtp_anosim.py \
		-o $orgn.wwtp_anosim.json \
		data/$orgn.compositional.tsv

	#break
done


## cross-hca
#for orgn1 in {acinetobacter,candidatus_accumulibacter,dechloromonas,tetrasphaera}; do
#	for orgn2 in {acinetobacter,candidatus_accumulibacter,dechloromonas,tetrasphaera}; do
#		if [[ $orgn1 != $orgn2 ]]; then
#			script/plot.oligo_abund_cross_hca.py \
#				-m euclidean \
#				-1 data/$orgn1.compositional.tsv \
#				--oligo-label-1 $orgn1 \
#				--oligo-list-1 data/$orgn1.abund_oligo.list \
#				-2 data/$orgn2.compositional.tsv \
#				--oligo-label-2 $orgn2 \
#				--oligo-list-2 data/$orgn2.abund_oligo.list \
#				-p cross_hca.$orgn1.$orgn2.png
#		fi
#	done
#done

