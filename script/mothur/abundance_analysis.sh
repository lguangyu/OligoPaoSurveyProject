#!/bin/bash

abund_dir="./abundance"
mkdir -p $abund_dir

# make abundance table and plot
for i in {phylum,class,order,family,genus}; do
	./analysis_script/get_mothur_taxonomy_abundance.py \
		-s mothur.output.shared \
		-t mothur.output.tax \
		-r $i \
	| tee $abund_dir/$i.tsv \
	| ./analysis_script/table_plot_mothur_taxonomy_abundance.py \
		--plot-percent \
		--max-n-taxa 40 \
		--plot-title $i \
		-p $abund_dir/$i.png
	# make transposed table
	./analysis_script/get_mothur_taxonomy_abundance.py \
		-s mothur.output.shared \
		-t mothur.output.tax \
		-r $i --transpose-output \
		-o $abund_dir/$i.transpose.tsv
done

# specifically looking for paos
for i in {genus,}; do
	./analysis_script/table_plot_mothur_taxonomy_abundance.py \
		--taxon-list ./analysis_script/data/pao.genus.stokholm_et_al_2017.txt \
		--plot-percent \
		--plot-title "putative PAO $i abundances" \
		-p $abund_dir/$i.tsv.pao.png \
		$abund_dir/$i.tsv
done
