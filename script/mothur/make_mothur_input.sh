#!/bin/bash

sample_list="samples.list"
input_list="mothur.input.list"
fastq_dir="mothur.input.fastq"

# clean old mothur.input.list
# this avoid deleting it
echo -n '' > $input_list
for i in $(cat $sample_list); do
	echo -e "$i\t$fastq_dir/${i}_R1.fastq\t$fastq_dir/${i}_R2.fastq" >> $input_list
done
