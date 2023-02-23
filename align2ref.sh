#!/usr/bin/env bash

projhome=<path to home directory>

proc=$projhome/04_processrad ## these are your clone_filtered, demultiplexed reads
al_out=$projhome/05_Alignments/<Project_Name> ## output directory
mkdir -p $al_out ## make a directory to hold the output if it doesn't already exist

module load bowtie2 ## I use bowtie2, you can use bwa or another aligner
module load samtools ## samtools converts to bam and sorts our alignments

# You already made indices from your reference genome
#bowtie2-build <reference>.fa <index-name>
bowtie_db=$projhome/05_Alignments/Reference_Genomes/<index-name>

samples=$(ls $proc/*rem.1.fq.gz) ## make a list of samples to loop through

rm $al_out/bowtie.log

echo 'Aligning samps to' $bowtie_db > $al_out/bowtie.log

for i in $samples; do
        samp=$(echo $i | cut -d / -f 7 | cut -f 1 -d .)
        echo 'Aligning' $samp

        echo $samp >> $al_out/bowtie.log

        ## I use --very-sensitive because divergent spp
        bowtie2 --very-sensitive -x $bowtie_db \
         -1 $proc/${samp}.1.fq.gz -2 $proc/${samp}.2.fq.gz \
        -S $al_out/${samp}.sam 2>> $al_out/bowtie.log
         samtools view -b $al_out/${samp}.sam | samtools sort --threads 8 \
         -o $al_out/${samp}.bam

done
