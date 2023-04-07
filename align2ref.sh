#!/bin/bash

##make sure to run this file in your folder that has your index reference genome that you will be aligning to
proc=~/Consulting_Projects/Eaker_Deer/03_processrad ##input directory for process_rad sample files
al_out=~/Consulting_Projects/Eaker_Deer/04_Alignment ## output directory
bowtie_db=Odocoileus_hemionus ##change to your index name

##change file names below to your sample names.  You can pull from you samples/barcodes file
files="MYZY6
MYZY7
MYZY8
MYZY9
MYZYA
MYZYB
MYZYC
MYZYD
MYZYE
MYZYF
MYZYG
MYZYH
MYZYI
MYZYJ
MYZYK
MYZYK_R"

module load bowtie2 ## I use bowtie2, you can use bwa or another aligner
module load samtools ## samtools converts to bam and sorts our alignments

echo 'Aligning samples to' $bowtie_db > $al_out/bowtie.log

for file in $files; 
do
        echo 'Aligning' $file

        echo ${file} >> $al_out/bowtie.log

        ## I use --very-sensitive because divergent spp
        #change the -x Odocoileus_hemionus to your index name
        bowtie2 --very-sensitive -x Odocoileus_hemionus -1 $proc/${file}.1.fq.gz -2 $proc/${file}.2.fq.gz -S $al_out/${file}.sam 2>>$al_out/bowtie.log
        samtools view -b $al_out/${file}.sam | samtools sort --threads 8 -o $al_out/${file}.bam

done
