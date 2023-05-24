#!/bin/bash
#SBATCH -J Bowtie2
#SBATCH --output=bowtie2-%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=staceyn@uidaho.edu

cd $SLURM_SUBMIT_DIR

projhome=~/Consulting_Projects/Swift_Fox
proc=$projhome/04_processrad ## these are your clone_filtered, demultiplexed reads
al_out=$projhome/05_Alignment ## output directory

module load bowtie2 ## I use bowtie2, you can use bwa or another aligner
module load samtools ## samtools converts to bam and sorts our alignments

# I made indices already from the Vulpes_vulpes genome
#bowtie2-build <reference>.fa <index-name>
bowtie_db=Vulpes_vulpes

#get file names
#ls *.rem.1.fq.gz|cut -d / -f 7 | cut -f 1 -d . > jobs ## make a list of samples to loop through
names=($(cat jobs))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

echo 'Aligning samples to' $bowtie_db > $al_out/bowtie.log

for i in $samples; do
        echo 'Aligning' $names

        echo ${names[${SLURM_ARRAY_TASK_ID}]} >> $al_out/bowtie.log

        ## I use --very-sensitive because divergent spp
        bowtie2 --very-sensitive -x Vulpes_vulpes -1 $proc/${names[${SLURM_ARRAY_TASK_ID}]}.1.fq.gz -2 $proc/${names[${SLURM_ARRAY_TASK_ID}]}.2.fq.gz -S $al_out/sam/${file}.sam 2>>$al_out/bowtie.log
        samtools view -b $al_out/sam/${names[${SLURM_ARRAY_TASK_ID}]}.sam | samtools sort --threads 8 -o $al_out/bam/${names[${SLURM_ARRAY_TASK_ID}]}.bam

done
