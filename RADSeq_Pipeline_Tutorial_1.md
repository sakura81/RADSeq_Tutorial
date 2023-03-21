# RADSeq Pipeline Tutorial #1

### Signing into server
1. Sign into Cisco Anywhere Client (VPN)
2. Sign into server using your terminal program of choice (i.e. iTerm)

```{bash}
ssh username@servername.hpc.uidaho.edu
```

*  	Then enter password

### Basic Data Analysis
Within your directory, create a 0-RAW folder. Skip this step if you already have one.

```{bash}
mkdir 0-RAW
```

Transfer data from sequencing servers using wget command or file transfer programs (i.e. cyberduck) into you 0-RAW folder

* Unzip your fastq files

```{bash}
gunzip <filename>.fastq.gz 	#individual file
gunzip *.fastq.gz 	#all fastq files in a folder
```
* Determine the number of reads in each file

```{bash}
cat <filename.fastq> | echo $((`wc -l`/4))
```
This will provide you information on the number of reads in your R1 and R2 files.  You will be able to compare your loss of reads from quality or barcodes later on.

* Make a new directory for FastQC files

```{bash}
mkdir 01_FastQC
```
* FastQC your raw data results.  This can be done on the RAW files or after process radtags on each individual file.  I recommend doing both.  This can take ~15min - 1 hour to complete depedning on the size of your file.

```{bash}
module load fastqc
fastqc <filename> #make sure to do for the R1 and R2 files
```
* You can download and view the html file to look over the quality of your fastq files. To download either use the scp command or use a file transfer program.
* Make a new directory for flipped files

```{bash} 
mkdir 02_Flipped
```
* Copy the `flip_trim_160301.pl` to this folder using scp or a file transfer program.
* Create a .txt file containing all the barcodes that you used on the plate and save to this folder. Below is an example. See `All_Barcodes_03202023.txt` file.

```
AAACATCG
AACAACCA
AACCGAGA
AACGCTTA
AACGTGAT
AAGACGGA
```
* Run the flip script

```{bash}
module load perl
perl flip_trim_160301.pl <path_to_file_contining_barcodes.txt> <path_to_R1.fastq file> <path_to_R2.fastq file> <path_and_name_of_file_for_R1_flipped.fastq> <path_and_name_of_file_for_R2_flipped.fastq> -> <name_of_flipped_output.txt>
````
```
The output text file will contain how many reads it found per barcode
AAACATCG	3532784
AACAACCA	584984
AACCGAGA	2493033
AACGCTTA	1964089
AACGTGAT	434452
AAGACGGA	215426
AAGGTACA	976262
AATGTTGC	2726770
ACAAGCTA	4832175
```
* Determine the number of reads in each flipped file

```{bash}
cat <filename.fastq> | echo $((`wc -l`/4))
```
* You can now compare what your orignal file had for reads (unflipped version) to the flipped version to determine the loss of reads due to incorrect barcoding issues
```
100 - ((Original_reads - Flipped_reads)/Original_Reads)
```
* Make a new directory for PCR duplicate removal using clone_filter in stacks

```{bash}
mkdir 03_clone_filter
```
Why clone filter? 

PCR can stochasticaly amplify one of the alleles more than the other, sometimes leading to individuals appearing as homozygotes rather than heterozygotes (Andrews et al 2016). For this reason, removing PCR duplicates likely reduces false locus assembly and the variance in number of loci across individuals (Díaz-Arce & Rodríguez-Ezpeleta 2019). Downstream, this can have differing effects on population strcture depending on the dataset (Díaz-Arce & Rodríguez-Ezpeleta 2019). 

* Run clone filter on flipped R1 and R2 fastq files.

```{bash}
module load stacks
clone_filter -P -p <path_to_flipped_R1andR2_fastq_files> -i fastq -o <path_to_output_directory> &> <path_to_output_directory>clone_filter.log
```

### Process_radtags
* Make a new directory for process radtags

```bash
mkdir 04_processrad
```
* Create a .txt file containing 2 columns separated by a tab.  Column 1 will have the barcodes used on the plate and Column 2 will have the sample names for that barcode. Do not create a header row.  See `Barcode_Sample.txt` file a an example.

```
AAACATCG	PYRA1
AACAACCA	PYRA2
AACCGAGA	WOLF1
AACGCTTA	WOLF4
AACGTGAT	DOG1
AAGACGGA	RAIL2
AAGGTACA	FOX1
AATGTTGC	FOX2
ACAAGCTA	WOLF5
```
* Run process_radtags

```{bash}
module load stacks
process_radtags -P -c -q -e sbfI -b <path_to_barcode-sample>file.txt> -i fastq -1 <path_to_clone_filtered_flipped_R1.fastq> -2 <path_to_clone_filtered_flipped_R2.fastq> -o <path_to_where_you_want_output_files_sent> -y gzfastq
```
	-P	pair-end reads
	-c 	clean data, remove any read with an uncalled base
	-q  discard reads with low quality scores
	-e	restriction-enzyme used
	-b	path to a file containing barcodes for this run
	-1	first input file in a set of paired-end sequences.
	-2	second input file in a set of paired-end sequences
	-i	input file type
	-y	output file type
	-o	path to output the processed files
* Check the output log for the number of reads per barcode that passed quality filtering.

```{bash}
less process_radtags.log
```

### Alignments
* Make a directory for alignments

```{bash}
mkdir 05_Alignments
```
* Within the 05_Alignments create another directory containing your reference genomes and one containing your project name.

```{bash}
mkdir Reference_Genomes ProjectName
```
* Determine which refrence genome you want to align your data to and download the .fna or .fa version.  Using scp or file transfer program, copy the .fa or .fna genome file into this folder
* Build index of genome using the following command:

```{bash}
module load bowtie2
bowtie2-build <genome_file_name.fna> <index-name>
```


* Use `cd` to get to your directory containing the `align2ref.sh` file.  If not already there, use scp or file transfer program to copy to `05_Alignment/ProjectName` folder.  Make sure to edit the file with your project's file path and directory names before running.

* Run alignment script

```{bash}
align2ref.sh
```
### Proceed to RADSeq Tutorial #2 for refernce aligned SNP genotyping
