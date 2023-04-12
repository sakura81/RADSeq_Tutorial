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
100 - ((Flipped Reads/Original_Reads)*100)
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
clone_filter -P -1 <path_to_flipped_R1.fastq> -2 <path_to_flipped_R2.fastq> -i fastq -y fastq -o <path_to_output_directory> &>> <path_to_output_directory>clone_filter.log
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
* Within the 05_Alignments create another directory containing your reference genomes for your project.

```{bash}
mkdir Reference_Genomes
```
* Determine which reference genome you want to align your data to and download the .fna or .fa version.  Using scp or file transfer program, copy the .fa or .fna genome file into this folder
* Build index of genome using the following command:

```{bash}
module load bowtie2
bowtie2-build <genome_file_name.fna> <index-name>
```

* Youl will now need to decide if you will run your alignments on the standalone servers or the cluster.  If you don't have very many files to run (<30), I would run on the standalone server.  If you have a large number of files (i.e. 2 plates worth), I would run on the cluster server.  They will use two different versions of the align2ref files so make sure you use the right one.

** STANDALONE SERVER OPTION

* Place the `align2ref.sh` file into `05_Alignment/Reference_Genomes` folder.  Make sure to edit the following locations before running
* 	proc= <path to your process_rad fastq files for each sample>
	al_out= <path to your output directory>
	bowtie_db= <index name of refernce genome you built in steps above>
	-x in bowtie2 command - change to the index name from above for your reference genome>
	
* Make sure the `align2ref.sh` file is executable (usually a different color in your terminal vs the files).  To do this:

```{bash}
chmod +x align2ref.sh
```
* Run alignment script

```{bash}
align2ref.sh
```
	
** CLUSTER SERVER OPTION

* First, generate a file jobs for the tasks by using the following command in your `~/04_processrad` folder
```{bash}
ls *.rem.1.fq.gz|cut -d / -f 7 | cut -f 1 -d . |cut -d_ -f1 > ~/05_Alignments/Reference_Genomes/jobs
```
* Now go to your	`~/05_Alignments/Reference_Genomes` folder and confirm that the `job` file is there.
* Upload the `align2ref.slurm` files into your `Reference_Genomes` folder.  Make sure to edit the following locations in the scrip:
```	
	#SBATCH --mail-user = <your email address>
	 projhome= ~/Consulting_Projects/Swift_Fox ##This is your project directory if you made one or it ca be '~'
	 proc=$projhome/04_processrad ## these are your clone_filtered, demultiplexed reads
	 al_out=$projhome/05_Alignment ## output directory
	 bowtie_db= <index name of refernce genome you built in steps above>
	 -x in bowtie2 command - change to the index name from above for your reference genome>
```
* Now sign into the cluster server
```{bash}
ssh username@fortyfive.hpc.uidaho.edu
```
* Enter your password
* Move to your `05_Alignments/Reference_Genomes` folder
* To begin running the file you will need to know how many samples are in your jobn file.  You can determine this by:
```{bash}
wc -l jobs
```	
To submit your job.  n = the # of samples in your job file:
```{bash}
sbacth -a 1-n align2ref.slurm
```
* The system will provide your with a job #.  You can check to make sure that your job is running by using
```{bash}
squeue --me
```
* You should see `n` jobs running and the length of time it has been running for.  Your files will be located in your `05_Alignments` folder.  There will be a `.sam` and a `.bam` file for each sample.
*You should receive an email once your job starts as well as when it finishes.  This will still take a lot of time (~10-12 hours)

### Proceed to RADSeq Tutorial #2 for refernce aligned SNP genotyping
