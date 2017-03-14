# Data analysis of cotton MNase-seq datasets
---

## Preprocessing of DNA sequencing datasets

### MNase-seq datasets
Short-read data will be deposited in the NCBI short read archive ([SRP??????](http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=SRP??????)), also as Biobroject [PRJNA??????](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA??????).

#### Access to fastq files
    cd cottonLeaf/rawfastq
    ln -s ~/jfw-lab-new/RawData/HGJ_leafMNase-seq/WTNHHW163125/data_release/raw_data/*gz .

#### Checking read quality with [FastQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)
    cd ..
    module load fastqc/0.11.3
    fastqc rawfastq/*     ####bookmark
    # move result files to QCreport folder
    mkdir QCreport
    mv fastq/*zip QCreport/
    mv fastq/*html QCreport/
    
#### Quality trimming and adaptor removal
We usually use [Sickle](https://github.com/najoshi/sickle) to trim off sequences below quality threshold, and another popular tool is  [Fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). One more alternative to remove adpaters or primers is [cutadapt](https://cutadapt.readthedocs.io/). Based FastQC results above, illumina universal adaptor removal is necessay, and [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) appears to be a really nice wrapper tool of FastQC and Cutadapt, quite easy to use! 

    module load python 
    # python needed for Cutadapt
    trim_galore_v0.4.2/trim_galore --paired -o trimmed/ fastq/SRR2542701_1.fastq fastq/SRR2542701_2.fastq
    # check results
    grep 'Total reads processed' trimmed/*report.txt >trimmed/summary.txt
    grep 'Reads with adapters' trimmed/*report.txt >>trimmed/summary.txt
    grep 'Total written' trimmed/*report.txt >>trimmed/summary.txt

### Maize B73 reference genome
Although (Rodgers-Melnick et al. PNAS 2016) used maize B73 AGPv3 genome assembly, the current version is AGPv4 as describer [here](http://www.maizegdb.org/assembly). After downloading "Zea_mays.AGPv4.dna.toplevel.fa.gz", make [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) reference index for chromosomes (ignoring organellar genomes and scaffolds).

    mkdir maizeRef
    cd maizeRef
    zcat Zea_mays.AGPv4.dna.toplevel.fa.gz | head -35105648 >Zea_mays.AGPv4.chr10.fa
    module load bowtie2
    bowtie2-build Zea_mays.AGPv4.chr10.fa maize

## Read mapping and calling of hypersensive sites
(Rodgers-Melnick et al. PNAS 2016): "After the computational trimming of adaptor sequences using CutAdapt (40), paired-end reads were mapped to the maize B73 AGPv3 reference genome, using Bowtie2 with options “no-mixed,” “no-discordant,” “no-unal,” and “dovetail” (41) for each replicate digest and for the genomic DNA. BED files were made from the resulting BAM files, using bedtools bamtobed, filtered for minimal alignment quality (≥10), and read coverage in 10-bp intervals was calculated using coverageBed (42). The DNS values were obtained by subtracting the mean normalized depth (in reads per million) of the heavy digest replicates from those of the light digest replicates. In this way, positive DNS values correspond to MNase hypersensitive footprints (as defined by ref. 8; and referred to here as MNase HS regions), whereas negative DNS values correspond to nuclease hyper-resistant footprints (MRF, as per ref. 8). A Bayes factor criterion was used to classify as significantly hypersensitive."

### Bowtie2 mapping
Default setting `-k 1` report 1 alignment for each read/pair) should work, while some might need to be modified as required by downstream tools.

    mkdir mapping
    bowtie2 -q -p 6 -t --no-mixed --no-discordant --no-unal --dovetail -x maizeRef/maize -1 trimmed/SRR2542701_1_val_1.fq -2 trimmed/SRR2542701_2_val_2.fq -S mapping/SRR2542701.sam 2>SRR2542701.log
    samtools view -bS SRR2531768.sam | samtools sort - -o SRR2531768.sort.bam ; samtools index SRR2531768.sort.bam

* `-x maize`: use ref maize genome
* `-1 trimmed/SRR2542701_1_val_1.fq`: paired end read 1
* `-2 trimmed/SRR2542701_2_val_2.fq`: paired end read 2
* `-q`: takes fastq files
* `-p 6`: use 6 thread
* `-t`: Print the amount of wall-clock time taken by each phase.
* `--no-mixed --no-discordant`: discard discordant mapping
* `--no-unal`: Suppress SAM records for reads that failed to align.
* `--dovetail`: allow pair to overlap and over extend

### Differential nucleosome occupancy analysis - [DANPOS2](https://sites.google.com/site/danposdoc/)

### Profiling Nucleopostioning - [nucleR](http://bioconductor.org/packages/release/bioc/html/nucleR.html) & [NUCwave](http://nucleosome.usal.es/nucwave/)
Briefly, we need to examine the wave-length patterns of mapping coverage on the reference genome, and then specifically locate peaks that fit the description of nucleosomes - proper length, non-overlapping, etc.

(Zhang et al. 2015): "Well-positioned and loosely positioned nucleosomes were identified using nucleR (Flores and Orozco, 2011). ... We used the filterFFT function of nucleR to remove noise and smooth the read count score of each position along chromosomes with the parameter pcKeepComp = 0.01. After noise removal, nucleosome peaks and centers/dyads were determined using the peakDetection function (threshold = 25%, score = true, width = 140). Overlapped peaks were merged into longer regions, which were defined as loosely positioned nucleosomes, and distinct individual peaks were defined as well-positioned nucleosomes. If the length of merged peaks is longer than 150 bp, this region is considered to contain more than two nucleosome dyads and thus, contains loosely positioned nucleosomes. If the length of merged peaks is shorter than 150 bp, this region is considered to contain a well-positioned nucleosome."

The phasogram and average distance between two adjacent nucleosomes were calculated using our previously reported methods (Zhang et al., 2013). The nucleosome occupancy change scores were calculated by DANPOS (Chen et al., 2013). Analyses of dinucleotide frequency followed previously published methods (Locke et al., 2010; Valouev et al., 2011).

## Visualization and other result presentation
