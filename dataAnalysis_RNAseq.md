# RNA-seq analysis corresponding to cotton MNase-seq datasets
---

## Preprocessing and mapping of RNA-seq datasets

### Locare FASTQ files
For both MNase-seq and RNA-seq experiments, mature leaf tissue was harvested from flowering branches at 5pm, and immediately flash frozen in liquid nitrogen and stored at −80°C. Four *Gossypium* accessions were used including a natural allopolyploid, *G. hirsutum* (AD1) cultivar Acala Maxxa, and models of its A- and D-genome diploid progenitors, *G. arboreum* (A2) and *G. raimondii* (D5), as well as the corresponding interspecific diploid F1 hybrid (A2 × D5).
    
    cd /home/hugj2006/jfw-lab/Projects/MNase-seq/cottonLeaf/RNAseq/rawfastq
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-Maxxa* .
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-D5* .
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-A2* .
    ls

### Quality trimming

    bash runTrimGalore.sh >runTrimGalore.071117.txt 2>&1
    grep -E "Total reads|Reads written" trimmed/*txt

### Cotton reference transcriptomes
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I will need 4 different reference genomes for AD1, A2, D5 and A2xD5. An improved AD1 reference became available on [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er).

    mkdir refTranscriptomes
    cd refTranscriptomes
    module load bowtie2
    ### D5_JGI - Paterson et al. 2012 Nature
    ln -s ~/jfw-lab/Projects/Eflen/seed_for_eflen_paper/bowtie2hylite/D5.transcripts* .
    # see: https://github.com/huguanjing/homoeologGeneExpression-Coexpression/blob/master/bowtie2hylite.sh
    
    
    ### A2_BGI - Li et al. 2014 Nature Genetics
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/
    ### AD1_NBI - Zhang et al, 2016 Nature biotechnology
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/TM1.fasta
    grep -n '>scaffold' TM1.fasta |head -10   #32244283:>scaffold27_A01
    head -32244282 TM1.fasta >TM1_26.fasta 
    ### make my own ref for A2xD5
    cat A2genome_13.fasta Dgenome2_13.fasta >F1_26t.fasta
    grep '>' F1_26t.fasta 
    sed 's/>Chr/>D5_chr/g' F1_26t.fasta >F1_26.fasta
    grep '>' F1_26.fasta
    rm F1_26t.fasta
    
    ### AD1_458 - Saski et al, 2017 (in revision)
    ln -s ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/v1.1/assembly/Ghirsutum_458_v1.0.fa.gz
    zcat Ghirsutum_458_v1.0.fa.gz |grep -n '>scaffold'|head -10 
    # 27134894:>scaffold_27
    zcat Ghirsutum_458_v1.0.fa.gz |head -27134893 >TM1new_26.fasta 
    
    # build bowtie2 ref
    bowtie2-build TM1_26.fasta TM1
    bowtie2-build TM1new_26.fasta TM1new
    bowtie2-build F1_26.fasta F1
    bowtie2-build A2genome_13.fasta A2
    bowtie2-build Dgenome2_13.fasta D5

### mapping

    bash runBowtie2.D.sh >runBowtie2.D.071417.txt 2>&1

????why does bowtie2 mapping rate so low???

## Expression profile of gene sets

