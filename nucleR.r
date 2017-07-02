## ------------ NucleR analysis ---------------
# http://bioconductor.org/packages/release/bioc/vignettes/nucleR/inst/doc/nucleR.pdf
# 1. data preprocessing: SE needs fragment length to correct strand-specific mapping that shifts reads, PE do not shift reads but filter out longer fragments; trim narrows fragments to nucleosome core, improving peak calling.
# 2. filtering noise for coverage profile using Fast Fourier Transform (FFT)
# 3. peak detection for nucleosome positioning: 25% threshold checked first quantile data to detect peaks
# 4. define well arranged and loosely positioned nucleosomes
# 5. export peaks and coverage
## --------------------------------------------

# https://github.com/gthar/nucleServ/blob/master/bin/nucleR.R

library(data.table)
library(GenomicRanges)
library(nucleR)
library(travis)
library(rtracklayer)

thread =4
sl<-read.table("chr.size.txt")
source("nucleR.source.r")

## load bed files, loop over files
beds <- grep("*_q20.bed",list.files(),value=TRUE)
for(bedfile in beds)
{
    print(bedfile)
    
    message("...read bed files into data frame")
    df<- fread(bedfile, data.table=FALSE, header=FALSE, select=c(1,2,3,6))
    names(df)<-c("space","start","end","strand")
    
    message("...convert data frame to genome ranges")
    rd<- as(df, "RangedData")
    # nuckeR requires "RangedData" from "IRanges" package, or "AlignedRead" object containing the position, length and strand of the sequence reads.
    
    message("...processing reads")
    # -- Dyad Length = 50, trim each read to 50bp around the dyad
    prep <- processReads(rd, type="paired", fragmentLen=260, trim=50)
    
    message("...calculate the coverage, directly in reads per million (r.p.m)")
    cover <- coverage.rpm(prep)
    
    message("...plot coverage summary by chromosome")
    # check raw depth, are there enough coverage from sequencing?
    depth<-coverage(rd)
    pdf(paste0(gsub(".bed","",bedfile),".checkRawDepth.pdf"))
    par(mfrow=c(1,2))
    for(i in 1:length(depth)){
        chr<-names(depth)[i]
        x<-as.data.frame(table(depth[[i]]))
        x$cumsum<- cumsum(x$Freq)
        x$cumfrac<- x$cumsum/sl$V2[sl$V1==chr]
        plot(x$Var1[-1], 1-x$cumfrac[-1], type='n', xlab="coverage", ylab="fraction of bases sampled >= coverage", ylim=c(0,1.0), xlim=c(1,10), main=paste0(chr,": Target Region Coverage"))
        plot(x$Var1[-1], x$Freq[-1]/sl$V2[sl$V1==chr], type='n', xlab="coverage", ylab="frequence of bases sampled >= coverage", xlim=c(1,10), main=paste0(chr,": Target Region Coverage"))
    }
    dev.off()
    
    message("...write coverage into WIG")
    export.bw(cover,paste0(gsub(".bed","",bedfile),".coverage.bw"))
    # export.wig(cover, paste0(gsub(".bed","",bedfile),".coverage"))
    
    message("...filtering noise for coverage profile")
    emptyHandler <- function (f)
    function (x, ...)
    if (length(x) > 0) {
        f(x, ...)
    } else {
        numeric()
    }
    # 2% component is suggested for NGS
    fft <- mclapply(cover, emptyHandler(filterFFT), pcKeepComp = 0.02, mc.preschedule = FALSE, mc.cores = thread)
    
    message("...peak detection for nucleosome positioning")
    peaks <- peakDetection(fft, width  =147, threshold = "35%", score = FALSE, mc.cores  = thread)
    # Wwidth =147, Size of each nucleosome, in bp, to be considered
    # threshold = "35%",: Minimum number of reads (Coverage) to call a nucleosome. Can be given as a percentage (i.e., "25%" means that the peaks with coverage in the 1st quartile of data won't be considered); or as an absolute coverage value (i.e., "20" means that the peaks with less than 20 reads per million of mapped reads won't be considered). Default = 35%.
    
    message("...scoring peaks")
    scores <- peakScoring(peaks, fft, threshold = "35%", dyad.length = 50, mc.cores = thread)
    # peak width score (score_w) indicates positioning, higher means sharper peak
    # peak height score (score_h) indicates coverage, higher means talled peak
    
    message("...merging peaks")
    merged <- mergeCalls(scores, min.overlap = 50, mc.cores  = thread)
    # min.overlap - Minimum overlap between two reads for merge them
    # discard.low=0.2 - Discard low covered calls (i.e. calls with score_h < discard.low will be discarded)
    
    message("...peak classifications")
    # Well-positioned (W): score_w>0.4 & score_h>0.6; Otherwise fuzzy (F)
    getType <- function(score_w, score_h, thresh_w, thresh_h, nmerge=1)
    {
        ifelse(nmerge > 2, "uncertain", ifelse(`&`(score_w > thresh_w, score_h > thresh_h),"W", "F"))
    }
    merged$class <- getType(merged$score_w, merged$score_h, 0.4, 0.6, merged$nmerge)
    
    message("...write nucleosome positions into GFF")
    merged <- as.data.frame(merged)
    names(merged)[names(merged) == "space"] <- "seqname"
    names(merged)[names(merged) == "score_h"] <- "score_height"
    names(merged)[names(merged) == "score_w"] <- "score_width"
    merged$nmerge <- NULL
    writeGff(df2gff(merged, source  = "nucleR", feature = "Nucleosome"), paste0(gsub(".bed","",bedfile),".nucleosome.gff") )
    
    message("...summarize nucleosome classes")
    print(table(merged$class ))
}

q()
n











