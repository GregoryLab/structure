
suppressPackageStartupMessages({
library("CSAR")
library("gplots")
})

#fix naming conventions 
ChIPseqScore = function (control, sample, backg = -1, file = NA, norm = 3 * 10^9, test = "Ratio", times = 1e+06, digits = 2) {
    if (length(intersect(as.character(control$filenames), as.character(sample$filenames))) > 
        0) {
        stop("Some files assigned to control and sample have the same name")
    }
    if (length(setdiff(sample$chr, control$chr)) != 0) {
        stop(paste(sample, "and", control, "have different chromosome names"))
    }
    normcontrol = 1
    normsample = 1
    if (norm != -1 & !is.na(norm)) {
        normcontrol = norm/sum(as.numeric(control$c1))
        normsample = norm/sum(as.numeric(sample$c1))
    }
    if (normcontrol < 1) {
        warning("The value of parameter \"norm\" is lower than w * number of mapped reads for the control.\nThis may decrease the range for the test score values, and cause problems calculating the FDR thresholds with permutations.\nPlease, consider to increase the value of the parameter \"norm\"")
    }
    nc <- sum(as.numeric(control$chrL))
    ns = sum(as.numeric(sample$chrL))
    ms <- sum(as.numeric(sample$c1))/ns * normsample
    mc <- sum(as.numeric(control$c1))/nc * normcontrol
    vc <- sum(as.numeric(control$c2))/(nc) * normcontrol^2 - 
        mc^2
    vs <- sum(as.numeric(sample$c2))/(ns) * normsample^2 - ms^2
    if (backg == -1) {
        backg <- ((sum(as.numeric(sample$c1))/sum(as.numeric(sample$chrL_0)) * 
            normsample - ms) * (vc/vs)^0.5) + mc
    }
    else {
        backg <- ((backg * normsample - ms) * (vc/vs)^0.5) + 
            mc
    }
    backg <- as.integer(round(max(1, backg)))
    filenames <- control$filenames
    for (i in 1:length(sample$chr)) {
        file1 <- file(description = paste(file, ".", sample$chr[i], ".CSARScore", sep = ""), "wb") #fixed this line
        filenames[i] <- paste(file, ".", sample$chr[i], ".CSARScore", 
            sep = "") #fixed this line
        con1 <- file(description = control$filenames[i], "rb")
        type = readBin(con = con1, what = "character")
        version = readBin(con = con1, what = "character")
        considerStrand = readBin(con = con1, what = "character")
        w = readBin(con = con1, what = "integer")
        uniquelyMapped = readBin(con = con1, what = "logical")
        uniquePosition = readBin(con = con1, what = "logical")
        chr = readBin(con = con1, what = "character")
        chrL = readBin(con = con1, what = "integer")
        con2 <- file(description = sample$filenames[i], "rb")
        type = readBin(con = con2, what = "character")
        version = readBin(con = con2, what = "character")
        considerStrand = readBin(con = con2, what = "character")
        w = readBin(con = con2, what = "integer")
        uniquelyMapped = readBin(con = con2, what = "logical")
        uniquePosition = readBin(con = con2, what = "logical")
        chr = readBin(con = con2, what = "character")
        chrL = readBin(con = con2, what = "integer")
        writeBin("CSARScore", con = file1)
        writeBin(version, con = file1)
        writeBin(considerStrand, con = file1)
        writeBin(w, con = file1)
        writeBin(uniquelyMapped, con = file1)
        writeBin(uniquePosition, con = file1)
        writeBin(as.character(chr), con = file1)
        writeBin(chrL, con = file1)
        j = 1L
        times <- as.integer(times)
        while (j <= chrL) {
            score1 = readBin(con = con1, n = times, what = "integer")
            score2 = readBin(con = con2, n = times, what = "integer")
            score1 <- as.numeric(score1 * normcontrol)
            score2 <- as.numeric(((score2 * normsample - ms) * 
                (vc/vs)^0.5) + mc)
            score1[score1 < backg] <- backg
            if (test == "Poisson") {
                score2 <- as.integer(round((-ppois(score2, score1, 
                  lower.tail = FALSE, log.p = TRUE)) * 10^digits))
            }
            if (test == "Ratio") {
                score2 <- as.integer(round((score2/score1) * 
                  10^digits))
            }
            writeBin(score2, con = file1)
            j <- j + times
        }
        close(con1)
        close(con2)
        close(file1)
        gc(verbose = FALSE)
        message(paste(sample$chr[i], " done..."))
    }
    info <- list(chr = sample$chr, chrL = sample$chrL, filenames = filenames, 
        digits = digits)
    rm(score1, score2)
    gc(verbose = FALSE)
    return(info)
}

args <- commandArgs(T)
if (length(args) < 7) {
	cat("USAGE: ./run_CSAR_saturation.R ds_tag ss_tag tag chr_len_fn thresholds_fn out_fn out_counts_fn\n")
	q()
}
ds_tag <- args[1]
ss_tag <- args[2]
tag <- args[3]
chr_len_fn <- args[4]
thresholds_fn <- args[5]
out_fn <- args[6]
out_counts_fn <- args[7]

# read in chromosome lengths
chr_len <- read.table(chr_len_fn, header=F, as.is=T, sep="\t")
colnames(chr_len) <- c("chr", "len")

# generate nhits data structures
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(ds_tag, "plus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsS.forward <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(ds_tag, "minus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsS.reverse <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(ss_tag, "plus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsC.forward <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(ss_tag, "minus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsC.reverse <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)

norm.forward <- sum(as.numeric(nhitsC.forward$c1))
norm.reverse <- sum(as.numeric(nhitsC.reverse$c1))
# run CSAR
score.forward <- ChIPseqScore(control=nhitsC.forward,sample=nhitsS.forward,file=sprintf("%s.plus",tag),norm=norm.forward,times=10000,test="Poisson")
score.reverse <- ChIPseqScore(control=nhitsC.reverse,sample=nhitsS.reverse,file=sprintf("%s.minus",tag),norm=norm.reverse,times=10000,test="Poisson")
win.forward<-sigWin(score.forward, t=1, g=1)
win.reverse<-sigWin(score.reverse, t=1, g=1)
win <- c(win.forward, win.reverse)

# read in FDR thresholds (from shuffled data)
data.thresholds <- read.table(thresholds_fn, header=F, as.is=T, sep="\t")
colnames(data.thresholds) <- c("FDR", "threshold", "shuffled_count")

# print out results to BED file
results <- data.frame(as.vector(seqnames(win)))
colnames(results) <- "chrom"
results$start <- start(ranges(win))-1
results$end <- end(ranges(win))
results$id <- unlist(lapply(1:length(win), function(x) sprintf("CSAR_peak_%d",x)))
results$score <- score(win)
results$strand <- c(rep("+",length(win.forward)), rep("-",length(win.reverse)))
write.table(results, file=out_fn, quote=F, sep="\t", row.names=F, col.names=F)

results.fdr05 <- subset(results, score>=data.thresholds$threshold[which(data.thresholds$FDR==0.05)])
write.table(results.fdr05, file=gsub(".bed$", "_FDR05.bed", out_fn), quote=F, sep="\t", row.names=F, col.names=F)

# count how many PPSs are real at each FDR threshold and print results
data.thresholds$PPS_count <- unlist(lapply(data.thresholds$threshold, function(x) length(which(results$score>=x))))
write.table(data.thresholds[,c("FDR","threshold","PPS_count")], file=out_counts_fn, quote=F, sep="\t", row.names=F, col.names=F)




