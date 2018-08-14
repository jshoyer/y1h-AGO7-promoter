
#' mailto:j.s.hoyer@wustl.edu
#' Started 2017-06-25.

library(Biostrings)

## Download upstream sequences
## 34 MB file.
## Not strictly necessary for this analysis, but convenient.
#ftppath <-
#    "ftp://ftp.arabidopsis.org/Sequences/blast_datasets/TAIR10_blastsets/upstream_sequences/TAIR10_upstream_"
#download.file(paste0(ftppath, "1000_20101104"), "../2014_Weirauch-etc/1000")
#' 3000 and 500 nt versions also available
#prox1000 <- readDNAStringSet("../2014_Weirauch-etc/1000")
#ago7num <- grep("AT1G69440", prox1000@ranges@NAMES)
#prox1kb <- prox1000[[ago7num]]
#prox1kbrc <- reverseComplement(prox1kb)

proxAGO7 <- readDNAStringSet("../2014_Weirauch-etc/AGO7-upstream-1000.fasta")
prox1kb <- proxAGO7[[1]]
prox1kbrc <- reverseComplement(prox1kb)

#' ** Helper functions
#' Models are provided as position-specific frequency matrices.
#' Frequencies at each position sum to 1.
read.pwm <-
    function(file)
    {
        ## they put letters as column names.  matchPWM wants the opposite.
        matrixTall <- read.table(file, header = TRUE)
        log(t(matrixTall[ , -1]))
    }
#' We log each per-base "probability" to make things additive.
#' This converts our PFM to a PWM.
#' We will reverse this at the end,
#' to get back to a score
#' that can be interpreted as a likelihood
#' (probability that each sequence position is a match,
#' conditional on the model).

read.pwmByID <-
    function(motif_id)
{
        read.pwm(paste0("../2014_Weirauch-etc/pwm/",
                        motif_id, ".txt"))
}

##

plotScores <- function(...)
{
    plot(#ylab = "Score",
         #xlab = "Position upstream of transcription start site (bp)",
         col = rgb(0.5, 0.5, 0.5, 0.75),
         pch = 19, las = 1,
         ...)
}

dashedLinesAtTwoPositions <- function()
{
    #abline(v = -425.5, lty = "dotted", col = "blue")
    abline(v = -456.5, lty = "dotted", col = "blue")
    #abline(v = -751.5, lty = "dotted", col = "blue")
    abline(v = -766.5, lty = "dotted", col = "blue")
}

scanBothStrandsPlot <-
    function(motifID, ##<< CisBP ID
             sideBySide  = FALSE, ##<< FWD/REV scans next to each other
             oneOverOther = FALSE,##<< FWD/REV scans one on top of the other
             dashedLines  = TRUE, ##<< Indicate top-scoring positions for TCPs?
             ...      ##<< Additional arguments to plot(),
                      ##   via plotScores() helper function
             )
{
    pwm <- read.pwmByID(motifID)
    modelWidth <- dim(pwm)[2]
    positions <- 1:(1000 - modelWidth + 1)
    centers   <- positions + (modelWidth - 1)/2
    scoresFWD <- PWMscoreStartingAt(pwm, prox1kb  , positions)
    scoresREV <- PWMscoreStartingAt(pwm, prox1kbrc, positions)
    probScoresFWD <- exp(scoresFWD)
    probScoresREV <- exp(scoresREV)

    ## We want tick marks to reflect X-dimension,
    ## which goes from -1000 to -1
    ## upstream of the annotated TSS.
    centersInverted <- rev(-centers)

    marginVector <- c(5, 4, 2, 3) + 0.1
    if (sideBySide) {
        par(mfrow = c(1, 2), mar = marginVector)
    } else if(oneOverOther) {
        par(mfrow = c(2, 1), mar = marginVector)
    }

    plotScores(y = probScoresFWD,
               x = centersInverted,
               ...)
    if (dashedLines) dashedLinesAtTwoPositions()
    plotScores(y = rev(probScoresREV),
               x = centersInverted,
               ...)
    if (dashedLines) dashedLinesAtTwoPositions()
}

######################################################################

#attr(scanBothStrandsPlot, 'ex') <- function() {  #  Examples:
    dev.new(width = 6.5, height = 7)
    par(mfrow = c(5, 2),
        mar = c(3, 4, 1, 3) + 0.1)
    ##          /\
    ## Too small---x-axis label cut off.
    scanBothStrandsPlot("M1648_1.02", # TCP2
                        xlab = "")
    scanBothStrandsPlot("M1645_1.02",
                        xlab = "",
                        ylim = c(0, 2.25 * 10^-3)) # TCP3

    motifID <- "M1646_1.02" # TCP4
    scanBothStrandsPlot(motifID,
                        xlab = "",
                        ylim = c(0, 9 * 10^-4))

    scanBothStrandsPlot("M1650_1.02", # TCP5   -- 8 bp wide
                        xlab = "",
                        ylim = c(0, 0.25))
    scanBothStrandsPlot("M1644_1.02",
                        xlab = "Position upstream of transcription start site (bp)",
                        ylab = "Score", # "Score (10^-4)"
                        ylim = c(0, 1.7 * 10^-4)) # TCP24

    scanBothStrandsPlot("M1562_1.02", # SPL11
                        ylab = "Score (10^-6)",  # <-- Adjust manually!
                        ylim = c(0, 1.65 * 10^-6),
                        xlab = "Position upstream of transcription start site (bp)",
                        dashedLines = FALSE)
#}
