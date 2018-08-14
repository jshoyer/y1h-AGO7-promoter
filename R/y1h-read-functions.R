
## JSH, Carrington Lab, DDPSC
## j.s.hoyer@wustl.edu
## Started Aug. 15 2013


#' ## Data from Gen5

#' For my experiments, I set Gen5 to automatically export files, which
#' by default go in order down each column of the plate, i.e.
#' column-by-column order/well-down order. (Gen5 refers to tall
#' columnar tables as "row-wise tables").
#' Export includes the number of the well (annoyingly prefixed with "SPL",
#' which could be confused with SQUAMOSA PROMOTER BINDING-LIKE factor
#' Arabidopsis gene codes!),
#' the alphanumeric well ID (A1, B1, ...),
#' an absorbance read at 600 nm,
#' and a luminescence read.

read.gen5 <- function(      # Reads plate data with current export format #2
                      file) ##<< The path to a single file to be read in
{
    plate <- read.table(file,
                        skip = 42, nrows = 384,
                        fill = TRUE,
                        col.names = c("well", "wellID", "A600", "Lum"))

    plate$well <- 1:384
    return(plate)
}


#' I have quite a bit of data is a similar text format,
#' from prior to 2013-08-09,
#' with one additional not particularly useful extra column,
#' and less info at the top of the file.

read.gen5.old <- function(      # Reads plate data with old export format #2
                          file) ##<< The path to a single file to be read in
{
    read.table(file,
               skip = 7, nrows = 384,
               col.names = c("Sample", "Well", "A600", "Lum", "above50RLU"))
}

#' Similar import function for A600 data:

readA600 <- function(filename) {
    read.delim(filename,
               skip = 28, header = FALSE,
               col.names = c("well", "wellID", "A600"))
}


#' Common text strings/Figure axis labels
rlu <- "Relative Luminescence Units"
activity <- "Luciferase activity (RLU per absorbance unit)"

numbycol <- "Well number (column-by-column order)"
numbyrow <- "Well number (row-by-row order)"
numbyzigzag <- "Well number (row-and-then-back zig-zag order)"
