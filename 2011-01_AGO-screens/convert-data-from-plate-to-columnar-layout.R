
#' * Conversion script, 2012
#' Steen Hoyer, Carrington Lab, DDPSC
#' Revised 2014-04-11.

#' Purpose:
#' Reformat data for data for cellHTS2, HTS Helper, etc.
#'
#' Input:
#' textfiles in 600nm
#'              420nm
#' These files start with a six line header,
#' then have a column label line (1, 2, .. 24),
#' followed by the plate reader data,
#' with a row label (A, B, .. P) as the first column
#' and the wavelength as the final column.
#'
#' Output:
#' textfiles in working-data-HTS-Helper,
#'              working-data-cellHTS2,
#'              and
#'              working-data-tall-420nm.txt
#'              working-data-tall-600nm.txt
#'
#' See also file:process-AGO-lacZ-screens-with-cellHTS2.R
#'
#' Order matters in processing these files,
#' because we need to apply the labels (short codes) in correct order.

### In ASCIIbetical order:
shortcodes <- scan("short-codes.txt", what = "character")

#' For the moment, for consistency with old scripts,
#' consider the AGOs in the order AGO1, AGO10, AGO7;
#' fragments in order from closest to the TSS to farthest.
shortcodes1107 <- shortcodes[c(5:8, 1:4, 9:12)]

## Convert from ASCIIbetical order to (AGO)numerical order: AGO1, AGO7, AGO10
#shortcodes1710 <- shortcodes[c(5:12, 1:4)]


dir.create("working-data-HTS-Helper")
dir.create("working-data-cellHTS2")

#' * Overly long function
#' The following just lets me go through the A600 and A400 files separately.
loop.through.plates <- function (directory, identifier)
{
    filenames <- dir(directory)

    number.of.plates <- length(filenames)
    total.number.of.wells <- number.of.plates * 384

    all.wells <- numeric()

                                        # 24 columns, 16 rows per plate
    well.IDs <- paste0(rep(LETTERS[1:16], each = 24), sprintf("%02d", 1:24))

    for (p in 1:number.of.plates) {
        ##print(paste("Plate", p, filenames[p]))
        plate.layout.data <- read.delim(paste0(directory, filenames[p]),
                                        skip=6, row.names = 1)
        ## for row-first columnar output, transpose the plate
        columnar.data <- c(as.vector(t(plate.layout.data[, -25])))
                                        # Omit the wavelength column
        all.wells <- c(all.wells, columnar.data)

        ## write data for cellHTS2
        dataframe.for.cellHTS2 <- cbind(filenames[p], well.IDs,
                                        columnar.data)
        ## 1. all data together
        write.table(dataframe.for.cellHTS2,
                    file = paste0("working-data-cellHTS2/",
                    filenames[p], ".txt"),
                    quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = FALSE)
        ## 2. individual screens - worth doing?
        ##system(paste0("mkdir -vp ../Data/split-for-cellHTS2/",
        ##filenames[p]
    }

    Plate  <- rep(1:60, each = 384)  ## 60 plates
    Row    <- rep(1:16, each = 24)   ## 24 rows on each plate
    Column <- rep(1:24)            ## 16 rows
    Value  <- all.wells
    HTS.helper.database.table <- cbind(Plate, Row, Column, Value)
    write.csv(HTS.helper.database.table,
              file = paste0("working-data-HTS-Helper/HTS_helper_database_table_",
              identifier, ".csv"),
              row.names = FALSE)

    write.table(cbind(rep(shortcodes1107, each=1920), well.IDs, all.wells),
                file = paste0("working-data-tall-",
                identifier, ".txt"),
                quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    ## Sanity check
    print(c("The value for the first well for each plate:",
            all.wells[(0:59) * 384 + 1]))
}

## For stepping through the above, use directory names below and
#' identifier <- "600nm"

directory1 <- "../2011-01_AGO-screens/600nm/"
loop.through.plates(directory1, "600nm")
directory2 <- "../2011-01_AGO-screens/420nm/"
loop.through.plates(directory2, "420nm")

getwd()
sessionInfo()
date()
