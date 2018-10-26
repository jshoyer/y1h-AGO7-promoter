
#' J. Steen Hoyer, Carrington lab, DDPSC.
#'
#' Input:
#'   Working data created with convert-data-from-plate-to-columnar-layout.R
#'   names_and_families.txt
#'   lists of empty wells (below)
#'
#' Output:
#'   A dataframe, y1h,
#'   that is *not currently saved*.
#'   Also save a CSV table of counts by family
#'   and sorted tables
#'   for all twelve promoter fragments.

library(dplyr)  # Only used at end of script.

A600 <- read.table("../2011-01_AGO-screens/working-data-tall-600nm.txt")
A420 <- read.table("../2011-01_AGO-screens/working-data-tall-420nm.txt")

AGOs <- c("AGO1","AGO10","AGO7")
AGO <- rep(AGOs, each = 7680)
frags <- c("-0585 to -1","-1170 to -536",
                  "-1755 to -1121","-2307 to -1706",
                  "-0520 to -1", "-1040 to -471",
                  "-1560 to -991","-2033 to -1511",
                  "-0495 to -1","-0990 to -446",
                  "-1485 to -941","-1931 to -1436")
fragment <- rep(frags, each = 1920)
row <- rep(1:16, each = 24)
column <- 1:24
namesFam <- read.table("../2011-01_AGO-screens/names_and_families.txt",
                       na.strings = "",
                        sep = "\t", blank.lines.skip = FALSE)

#' head(namesFam)
#' Long list:
table(namesFam$V2)
#' There were supposed to be 1678 full wells:
length(na.omit(namesFam$V1))  # AGI
length(na.omit(namesFam$V2))
96 * 20
96 * 17

y1h <- data.frame("AGI" = namesFam[ , 1], "Family" = namesFam[ , 2],
                  plateNum = as.factor(rep(1:5, 12, each = 384)),
                  "row" = as.factor(row), "column" = as.factor(column),
                  "AGO" = AGO, "fragment" = as.factor(fragment),
                  "Well" = A420[ , 2], "A600" = A600[ , 3],
                  "A420" = A420[ , 3],
                  "betagal" = A420[ , 3]/A600[ , 3],
                  stringsAsFactors = FALSE)
y1h$fullfragment <- as.factor(paste(y1h$AGO, y1h$fragment))
y1h$state <- "TF-AD"
y1h$state[
    (y1h$plateNum == 1 & y1h$Well %in% c("C01","C16","N13","B12"))|
    (y1h$plateNum == 2 & y1h$Well %in% c("K13","E08","F17","P16"))|
    (y1h$plateNum == 3 & y1h$Well %in% c("P23","B18"))|
    (y1h$plateNum == 4 & y1h$Well %in% c("C01","L23","P10"))|
    (y1h$plateNum == 5 & y1h$Well %in% c("G03","O16"))] <- "pEXP-AD502"
identical(which(y1h$state == "pEXP-AD502"),
          which(y1h$AGI ==   "pEXP-AD502"))

empty1 <- read.table("../2011-01_AGO-screens/listEmptyWells/Plate1.txt")
empty2 <- read.table("../2011-01_AGO-screens/listEmptyWells/Plate2.txt")
empty3 <- read.table("../2011-01_AGO-screens/listEmptyWells/Plate3.txt")
empty4 <- read.table("../2011-01_AGO-screens/listEmptyWells/Plate4.txt")
empty5 <- read.table("../2011-01_AGO-screens/listEmptyWells/Plate5.txt")

#' 174 strains did not grow:
length(c(empty1$V1, empty2$V1, empty3$V1, empty4$V1, empty5$V1))
#' Number of strains that did grow when spotted before experiment:
#' 1678 - 174
#' (Not a unique count---includes a spot for the control on each plate.
#' There were 15 "no DNA-binding domain" control wells,
#' as can be seen above.
length(which(y1h$AGI ==   "pEXP-AD502"))/12
#' All grew when spotted---cf. below.
#'
#' The following count is off by eight,
#' relative to the count below.
#' 1678 - 174 - 15

y1h$didntgrow <-
    (y1h$plateNum == 1 & y1h$Well %in% empty1[, 1]) |
    (y1h$plateNum == 2 & y1h$Well %in% empty2[, 1]) |
    (y1h$plateNum == 3 & y1h$Well %in% empty3[, 1]) |
    (y1h$plateNum == 4 & y1h$Well %in% empty4[, 1]) |
    (y1h$plateNum == 5 & (y1h$Well %in% empty5[, 1] |
     y1h$row %in% seq(2, 16, by = 2))) # even numbered rows empty

y1h$state[y1h$didntgrow] <- "Empty"
y1h$betagal[y1h$didntgrow] <- NA

firstFivePlates <- y1h[1:(5 * 384), ]
tfActuallyGrew  <- firstFivePlates[!firstFivePlates$didntgrow, ]
#' str(tfActuallyGrew)
#' tail(tfActuallyGrew)
#' write.csv(tfActuallyGrew, "working-data-strains-that-grew.csv")
tfActuallyGrewFiltered <- dplyr::filter(tfActuallyGrew, !is.na(AGI))
#' str(tfActuallyGrewFiltered)
#' 44 wells (in odd-numbered rows)
#' came from half-empty final 96-well plate:
#' 1556 - 1512
actualTable <-
    table(tfActuallyGrewFiltered$Family, useNA = "always")
#' str(actualTable)
write.csv(actualTable,
          "../2011-01_AGO-screens/working-data-counts-by-family-not-including-strains-that-did-not-grow.csv")

y1h$betagal[y1h$A600 < 0.2] <- NA
y1h$state[y1h$A600 < 0.2 & ! y1h$state %in% c("pEXP-AD502","Empty")] <- "Masked"


y1h$logbetagal <- log(y1h$betagal)

plateMedians <- vector(length = 60)
for (i in 0:59) {
    plateMedian <- median(y1h$betagal[i*384 + 1:384], na.rm = TRUE)
    plateMedians[i + 1] <- plateMedian
    y1h$diffMedian[i*384 + 1:384] <- y1h$betagal[i*384 + 1:384] - plateMedian
    y1h$foldChange[i*384 + 1:384] <- y1h$betagal[i*384 + 1:384]/plateMedian
}

MADs <- vector(length = 60)
for (i in 0:59) {
    MAD <- median(abs(y1h$diffMedian[i*384 + 1:384]), na.rm = TRUE)
    MADs[i+1] <- MAD
    y1h$MADsAboveMedian[i*384 + 1:384] <- y1h$diffMedian[i*384 + 1:384] / MAD
}
sorted <- y1h[order(y1h$MADsAboveMedian, decreasing = TRUE), ]

medianAcrossScreens <- vector(length = 1920)
for (i in 1:1920) {
    medianAcrossScreens[i] <- median(y1h$foldChange[i + 384 * 5 * 0:11],
                                 na.rm = TRUE)
}

#' Probably a noisy measure:
y1h$medianAcrossScreens <- medianAcrossScreens
y1h$crossScreenMedianNormalized <-
    y1h$foldChange / medianAcrossScreens

nonspecificAutoactivatorWells <-
    filter(y1h, medianAcrossScreens > 2)
y1hSansNonspecificAutoactivators <-
    filter(y1h, medianAcrossScreens < 2)

#' library(dplyr)
allSorted <- arrange(
    y1hSansNonspecificAutoactivators,
    desc(MADsAboveMedian))
# Cf. 'sorted' data frame, just above.
allSorted <- filter(allSorted,
                    state %in% c("TF-AD", "pEXP-AD502"))
allSorted <- filter(allSorted,
                    ! is.na(betagal))
# pEXP-AD502 wells for which A600 < 0.2
allSorted <- select(allSorted,
                    -state, -didntgrow)

fragments <- paste(
    rep(AGOs, each = 4), frags)

dir.create("working-data")

# Twelve sorted tables for supplemental Excel file,
# deposited in Zenodo:
# https://doi.org/10.5281/zenodo.1472235
for (i in 1:12) {
    tableForFragment <- filter(
        allSorted, fullfragment == fragments[i])
    ## Column stuck in the middle,
    ## such that it is easiest to just remove it here:
    tableForFragment <- select(tableForFragment, -fullfragment)
    write.csv(tableForFragment,
              paste0("working-data/", fragments[i], ".csv"))
}
