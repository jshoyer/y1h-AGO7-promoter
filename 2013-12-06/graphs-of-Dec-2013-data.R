
#' J. Steen Hoyer, Carrington lab, DDPSC.
#' Analysis started 2014-01-07
#'
#' This script is mainly meant
#' to concretely reproduce the figure in the paper.
#' I have some other exploratory analysis code --
#' email me if you have questions.

#' Eight Y1H assay plates scanned on 2013-12-06.
library(lattice)
source("../R/y1h-read-functions.R")
source("../R/plate-indexing-and-layouts.R")

relpath <- "../2013-12-06/"
filenames <- paste0(relpath, "plate0017", 1:8, ".txt")

table <- read.gen5(filenames[1])
for (i in 2:8) {
    table <- rbind(table, read.gen5(filenames[i]))
}


#' Two plate layouts,
#' four growth blocks processed for each.
#' The D block (for each layout) was not fed,
#' the A block was fed for 3 hours,
#' and the B and C blocks were fed for 5 hours.
table$layout <- rep(c(1, 2), each = 384)
table$block <- rep(c("D", "A", "B", "C"), each = 384 * 2)
table$feed <- c(rep(c("No feed", "3 hours"), each = 384 * 2),
                rep("5 hours", 384 * 4))

prey <- c(
    readLines("../2013-12-06/prey1.csv"),
    readLines("../2013-12-06/prey2.csv")
)
table$prey <- prey


#' Normalize luminescence measurements by cell density:
bgA600 <- median(table$A600[table$A600 < 0.1])
#' The median for the latter group of wells is `r bgA600`,
#' Subtract this value from A600 values to correct for background
#' absorbance, and mask the wells with A600 < 0.2
#' so that the low absorbance does not result in
#' disproportionately high gLUC activity values (RLU per absorbance unit),
#' from dividing by a very small number.
table$A600BS <- table$A600 - bgA600
table$A600BS[table$A600 < 0.2] <- NA

bgLum <- 7
table$LumBS <- table$Lum - bgLum
table$LumBS[table$Lum < 7] <- NA
table$activity <- table$LumBS / table$A600BS

table5 <- table[table$feed == "5 hours", ]
bymean <- with(table5, reorder(prey, activity, mean, na.rm = TRUE))

dotplot(bymean ~ activity, data = table5, xlab = activity, alpha = 0.25)
