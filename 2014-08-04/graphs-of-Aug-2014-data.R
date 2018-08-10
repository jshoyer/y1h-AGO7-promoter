
#' ** 2014-08-04 Y1H reads
#' J. Steen Hoyer, Carrington lab, DDPSC.
#' Started 2014-08-04.
#'
#' Summary: results consistent with the hypothesis that SPL11 binds
#' both 'GTAC' sites in our AGO7 promoter fragment of interest,
#' though muddled by the fact that deleting a 'TGGTCC' site reduces
#' activity too, albeit to a lesser extent.
#'
#' This script is mainly meant
#' to concretely reproduce the figure in the paper.
#' I have some other exploratory analysis code --
#' email me if you have questions.

library(lattice)
library(dplyr)
source("../R/y1h-read-functions.R")
source("../R/plate-indexing-and-layouts.R")

dir("../2014-08-04", full.names = TRUE)

# first read: 4.5 hour feed
p531 <- read.gen5("../2014-08-04/plate00531.txt")
# second read: 8.5 hour feed---much stronger luminescence
p532 <- read.gen5("../2014-08-04/plate00532.txt")
prey <- readLines("../2014-08-04/prey.txt")
bait <- readLines("../2014-08-04/bait.txt")
                                        # see labels.R
prey2 <- readLines("../2014-08-04/prey2.txt")
bait2 <- readLines("../2014-08-04/bait2.txt")
p531$bait <- bait
p531$prey <- prey
p532$bait <- factor(bait2,
                    ## The following is the order for graphs.
                    ## Use rev(), because lattice orders things upside-down
                    levels = rev(c(
                        "pAGO7 -750/-501",
                        "pAGO7 -750/-476",
                        "WT",
                        "TCPbs",
                        "SPLbs2",
                        "SPLbs1",
                        "SPLbs",
                        "8-bp")))
p532$prey <- prey2

#' Superdense well: define a subset without it.
p531sub <- p531[p531$A600 < 0.5, ]
p531[p531$A600 > 0.5, ]

#' background values
bgA600 <-  0.041
bgLum <- 4

p531$A600BS <- p531$A600 - bgA600
p531$LumBS <- p531$Lum - bgLum
p531$activity <- p531$LumBS / p531$A600BS
p532$A600BS <- p532$A600 - bgA600
p532$LumBS <- p532$Lum - bgLum
p532$activity <- p532$LumBS / p532$A600BS


#' For whatever reason,
#' luminescence values were unusually low for SPL10 wells
#' (in contrast to SPL11)
#' in this particular experiments.
p532sub <- filter(p532, prey != "SPL10")

dotplot(bait ~ activity, data = p532sub, groups = prey,
        auto.key = TRUE, xlab = activity)

#' I'm not quite as happy with the following kind of plot as I was for
#' the last experiment:
dotplot(bait ~ activity | prey, data = p532sub, auto.key = TRUE, xlab = activity)
