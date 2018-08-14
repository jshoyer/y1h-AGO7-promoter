
#' J. Steen Hoyer, Carrington lab, DDPSC.
#' Started 2014-07-24 (experiment end date; only data collection day).
#'
#' This script is mainly meant
#' to concretely reproduce the figure in the paper.
#' I have some other exploratory analysis code --
#' email me if you have questions.

library(lattice)
source("../R/y1h-read-functions.R")
source("../R/plate-indexing-and-layouts.R")

prey1 <- readLines("../2013-12-06/prey1.csv")
prey2 <- readLines("../2013-12-06/prey2.csv")
bait <- readLines("../2014-07-24/bait-layout.csv")
baitClone <- readLines("../2014-07-24/bait-layout-by-clone.csv")

dir("../2014-07-24/", pattern = "plate")
p523 <- read.gen5("../2014-07-24/plate00523.txt")
p524 <- read.gen5("../2014-07-24/plate00524.txt")
p525 <- read.gen5("../2014-07-24/plate00525.txt")
p526 <- read.gen5("../2014-07-24/plate00526.txt")

p523$prey <- prey1
p524$prey <- prey1
p525$prey <- prey2
p526$prey <- prey2
p523$bait <- bait
p524$bait <- bait
p525$bait <- bait
p526$bait <- bait
p523$baitClone <- baitClone
p524$baitClone <- baitClone
p525$baitClone <- baitClone
p526$baitClone <- baitClone
p523$plate <- as.factor(523)
p524$plate <- as.factor(524)
p525$plate <- as.factor(525)
p526$plate <- as.factor(526)
p523$preyLayout <- as.factor(1)
p524$preyLayout <- as.factor(1)
p525$preyLayout <- as.factor(2)
p526$preyLayout <- as.factor(2)

plates7.25 <- rbind(p523, p524, p525, p526)

noBait <- plates7.25[plates7.25$bait == "No bait cells", ]
matingControl <- plates7.25[plates7.25$bait == "YM4271", ]

bgA600 <- 0.039
bgLum <- median(c(noBait$Lum, matingControl$Lum))

plates7.25$A600BS <- plates7.25$A600 - bgA600
plates7.25$LumBS <- plates7.25$Lum - bgLum
plates7.25$A600BS[plates7.25$A600 < 0.2] <- NA
plates7.25$LumBS[plates7.25$A600 < 0.2] <- NA
plates7.25$activity <- plates7.25$LumBS / plates7.25$A600BS

p523 <- plates7.25[plates7.25$plate == 523, ]
p524 <- plates7.25[plates7.25$plate == 524, ]
p525 <- plates7.25[plates7.25$plate == 525, ]
p526 <- plates7.25[plates7.25$plate == 526, ]

######################################################################

nonautoactive <- plates7.25[plates7.25$baitClone != "SPLbs2 clone 2" &
                            plates7.25$baitClone != "pAGO7 -990/-446 clone 1" &
                            plates7.25$bait != "YM4271" &
                            plates7.25$bait != "No bait cells"
                          , ]

simpleTCPsubset <- nonautoactive[
                                 nonautoactive$bait == "pAGO7 -750/-476"|
                                 nonautoactive$bait == "pAGO7 -750/-501"|
                                 nonautoactive$bait == "pAGO7 -531/-446"|
                                 nonautoactive$bait == "TCPbs"          |
                                 nonautoactive$bait == "pAGO7 -990/-446",
                                 ]

dotplot(bait ~ activity, data = simpleTCPsubset[simpleTCPsubset$prey == "TCP2", ],
        main = "TCP2", xlab = activity, xlim = c(-100, 2500),)

dotplot(bait ~ activity, data = simpleTCPsubset[simpleTCPsubset$prey == "Empty", ],
        main = "Empty vector", xlab = activity, xlim = c(-100, 2500),)
