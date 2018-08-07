
#' J. Steen Hoyer, Carrington lab, DDPSC.
#'
#' Input:
#'   data frame 'y1h'
#'
#' Output:
#'   Graphs, in PDF files


source("convert-data-from-plate-to-columnar-layout.R")
#' Or
#' knitr::spin("convert-data-from-plate-to-columnar-layout.R")
source("workup_Feb2012.R")

promoterFragment <- paste(rep(AGOs, each = 4), frags)

#---------------------------------------------------------------------
# mark autoactivators in gray
colors <- rep("black", 1920)   # 384 wells * 5 plates = 1920
colors[medianAcrossScreens > 2] <- "gray"

plot.2D <- function(fragNum, ##<< 1 to 12
                    ...      ##<< Additional args for plot()
                    )
{
    plot(x = medianAcrossScreens,
         y = y1h$foldChange[1:1920 + 1920 * (fragNum - 1)],
         log = "xy", pch = 20, cex = 0.5,
         xlim = c(0.3, 5), ylim = c(0.3, 7),
         col = colors,
         las = 1,
         asp = 1,
         ylab = "", xlab = "", bty = "l",
         ...)
    abline(v = 2, lty = "dashed")
    segments(y0 = 2, x0 = 0.1, x1 = 2, lty = "dashed")
}

plotTF <- function(num, plate, row, col, name, colour) {
    points(medianAcrossScreens[384*(plate-1)+24*(row-1)+col],
           y1h$foldChange[num],
           pch = 20, cex = 0.5, col = colour)
    text(medianAcrossScreens[384*(plate-1)+24*(row-1)+col],
         y1h$foldChange[num],
         name, col = colour)
}

dir.create("figures")


# Used "darkred" and "darkblue" for slides.  Use lighter colors for poster:
splColor <- "red"
tcpColor <- "blue"


pdf("figures/2Dplot_SPL_TCP.pdf", width = 6.5, height = 6.5)

plot.2D(10) # fragment 10: pAGO7 part 2
plotTF(17668,2 ,  1,      4, "TCP2", tcpColor)
plotTF(17417         ,1   ,6     ,17,"TCP4", tcpColor)
plotTF(18102         ,3,   3      ,6, "TCP10", tcpColor)
plotTF(17284         ,1   ,1 , 4, "SPL4", splColor)
plotTF(17625         ,1  ,15      ,9, "SPL11", splColor)
plotTF(17371                ,1,   4     ,19, "SPL10", splColor)
plotTF(17346,        1,   3,     18, "SPL7", splColor)
## plotTF(18381,        3,  14,     21, "TCP18", tcpColor)

dev.off()

#---------------------------------------------------------------------
pdf("figures/TCP_SPL.pdf", width = 6.5, height = 6.5)

matplot(y1h$betagal[17281:19200], pch=20, log = "y", ylab = "",
        xaxp = c(0, 1920, 5), main = "",
        las = 1,
        xlab = "", bty = "l")
for (i in 46:50) {
    segments(x0 = (i - 46)*384 + 1, x1 = (i-45)*384, y0 = plateMedians[i])
    segments(x0 = (i - 46)*384 + 1, x1 = (i-45)*384,
             y0 = plateMedians[i] + 6*MADs[i], lty = "dashed")
}

plotTF2 <- function(num, name, color) {
    text(num - 17280, y1h$betagal[num],
         name, col = color)
    points(num - 17280, y1h$betagal[num],
           pch = 20, col = color)
}

plotTF2(17668,"TCP2", tcpColor)
plotTF2(17284, "SPL4", splColor)
plotTF2(17417,"TCP4", tcpColor)
plotTF2(17625, "SPL11", splColor)
plotTF2(18102, "TCP10", tcpColor)
plotTF2(17371, "SPL10", splColor)
plotTF2(17346, "SPL7", splColor)

dev.off()
#plotTF2(18381, "TCP18")


######################################################################

for (j in 1:3) {
    pdf(paste0("figures/supp-scatterplot-", j, ".pdf"),
        width = 6, height = 6)
    par(mfrow = c(2, 2),
        mar = c(5, 4, 2, 2) + 0.1)
    ##        c(5, 4, 4, 2) + 0.1)
    for (i in 1:4) {
        index <- i + (j - 1) * 4
        plot.2D(index, main = paste0(LETTERS[i], ". ",
                                 promoterFragment[index]),
                font.main = 1, cex.main = 1)
    }
    dev.off()
}

#' One by one:
plot.2D(1)
plot.2D(2)
plot.2D(3)
plot.2D(4)

plot.2D(5)
plot.2D(6)
plot.2D(7)
plot.2D(8)

plot.2D(9)
plot.2D(10)
plot.2D(11)
plot.2D(12)
