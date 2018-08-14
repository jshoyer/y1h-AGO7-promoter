
#' J. Steen Hoyer, Carrington lab, DDPSC.

#' To check these vectors in a tabular format, see
#' ../results/all-experiments/write-plate-layouts.R
#' and perhaps
#' ../vignettes/

## See also cellHTS2 getAlphaNumeric()

colbycol384 <- matrix(1:384, nrow = 16, ncol = 24)
rowbyrow384 <- matrix(1:384, nrow = 16, ncol = 24, byrow = TRUE)
#rowbyrow384[1:192]   ## Indexing is col-by-col

rowbyrow2colbycol <- as.vector(rowbyrow384)
colbycol2rowbyrow <- as.vector(t(matrix(1:384, nrow = 16, ncol = 24)))

#' Travel down a row, and then back along the next row, and then down the next,
#' i.e. the way the injector and detector travels,
#' i.e. A1, A2, ... A24, B24, B23, ... B1, C1, etc.
#' similar to 96-well version in
#' file:../results/all-experiments/barcoded-box1.R
downOneRowAndThenBack <- c(1:24, 48:25)
# This pattern repeats 8 times -- 8 pairs of rows, 2 * 24 wells per pair of rows
offsetsby2rows <- rep(0:7 * 48, each = 48)
rowandthenback <-
    matrix(downOneRowAndThenBack + offsetsby2rows, ncol = 24, byrow = TRUE)

#' Check with
#' options(width = 300)

colbycol2rowandthenback <- order(
    as.vector(rowandthenback)
    )
#' Can check with
#' cat(colbycol384[colbycol2rowandthenback], fill = 3)
#' but a more intuitive check with alphanumeric indices is below.

colbycol96 <- matrix(1:96, nrow = 8, ncol = 12)
rowbyrow96 <- matrix(1:96, nrow = 8, ncol = 12, byrow = TRUE)


### 384-well plates
                                        # Alphanumeric well IDs
                                        # 24 columns, 16 rows per plate
well384 <- paste0(rep(LETTERS[1:16], each=24) , sprintf("%02d", 1:24))
                                        # row-by-row order
mat384 <- data.frame(matrix(well384, nrow = 16, ncol = 24, byrow = TRUE))

well384byCol <- paste0(LETTERS[1:16],
                       rep(sprintf("%02d", 1:24), each = 16))

well384[rowbyrow2colbycol]
well384byCol[colbycol2rowbyrow]
well384byCol[colbycol2rowandthenback]

######################################################################
### 96-well plates
                                        # 12 columns, 8 rows per plate
well96 <- paste0(rep(LETTERS[1:8], each=12) , sprintf("%02d", 1:12))
mat96 <- data.frame(matrix(well96, nrow = 8, ncol = 12, byrow = TRUE))

well96byCol <- paste0(LETTERS[1:8],
                      rep(sprintf("%02d", 1:12), each = 8))


## Positions (column-by-column order) occupied by 96-well plates
## transferred into different quads of a 384-well plate
# A1: odd, odd
# A2: odd, even
# B2: even, even
# B1: even, odd
# A1 --> A2
#        |
#        V
# B1 <-- B2
oddRows <- seq(1, 15, 2)
oddCols <- seq(1, 23, 2)
quadA1 <- oddCols + rep((oddRows - 1) * 24, each = 12)
                                        # Row-by-row order
# A1  A3  A5
#  1   3   5
quadB1 <- quadA1 + 1
quadA2 <- quadA1 + 16
quadB2 <- quadA2 + 1

quadsA1andB1 <- rep(quadA1, each = 2) + c(0, 1)
quadsA1andA2 <- seq(1, 383, 2)

quadA1byCol <- rep((oddCols - 1) * 16, each = 8) + oddRows
                                        # Column-by-column order
quadB1byCol <- quadA1byCol + 1
quadA2byCol <- quadA1byCol + 16
quadB2byCol <- quadA2byCol + 1

######################################################################
#' # 384-well splots::plotScreen wrappers
columnarData <- 1:384
columnarData[colbycol2rowbyrow]


#' plotScreen help page: "values are assumed to come in row-by-row
#' order, e.g. A1, A2, A3, ..., B1, B2, ...".

### splots::plotScreen() expects a list of numerical vectors in
### row-by-row order.  Reformat column-by-column vector accordingly
### and plot.
splotWrapVector <-
    function(
        columnarData,  ##<< Down-each-column order
        ... ##<< Other arguments for plotScreen
        )
{
    require(splots)
    plotScreen(list(columnarData[colbycol2rowbyrow]),
               fill=c("white", "darkblue"), ncol = 1,
               ...)

    ## Returns nothing; used for side effect
}

## Example:
attr(splotWrapVector, 'ex') <- function() {
    originalOrder <- 1:384
    splotWrapVector(originalOrder)
}


### splots::plotScreen() expects a list of numerical vectors in
### row-by-row order.  Reformat a list of column-by-column vectors
### accordingly and plot.
splotWrapList <-
    function(
        columnarDataList,  ##<< Down-each-column order
        ... ##<< Other arguments for plotScreen
        )
{
    rowDataList <- lapply(columnarDataList,
                               function(i) i[colbycol2rowbyrow])
    require(splots)
    plotScreen(rowDataList,
               fill=c("white", "darkblue"),
               ...)

    ### Returns nothing; used for side effect
}

## Example:
## Looks best in a graphics device 3 * 24 / 18 = 4 times wider than tall
attr(splotWrapList, 'ex') <- function() {
    l <- list(1:384, 385:768, 769:1152)
    splotWrapList(l, ncol = 3)

    l <- list(c(1:383, NA), 385:768, 769:1152)
    splotWrapList(l, ncol = 3)
}
