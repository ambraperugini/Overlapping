rm(list=ls())
library(overlapping)
library(datawizard)

DATA<-read.csv("EngTurk.csv")


table(DATA$item[DATA$ItemType == "High_freq" & DATA$Language == "English"])
table(DATA$item[DATA$ItemType == "Low_freq" & DATA$Language == "English"])

low_frequency <-DATA$ReactionTime[DATA$item == "poor children"]
high_frequency <- DATA$ReactionTime[DATA$item == "blue eyes"]

xList <- list( x1 = low_frequency, x2 = high_frequency ) 
length(unlist(xList))

overlap( xList, plot = T )

TTEST <- with( xList, t.test( x1, x2, var.equal = TRUE ) )
TTEST

obsz <- 1 - overlap( xList )$OV
obsz
MX <- lapply(xList, mean)
SDX <- lapply(xList, var)
SKX <- lapply(xList, skewness)

B <- 2e3
n <- length(xList[[1]])
zperm <- sapply( 1:B, function(x){
  xperm <- sample( unlist( xList ) )
  xListperm <- list( x1 = xperm[1:n], x2 = xperm[(n+1):(n*2 )] )
  1 - overlap( xListperm )$OV
})


PVAL <- (sum( zperm > obsz )+1) / (length( zperm )+1)
PVAL

TTEST <- with( xList, t.test( x1, x2, var.equal = TRUE ) )
TTEST

## DIVERSE PROVE

low_frequency <-DATA$ReactionTime[DATA$item == "difficult life"]
high_frequency <- DATA$ReactionTime[DATA$item == "foreign business"]

# caso con significativo <.001 per-t ma non t-test
low_frequency <-DATA$ReactionTime[DATA$item == "great service"]
high_frequency <- DATA$ReactionTime[DATA$item == "blue eyes"]

low_frequency <-DATA$ReactionTime[DATA$item == "only friend"]
high_frequency <- DATA$ReactionTime[DATA$item == "blue eyes"]

# BINGO
low_frequency <-DATA$ReactionTime[DATA$item == "poor children"]
high_frequency <- DATA$ReactionTime[DATA$item == "blue eyes"]
