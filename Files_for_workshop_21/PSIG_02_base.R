##################################################################
#
# Script to apply LCC to all sites using Base R
#
##################################################################

suppressWarnings(library(tidyverse))
suppressWarnings(library(readxl))
suppressWarnings(library(rio))
suppressWarnings(library(rioja))
source("PSIG_Functions.R")

sites <- excel_sheets("Woodbridge_et_al_2014_Data.xlsx")

sites

LCC2 <- NULL

for (i in sites) {
  print(i)
  d <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet=i)
  tmp <- fun_LCC2_base(d, LCC_taxon_list=LCC_taxon_list, 
                        LCC_taxon_classes=LCC_taxon_classes)
  tmp <- data.frame(Site_code=i, tmp)
  LCC2 <- rbind(LCC2, tmp)   
}

# Create a new variable that classifies each data into a 200 year date range
cuts <- seq(0, 20000, by=200)
Age2 <- cut(LCC2$Age_BP, breaks=cuts, labels=FALSE, include.lowest=TRUE) * 200 - 200

# use table() to cross-tab number of each LCC2 in each age slice
LCC2_count <- with(LCC2, table(Age2, LCC2_ID))

# normalise numbers of LCC2 for each age slice
LCC2_percent <- LCC2_count[, -1] / rowSums(LCC2_count[, -1]) * 100

# convert table to a data frame
LCC2_percent <- as.data.frame.matrix(LCC2_percent)

# merge with Age2
LCC2_percent <- data.frame(Age2=as.integer(rownames(LCC2_count)), LCC2_percent, 
                           check.names=FALSE)

# filter the data to exclude ages > 9000 
LCC2_percent <- subset(LCC2_percent, Age2 <= 9000)

# create vector of LCC2 class names
LCC2_names <- LCC2_assem_classes$LCC2_name[as.integer(colnames(LCC2_percent)[-1])]

# and plot...
strat.plot(LCC2_percent[, -1], yvar=LCC2_percent$Age2+100, 
         scale.percent=TRUE, y.rev=TRUE, x.names=LCC2_names, 
         cex.xlabel=1, yTop=0.7, plot.bar=TRUE, plot.line=FALSE, 
         lwd.bar=8, ylim=c(0, 9000), y.tks=seq(0, 9000, by=500), col.bar="black")

