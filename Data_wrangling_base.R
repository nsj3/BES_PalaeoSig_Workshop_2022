##################################################################
#
# Wrangling palaeodata with the Tidyverse 
# Part 1: Script to convert raw pollen counts to LCC classes and 
# aggregate data over all sites using Base R
#
##################################################################

library(readxl)
library(riojaPlot)

# import pollen data from Red Mere, peat bog in English Fenland
# 48 levels spanning 700-4296 yr BP.

polldata <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")
polldata

lcc_lookup <- read_excel("LCC_info.xlsx", sheet="LCC_Lookup")
lcc_lookup

# create character vector of non-pollen variables to remove
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

# remove unwanted columns (ie. those in non_pollen)
del <- colnames(polldata) %in% non_pollen
del
polldata <- polldata[, !del]

# split data into depth/age & species parts
depth_age <- subset(polldata, select=c(`Depth (cm)`, `Cal. yr. BP`))
poll_count <- subset(polldata, select=-c(`Depth (cm)`, `Cal. yr. BP`))

# Data in csv format?
d <- read.csv("Data/Redmere.csv")
colnames(d)
d <- read.csv("Data/Redmere.csv", check.names=FALSE)
colnames(d)
d <- readr::read_csv("Data/Redmere.csv", show_col_types = FALSE)
colnames(d)

#rename depth / age columns
colnames(depth_age) <- c("Depth", "Age_BP")

# transform counts to percentages
poll_pc <- poll_count / rowSums(poll_count) * 100

# transform to sqrt
poll_sqrt <- sqrt(poll_pc)

# aggregate columns into LCC groups
# First create a vector of LCC groups corresponding to taxon names
taxa_lcc <- merge(data.frame(VarName=colnames(poll_sqrt)), 
             lcc_lookup, by="VarName", all.x=TRUE, all.y=FALSE)

# Unfortunately merge sorts the resulting table alphabetically
# We need them in the original order so they match the order of the names
# in the pollen data
# Use match to 

m <- match(colnames(poll_sqrt), taxa_lcc$VarName)
taxa_lcc <- taxa_lcc[m, ]
head(taxa_lcc)
colnames(poll_sqrt)

# calculate column sums
# no easy way to calc grouped columns sums directly
# so we transpose data, calculate grouped rowsums then transpose 
# result back. 

poll_lcc <- t (rowsum(t(poll_sqrt), group=taxa_lcc$LCC_name))

# Yes, that's right, there are 2 functions to perform rowsums in base R!
# rowSums for non-grouped data and rowsum for grouped data. Obvious!

# Transpose returns a matrix so we convert back to a data frame.  
# so we use check.names=FALSE to prevent R from replacing spaces 
# in the names with underscore
poll_lcc <- data.frame(poll_lcc, check.names=FALSE)

# we have sums of sqrt data, so normalise data
poll_norm <- poll_lcc / rowSums(poll_lcc) * 100

riojaPlot::riojaPlot(poll_norm, depth_age, 
         scale.percent=TRUE,
         yvar.name="Age_BP")

dominant_class <- apply(poll_norm, 1, which.max)
poll_norm$lcc_class <- factor(colnames(poll_norm)[dominant_class])
poll_norm

# Calculate the sum of tree (A = LCC classes 1-3) 
# and open (C = LCC classes) and  Affinity score (A-C)
# Affinity gives an indication of "openness"
# High +ve Affinity = wooded, high -ve = open, 
# low +ve or -ve scores = semi-open
sum_arboreal <- rowSums(poll_norm[, c("Coniferous woodland", 
                                     "Deciduous woodland", 
                                     "Wet/fen woodland")])  # arboreal classes
sum_open <- rowSums(poll_norm[, c("Heath", 
                                     "Pasture/meadow", 
                                     "Arable indicators")])  # open classes
poll_norm$Affinity <- sum_arboreal - sum_open

poll_norm

# Woodbridge et al. use affinity to assign the land cover classes above to landscape types to 
# further split the woodland and pasture types into semi-open woodland and semi-open pasture
# For simplicity we will omit this step and just use the dominant land cover class for each 
# level in the core.

plot(depth_age$Age_BP, poll_norm$Affinity, type="l", 
     xlab="Age (years BP)",
     ylab="Affinity")
points(depth_age$Age_BP, poll_norm$Affinity, pch=19, cex=1, col=as.integer(poll_norm$lcc_class))
legend("topleft", pch=19, col=1:7, legend=levels(poll_norm$lcc_class))

# Now your turn!!

#####################################################

# Before you start, a few guiding principles

# 1. Coding R is hard, lots of details to learn.
# 2. Don't get bogged down in the details of individual functions and operations.  
# 3. Focus on the big picture - the details will come with practice.

#####################################################

# Task 1

# poll_sp contains full % data.  We can use this to plot a full diagram:

riojaPlot(poll_pc, depth_age, 
         scale.percent=TRUE,
         yvar.name="Age_BP")

# Use whatever method you want, plot the stratigraphic diagram with taxa
# less that 5% max abundance omitted
# Hint: in base R I would use apply or sapply to calculate the max of 
# each column (ie. taxon) and then use a condition to either create a
# logical vector or vector of names to index the required columns


##################################################################
#
# Wrangling palaeodata with the Tidyverse 
# Part 3: Apply base R code above to all 41 sites and aggregate 
# data in 200-year intervals
#
##################################################################

source("Data_wrangling_functions.R")

sites <- excel_sheets("Woodbridge_et_al_2014_Data.xlsx")

sites

LCC_list <- NULL
for (i in sites) {
  print(i)
  d <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet=i)
  tmp <- fun_LCC_base(d, lcc_lookup, non_pollen)
  tmp <- data.frame(Site_code=i, tmp)
  LCC_list <- rbind(LCC_list, tmp)   
}

LCC_list

# Create a new variable that classifies each data into a 200 year date range
cuts <- seq(0, 20000, by=200)
LCC_list$Age2 <- cut(LCC_list$Age_BP, breaks=cuts, labels=FALSE, include.lowest=TRUE) * 200 - 200

# use table() to cross-tab number of each LCC2 in each age slice
# and convert table to a data frame
LCC_count <- as.data.frame.matrix(with(LCC_list, table(Age2, lcc_class)))

# normalise numbers of LCC2 for each age slice
LCC_percent <- LCC_count / rowSums(LCC_count) * 100

# merge with Age2
LCC_percent <- data.frame(Age2=as.integer(rownames(LCC_count)), LCC_percent, 
                           check.names=FALSE)

# filter the data to exclude ages > 9000 
LCC_percent <- subset(LCC_percent, Age2 <= 9000)

# and plot...
yvar <- data.frame(AgeBP=LCC_percent$Age2+100)

riojaPlot(LCC_percent[, -1], yvar, 
         scale.percent=TRUE, 
         cex.xlabel=1, plot.bar=TRUE, 
         plot.line=FALSE, plot.poly=FALSE,
         lwd.bar=8, ymin=0, ymax=9000, yinterval=500,
         col.bar="black")


