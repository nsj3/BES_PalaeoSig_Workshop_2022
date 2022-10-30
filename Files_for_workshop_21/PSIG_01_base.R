##################################################################
#
# Script to convert raw pollen counts to LCC classes using Base R
#
##################################################################

suppressWarnings(library(readxl))
suppressWarnings(library(rioja))
suppressWarnings(library(ggplot2))
suppressWarnings(library(kit))
source("PSIG_Functions.R")

# import pollen data from Red Mere, peat bog in English Fenland
# 48 levels spanning 700-4296 yr BP.

poll <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")

# remove unwanted columns (ie. those in non_pollen)
del <- colnames(poll) %in% non_pollen
poll <- poll[, !del]

# split data into depth/age & species parts
depth_age <- subset(poll, select=c(`Depth (cm)`, `Cal. yr. BP`))
spec <- subset(poll, select=-c(`Depth (cm)`, `Cal. yr. BP`))

#rename depth / age columns
colnames(depth_age) <- c("Depth", "Age_BP")

# transform counts to percentages
spec_pc <- spec / rowSums(spec) * 100

# transform to sqrt
spec_sqrt <- sqrt(spec_pc)

# aggregate columns into LCC groups
# create vector of LCC groups corresponding to columns names
# first, look up index of colnames in LCC taxon list
sel <- match(colnames(spec_sqrt), LCC_taxon_list$VarName)

sel

# extract LCC class for each taxon
LCC_code <- LCC_taxon_list$LCC_ID[sel]

LCC_code

# calculate column sums
# no easy way to calc grouped columns sums directly
# so we transpose data, calculate grouped rowsums then transpose 
# result back. 
spec_LCC <- t (rowsum(t(spec_sqrt), group=LCC_code))

# Yes, that's right, there are 2 functions to perform rowsums in base R!
# rowSums for non-grouped data and rowsum for grouped data. Obvious!

# transpose returns a matrix so we convert back to a data frame.  
# Colnames are numbers (LCC numeric classes)
# so we use check.names=FALSE to prevent R from prefixing 
# colnames with "x"
spec_LCC <- data.frame(spec_LCC, check.names=FALSE)

# we have sums of sqrt data, so normalise data
spec_LCC <- spec_LCC / rowSums(spec_LCC) * 100

spec_LCC

# Calculate the sum of tree (A = LCC classes 1-3) 
# and open (C = LCC classes) and  Affinity score (A-C)
spec_LCC$A <- rowSums(spec_LCC[, c("1", "2", "3")])  # arboreal classes
spec_LCC$C <- rowSums(spec_LCC[, c("5", "6", "7")])  # open classes
spec_LCC$Affinity <- spec_LCC$A - spec_LCC$C

# Create land cover class for each level in the core (LCC2)
# if Affinity > 20, LCC2 = max class max from 1-3 (woodland types)
# if Affinity < 20, LCC2 = max class  class from 5-7 (open types)
# if Affinity -20 to +20, if LCC = woodland type, then LCC2 = semi-open arboreal (8)
# if Affinity -20 to +20, if LCC = open type, then LCC2 = semi-open type (9-11)

# This is horrible to code in base R!  
# We could do it in a loop using switch() on each row, or vectorise it, 
# or use nested ifelse() but this is really horrible... 
# Help is on hand with a vectorised nested if function (nif) in package kit:
# BTW, this package has been around for 6 months, lucky us... 

# First we make vectors of the column names for our LCC groups 
# (A= 1-3 =arboreal; B = 4-6 = open, AC = both)
A_nms <- c("1", "2", "3")
C_nms <- c("5", "6", "7")
AC_nms <- c(A_nms, C_nms)

LCC2_ID <- kit::nif(
  spec_LCC$Affinity > 20, as.integer(A_nms[apply(spec_LCC[, A_nms], 1, which.max)]),
  spec_LCC$Affinity < -20, as.integer(C_nms[apply(spec_LCC[, C_nms], 1, which.max)]),
  default = ifelse(as.integer(AC_nms[apply(spec_LCC[, AC_nms], 1, 
                                           which.max)]) < 4, 
                   8L,
                   as.integer(AC_nms[apply(spec_LCC[, AC_nms], 1, which.max)]) + 4L)
)

# add the result to our df
spec_LCC$LCC2_ID <- LCC2_ID

# merge with depths
LCC2 <- cbind(depth_age, spec_LCC)

# merge with LCC2 lookup table to add LCC2 names to rows for plotting
LCC2 <- merge(LCC2, LCC2_assem_classes, by="LCC2_ID", all.x=TRUE, sort=FALSE)

# Output from merge is either sorts by the key (ie. LCC, which we don't want), or
# in an unspecified order (yes, the help really does say this).
# so re-sort data by Depth
LCC2 <- LCC2[order(LCC2$Depth), ]

LCC_names <- LCC_taxon_classes$LCC_name[as.integer(AC_nms)]

strat.plot(LCC2[, AC_nms], yvar=LCC2$Age_BP, 
         scale.percent=TRUE, y.rev=TRUE, x.names=LCC_names, 
         cex.xlabel=1, yTop=0.7)

ggplot(LCC2, aes(x=Age_BP, y=Affinity, col=LCC2_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=500)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top") +
  guides(col=guide_legend(nrow=2))

# Now your turn!!

#####################################################

# Before you start, a few guiding principles

# 1. Coding R is hard, lots of details to learn.
# 2. Don't get bogged down in the details of individual functions and operations.  
# 3. Focus on the big picture - the details will come with practice.

#####################################################

# Task 1
# Work through the above code using a different site.
# Function excel_sheets returns a vector of sheet names. Select another site, 
# edit line 16 to import this and work through the code with the new site.

if (0) {

sites <- excel_sheets("Woodbridge_et_al_2014_Data.xlsx")
print(sites)

}

# Task 2
# spec_sp contains % data.  We can use this to plot a full diagram:
if (0) {
   strat.plot(spec_pc[, -c(1:2)], yvar=depth_age$Age_BP, 
      scale.percent=TRUE, y.rev=TRUE, 
      cex.xlabel=0.8, yTop=0.7, plot.poly=TRUE,
      col.poly="darkgreen")
}

# Use whatever method you want, plot the stratigraphic diagram with taxa
# less that 2% max abundance omitted
# Hint: in base R I would use apply or sapply to calculate the max of 
# each column (ie. taxon) and then use a condition to either create a
# logical vector to index the required columns
# vector using a or a numeric vector 

