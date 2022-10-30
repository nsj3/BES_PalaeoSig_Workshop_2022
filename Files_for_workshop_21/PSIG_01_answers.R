##################################################################
#
# Solutions for PSIG01_base.R and PSIG01_tidy.R
#
##################################################################

# Base R:

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
# logical vector to index the required columns or a numeric vector of 
# column numbers

sp_max <- apply(spec_pc, 2, max)
#or 
sp_max <- sapply(spec_pc, max)

# Then use the max values as a condition to subset columns
tmp <- spec_pc[, sp_max > 2]

strat.plot(tmp[, -c(1:2)], yvar=depth_age$Age_BP, 
    scale.percent=TRUE, y.rev=TRUE, 
    cex.xlabel=0.8, yTop=0.7, plot.poly=TRUE,
    col.poly="darkgreen")

############################################################

# Tidyverse

############################################################
# Task 1. Using tidy methods and function strat.plot as above, 
# reproduce the stratigraphic diagram of pollen types omitting taxa 
# with a maximum abundance of less than 2%
# Hint: poll_long_trans contains the % data.  Start with this and calculate 
# the max abundance of each taxon and use this to filter out those < 2%

tmp <- poll_long_trans %>%
  group_by(VarName) %>%
  mutate(max = max(Percent)) %>%
  filter(max > 2) %>% 
  pivot_wider(id_cols=c(Depth, Age_BP), 
      names_from=VarName, values_from=Percent, values_fill=0)

# A more concise and marginally faster but slightly less intuitive way is to create a function that calculates max value of a vector and compares to a cut value, (2%), then call that function in filter to select rows where max > 2.

fun <- function(x, cut) { max(x) > cut }
tmp <- poll_long_trans %>%
  group_by(VarName) %>%
  filter(fun(Percent, 2)) %>% 
  pivot_wider(id_cols=c(Depth, Age_BP), 
      names_from=VarName, values_from=Percent, values_fill=0)

strat.plot(tmp[,-c(1:2)], yvar=tmp$Age_BP, scale.percent=TRUE, 
           y.rev=TRUE, yTop=0.7, cex.title=1.2, cex.xlabel=0.8, 
           plot.poly=TRUE, col.poly="darkgreen", ylabel="Age (yr BP)", 
           srt.xlabel=45) 

# There is a way to accomplish the above without converting to long format
# But it is less intuitive and I don't know a simple way of summing across 
# LCC classes with the data in wide format. We have to use long format 
# form summing LCC classes so it is more efficient to use it for all 
# transformations.
# But if you are curious, here is how to convert to % and remove 
# taxa < 2% max without converting to long format
fun <- function(x, cut) { max(x) > cut }
poll2 <- poll %>% 
  select(!contains(non_pollen)) %>%
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  rowwise() %>%
  mutate(total=sum(across(-c(Depth, Age_BP)))) %>%
  mutate(across(-c(Depth, Age_BP, total), ~.x/total*100)) %>%
  select(-total) %>%
  select(where(~fun(., 2)) | contains(c("Depth", "Age_BP"))) %>%
  ungroup()

# Task 2: How does palynological richness vary with time?
# A simple measure is the number of taxa per level (T), but N depends on the 
# total pollen count, which also varies among levels.  So use rarefaction to 
# calculate the expected number of taxa (t) in a sample of fixed size (n) taken 
# from a larger sample (N) containing (T) taxa.
# Hint: use function rarefy in package vegan with a sample size (n) of 200
# Usage is rarefy(count_data, n)
# For this we can start with our raw count data in df poll.  
# Just remove unwanted columns, then use 
# summarise(Depth, Age_BP, rare=rarefy(.[, -c(1:2)], 200))
# this applies summarise to every row, and the .[, -(1:2)], passes a df of 
# each row, and removes cols 1 & 2 (ie. Depth and AgeBP before passing the 
# count data to rarefy).  Some pollen counts are decimals (N.5) due to 
# counting of fragments, so we surround the data argument with round() to 
# convert to integer.

poll_rich <- poll %>% select(!contains(non_pollen)) %>%
    rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%
    summarise(Depth, Age_BP, richness=rarefy(round(.[, -(1:2)]), 200)) 

poll_rich

# Ignore warning message about sample size less than 200

# Task 3: Plot the rarefied scores vs. Age using ggplot

ggplot(poll_rich, aes(Age_BP, richness)) + 
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=500)) +
  labs(x="Years BP") 

