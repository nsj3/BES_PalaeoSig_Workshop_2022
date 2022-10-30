##################################################################
#
# Script to convert raw pollen counts to LCC classes using Tidyverse
#
##################################################################

suppressWarnings(library(tidyverse))
suppressWarnings(library(readxl))
suppressWarnings(library(vegan))
suppressWarnings(library(rioja))
source("PSIG_Functions.R")

# import pollen data from Red Mere
poll <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")

# It is possible to do MOST of what we need with the data in wide format.
# But it is not tidy, and tidyverse works best with tidy data (long format).
# 1. Code for tidy data is more intuitive 
# 2. I don't know a tidy easy way to group into LCC classes in wide format

# Rename age and depth columns
# Convert to long format, using Depth and Age_BP as our indicator columns
# and VarName and Count as the key / variable columns
# Then filter to remove non-pollen variables
poll_long <- poll %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  pivot_longer(cols=-c("Depth", "Age_BP"),                   
               names_to="VarName", values_to="Count") %>%    
  filter(!(VarName %in% non_pollen))                         

# Transform counts to percentage and sqrt percent
# Group by Depth 
# Mutate creates a new variable Percent which is the count for each row, 
# divided by the sum of count for each depth (ie. each core level) * 100
# Then ungroup to tidy up

poll_long_trans <- poll_long %>% 
  group_by(Depth) %>%
  mutate(Percent = Count / sum(Count) * 100, 
         SQRT_PC = sqrt(Percent)) %>%
  ungroup() 

# Group taxa into LCC classes and sum sqrt_percent data over each class
# Merge the transformed data with the list of taxon codes to link to 
# the LCC class for each taxon using left_join
# Then group by Depth, Age and LCC class, and sum Percent and SQRT_PC 
# columns over each group
# Normalise the sqrt_percent column to 100% 
# Finally, join data to lookup with LCC class names 

poll_long_LCC <- poll_long_trans %>% 
  left_join(LCC_taxon_list, by="VarName") %>%
  group_by(Depth, Age_BP, LCC_ID) %>%
  summarise(Percent=sum(Percent), 
            SQRT_PC=sum(SQRT_PC), 
            .groups="drop_last") %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100) %>%
  ungroup() %>%
  left_join(LCC_taxon_classes, by="LCC_ID") 


# Plot the data for LCC classes
ggplot(poll_long_LCC, aes(y=Age_BP, x=Percent)) +
  geom_line(orientation="y") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 16000, by=1000)) +
  labs(x="%", y="Years BP")


# Calculate the sum of tree (A = LCC classes 1-3) 
# and open (C = LCC classes) and  Affinity score (A-C)
# this is just a summation over each level (ie. unique Depth / Age)
# Then select the row with the maximum value of normalised sqrt percent for 
# each level using slice_max

# Then calculate land cover class for each level (LCC2)
# if Affinity > 20, LCC2 = max class max from 1-3 (woodland types)
# if Affinity < 20, LCC2 = max class  class from 5-7 (open types)
# if Affinity -20 to +20, if LCC = woodland type, then LCC2 = semi-open arboreal (8)
# if Affinity -20 to +20, if LCC = open type, then LCC2 = semi-open type (9-11)
# Do this using the case_when function, which is dplyr's 
# vectorised / nested version of ifelse
# Finally join with lookup table to add LCC2 class names
poll_long_LCC2 <- poll_long_LCC %>%
  group_by(Depth, Age_BP) %>%
  mutate(A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
         Affinity=A-C) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(LCC2_ID=case_when(
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
               TRUE ~ LCC_ID
               )
         ) %>%
  left_join(LCC2_assem_classes, by="LCC2_ID")


# Plot the data - here Affinity score vs Age, with LCC2 classes
ggplot(poll_long_LCC2, aes(x=Age_BP, y=Affinity, col=LCC2_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=500)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top") +
  guides(col=guide_legend(nrow=2))

# But what about analysis that needs the data in wide format e.g. plotting 
# with rioja, ordination or rate of change analysis?
# We have a cleaned version of the percentage data in poll_long_trans 
# so just convert this to wide format:
poll_wide <- poll_long_trans %>% 
  pivot_wider(id_cols=c(Depth, Age_BP), 
  names_from=VarName, values_from=Percent, values_fill=0)

# Then plot...
strat.plot(poll_wide[,-c(1:2)], yvar=poll_wide$Age_BP, scale.percent=TRUE, 
           y.rev=TRUE, yTop=0.7, cex.title=1.2, cex.xlabel=0.8, 
           plot.poly=TRUE, col.poly="darkgreen", ylabel="Age (yr BP)", 
           srt.xlabel=45) 

# save Percent data as Excel file
if (0) {

   openxlsx::write.xlsx(poll_wide, "Red_Mere_Percent.xlsx")

}

# Your turn!!!

# Task 1. Using tidy methods and function strat.plot as above, 
# reproduce the stratigraphic diagram of pollen types omitting taxa 
# with a maximum abundance of less than 2%.
# Hint: poll_long_trans contains the % data.  Start with this and calculate 
# the max abundance of each taxon and use this to filter out those < 2%

# Task 2: How does palynological richness vary with time?
# A simple measure is the number of taxa per level (T), but N depends on the 
# total pollen count, which also varies among levels.  We therefore use 
# rarefaction to # calculate the expected number of taxa (t) in a 
# sample of fixed size (n) taken from a larger sample (N) containing (T) taxa.
# Hint: use function rarefy in package vegan with a sample size (n) of 200
# Usage is rarefy(count_data, n)

# For this we can start with our raw count data in df poll.  
# Just remove unwanted columns, then use:

# summarise(Depth, Age_BP, richness=rarefy(round(.[, -c(1:2)]), 200))

# This applies summarise to every row, and the .[, -(1:2)] passes each 
# row and removes cols 1 & 2 (ie. Depth and AgeBP before passing the 
# count data to rarefy), Some pollen counts are decimals (*.5) due to 
# counting of fragments, so we surround the data argument with round() to 
# convert to integer.

# Task 3: Plot the rarefied scores vs. Age using ggplot


