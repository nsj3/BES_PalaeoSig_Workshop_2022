##################################################################
#
# Wrangling palaeodata with the Tidyverse 
# Part 2: Script to convert raw pollen counts to LCC classes 
# using Tidyverse
#
##################################################################

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(readxl)
library(riojaPlot)

# import pollen data from Red Mere
polldata <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")
lcc_lookup <- read_excel("LCC_info.xlsx", sheet="LCC_Lookup")

# It is possible to do MOST of what we need with the data in wide format.
# But it is not tidy, and tidyverse works best with tidy data (long format).
# 1. Code for tidy data is more intuitive 
# 2. I don't know a tidy easy way to group into LCC classes in wide format

# create character vector of non-pollen variables to remove
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

# Rename age and depth columns
# Convert to long format, using Depth and Age_BP as our indicator columns
# and VarName and Count as the key / variable columns
# Then filter to remove non-pollen variables
polllong <- polldata %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  pivot_longer(cols=-c("Depth", "Age_BP"),                   
               names_to="VarName", values_to="Count") %>%    
  filter(!(VarName %in% non_pollen))                         

# Transform counts to percentage and sqrt percent
# Group by Depth 
# Mutate creates a new variable Percent which is the count for each row, 
# divided by the sum of count for each depth (ie. each core level) * 100
# Then ungroup to tidy up

polllong_sqrt <- polllong %>% 
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

polllong_norm <- polllong_sqrt %>% 
  left_join(lcc_lookup, by="VarName") %>%
  group_by(Depth, Age_BP, LCC_name, LCC_group) %>%
  summarise(SQRT_PC=sum(SQRT_PC), .groups="drop") %>%
  group_by(Depth, Age_BP) %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100) %>%
  ungroup() %>%
  select(-SQRT_PC) 

polllong_norm

# Plot the data for LCC classes
theme_set(theme_bw(base_size=16))
ggplot(polllong_norm, aes(y=Age_BP, x=Norm_SQRT_PC)) +
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

polllong_lcc <- polllong_norm %>%
  group_by(Depth, Age_BP) %>%
  mutate(A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
         Affinity=A-C) %>%
  select(-c(A, C)) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  select(-c(Norm_SQRT_PC)) %>%
  ungroup() 

# Plot the data - here Affinity score vs Age, with LCC2 classes
ggplot(polllong_lcc, aes(x=Age_BP, y=Affinity, col=LCC_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=500)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top") +
  guides(col=guide_legend(nrow=2))

# But what about analysis that needs the data in wide format e.g. plotting 
# with rioja, ordination or rate of change analysis?
# We have a cleaned version of the percentage data (and sqrt data) in polltidy_sqrt
# so just convert this to wide format:
poll_wide <- polllong_sqrt %>% 
  pivot_wider(id_cols=c(Depth, Age_BP), 
  names_from=VarName, values_from=Percent, values_fill=0)

# Then plot...
riojaPlot(poll_wide[,-c(1:2)], poll_wide[, 1:2], 
          scale.percent=TRUE, 
          srt.xlabel=45) 

# save Percent data as Excel file
if (0) {
   openxlsx::write.xlsx(poll_wide, "Red_Mere_Percent.xlsx")
}

# Your turn!!!

# Task 1. Using tidy methods and function riojaPlot 
# reproduce the stratigraphic diagram of pollen types omitting taxa 
# with a maximum abundance of less than 5%.
# Hint: polllong_sqrt contains the % data.  Start with this and calculate 
# the max abundance of each taxon and use this to filter out those < 5%


##################################################################
#
# Wrangling palaeodata with the Tidyverse 
# Part 4: Script to apply LCC to all sites and aggregate results at
# 200-year intervals using Tidyverse
#
##################################################################

source("Data_wrangling_functions.R")

allpoll_list <- rio::import_list("Woodbridge_et_al_2014_Data.xlsx")

allpoll_nested <- tibble(Site=names(allpoll_list), polldata=allpoll_list) 

allpoll_long <- allpoll_nested %>% 
  mutate(polldata, polldata=map(polldata, ~ fun_long_tidy(.x, non_pollen=non_pollen))) %>%
  unnest(cols=polldata)

allpoll_long

allpoll_lcc <- allpoll_long %>% 
  group_by(Site, Depth, Age_BP) %>%
  mutate(Percent = Count / sum(Count) * 100, 
         SQRT_PC = sqrt(Percent)) %>%
  ungroup()  %>% 
  left_join(lcc_lookup, by="VarName") %>%
  group_by(Site, Depth, Age_BP, LCC_name, LCC_group) %>%
  summarise(SQRT_PC=sum(SQRT_PC), .groups="drop") %>%
  group_by(Site, Depth, Age_BP) %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100) %>%
  ungroup() %>%
  select(-SQRT_PC)  %>%
  group_by(Site, Depth, Age_BP) %>%
  mutate(A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
         Affinity=A-C) %>%
  select(-c(A, C)) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  select(-c(Norm_SQRT_PC)) %>%
  ungroup() 

allpoll_agg <- allpoll_lcc %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC_name) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(id_cols=c(Age2), names_from=LCC_name, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Age2), names_to="LCC_name", values_to="PC") 

# We then plot the data!
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP") +
  theme(strip.text=element_text(angle=90, size=10, hjust=0, vjust=1))

# Add a blank geom to force graphs to be at least 10 units wide

tmp <- data.frame(PC=10, Age2=2000)
# We then plot the data!
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP") +
  theme(strip.text=element_text(angle=90, size=10, hjust=0, vjust=1)) +
  geom_blank(data=tmp)

##################################################################
#
# Wrangling palaeodata with the Tidyverse 
# Part 5: Script to apply modelling functions to all pollen sites
# Calculating palynological richness
#
##################################################################

# How does palynological richness vary with time?
# A simple measure is the number of taxa per level (T), but N depends on the 
# total pollen count, which also varies among levels.  We therefore use 
# rarefaction to calculate the expected number of taxa (t) in a 
# sample of fixed size (n) taken from a larger sample (N) containing (T) taxa.

# The rarefy function in the vegan package takes a matrix or data frame of taxon
# counts and returns a vector of rarefied number of taxa.

# We wrap the rarefy function in a wrapper function that also removes non-pollen
# variables and renames columns of the original data, as before.

library(vegan) # for rarefy function

fun_rare <- function(x) {
   x %>% select(!contains(non_pollen)) %>%
   rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%
   summarise(Depth, Age_BP, richness=rarefy(round(.[, -(1:2)]), 200)) %>%
   filter(Age_BP < 9000)
}

# apply fun_rare to our nested list using mutate / map
rich <- allpoll_nested %>% 
  mutate(polldata, rare=map(polldata, ~fun_rare(.x))) %>%
  unnest(rare) %>%
  select(-polldata)

## Ignore warnings about count < 200

rich
  
ggplot(rich, aes(Age_BP, richness, col=Site)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 20000, by=2000)) +
  theme(legend.position="none")

site_info <- read_excel("LCC_info.xlsx", sheet="Site_Info")

rich2 <- rich %>%
  left_join(site_info, by=c("Site"="Site_code"))

# Plot richness by region
ggplot(rich2, aes(Age_BP, richness, col=Site)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 20000, by=2000)) +
  facet_wrap(~Region) +
  labs(x="Years BP") +
  theme(legend.position="none")

# Add a regression line to identify trends
ggplot(rich2, aes(Age_BP, richness, col=Site)) +
  geom_line() +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_continuous(breaks=seq(0, 20000, by=2000)) +
  facet_wrap(~Region) +
  labs(x="Years BP") +
  theme(legend.position="none")

# Which trends are significant?
# To answer this we need to perform separate regressions using lm and 
# extract the slopes and p-vaues
# A job for summarise on grouped data:
# First we create a new age variable (that )-ve Age_BP) so time runs forward

rich3 <- rich2 %>%
  mutate(neg_Age=-Age_BP) 

if (0) {
rich3 %>%
  group_by(Site, Region) %>%
  summarise(lm(richness ~ neg_Age), .groups="drop") 
}

RoC <- rich3 %>%
  group_by(Site, Region) %>%
  summarise(broom::tidy(lm(richness ~ neg_Age)), .groups="drop") 

RoC <- rich3 %>%
  group_by(Site, Region) %>%
  summarise(broom::tidy(lm(richness ~ neg_Age)), .groups="drop") %>% 
  filter(term=="neg_Age")

RoC

ggplot(RoC, aes(Region, estimate)) +
  geom_boxplot() +
  coord_flip()

# Join with site info and sort region by latitude
RoC %>% left_join(site_info, by=c("Region", "Site"="Site_code"))  %>%
  mutate(Region=fct_reorder(Region, Latitude, mean)) %>%
  ggplot(aes(Region, estimate)) +
     geom_boxplot() +
     coord_flip()

# Add dotplot to show data points (with small dataset boxes can be misleading)

RoC %>% left_join(site_info, by=c("Region", "Site"="Site_code"))  %>%
  mutate(Region=fct_reorder(Region, Latitude, mean)) %>%
  ggplot(aes(Region, estimate)) +
     geom_boxplot(col="grey") +
     geom_dotplot(stackdir="center", binaxis="y", dotsize=0.4) +
     coord_flip()

# plot rate of change of palynological richness vs latitude
RoC %>% 
  left_join(site_info, by=c("Site"="Site_code")) %>%
  ggplot(aes(Latitude, estimate * 1000)) +
  geom_point() +
  geom_smooth() +
  labs(y="Rate of change of pollen richness\n(taxa per 1000 years)")

# Is relationship significant?
# Crude test using rank correlation

RoC %>% 
  left_join(site_info, by=c("Site"="Site_code")) %>%
   summarise(broom::tidy(cor.test(estimate, Latitude, method="spearman")))

# How to save multiple diagrams?

fun_plot <- function(x, y) {
   print(y)
   x <- x %>% select(!contains(non_pollen)) %>%
   rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%
   filter(Age_BP < 9000) 
   spec <- x[, -c(1:2)]
   spec <- spec / rowSums(spec) * 100
   spec <- spec[, apply(spec, 2, max) > 2]
   ymin <- min(x$Age_BP)
   ymax <- max(x$Age_BP)
   yinterval <- 200
   ytks1 <- seq(0, 9000, by=200)
   riojaPlot(spec, x[, 1:2], 
             yvar.name="Age_BP",
             plot.groups=TRUE,
             plot.cumul=TRUE,
             scale.percent=TRUE, ymin=ymin, ymax=ymax, 
             yinterval=200, ytks1=ytks1,
             srt.xlabel=45)
   mtext(y, side=3, outer=F, line=3, adj=0)
}

# We then use walk2 to call fun_plot to plot the diagrams.
# walk2 is a variant of walk that takes 2 arguments and passes them 
# to the mapped function.  Here the first argument to walk is .x and 
# we supply the list-column of data, the second, .y is the Site name 
# that we use to add a title to the plot.
# We surround the call with pdf() and dev.off() to send all graphic 
# output within the call to a pdf file.
if (0) {

  pdf("Pollen_plots.pdf", width=11, height=9, paper="a4r")
  allpoll_nested %>%
    walk2(.x=.$polldata, .y=.$Site, .f=~fun_plot(.x, .y))
  
  dev.off()
}

##########################################################################
#
#  Your Turn!
#
##########################################################################

# Task 3
# How does richness vary with Affinity (degree of Openness of landscape)?
#
# Join the table of palynolical richness (rich3) to the Affinity scores (allpoll_lcc) and 
# plot richness against affinity.  Is there a relationship for the whole dataset?
# Is the relationship different for different sites.  Is there a regional pattern in the relationship?


# Regress richness against Affinity for each site.  Does the slope of this relationship 
# vary by region and latitude?

