##################################################################
#
# Script to apply modelling functions to all pollen sites
#
##################################################################

suppressWarnings(library(tidyverse))
suppressWarnings(library(rio))
suppressWarnings(library(readxl))
suppressWarnings(library(broom))
suppressWarnings(library(vegan))
suppressWarnings(library(rioja))
source("PSIG_Functions.R")

# import all worksheets in file as a named list of data frames
# using function import_list in package rio
allpoll_list <- import_list("Woodbridge_et_al_2014_Data.xlsx")

# Convert to a tibble with 2 cols: Site_code and the list-column of data frames
# That is, the dataset for each site nested within the tibble
allpoll_nested <- tibble(Site_code=names(allpoll_list), data=allpoll_list) 

# Function to calculate rarefied pollen richness from original pollen data
# stored in our nested table.
fun_rare <- function(x) {
   x %>% select(!contains(non_pollen)) %>%
   rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%
   summarise(Depth, Age_BP, richness=rarefy(round(.[, -(1:2)]), 200)) %>%
   filter(Age_BP < 9000)
}

# apply fun_rare to our nested list using mutate / map
allpoll_rich <- allpoll_nested %>% 
  mutate(data, data=map(data, ~fun_rare(.x))) %>%
  unnest(data) %>%
  left_join(Site_info, by="Site_code") 

allpoll_rich
  
## Ignore warnings about count < 200

# Plot richness by region
ggplot(allpoll_rich, aes(Age_BP, richness, col=Site_code)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 20000, by=2000)) +
  facet_wrap(~Region) +
  labs(x="Years BP") +
  theme(legend.position="none")

# Question: How do affinity and richness co-vary?
# We can quantify / test the strength of relationship between 
# these two variables by calculating their correlation at each site.  
# There are a number of reasons why this is not the best way to answer 
# this question but it is useful for demonstration!
# We are also committing the schoolboy error of modelling before plotting 
# the data...

# To do this we need merge richness with the Affinity scores we calculated 
# in Part 2 (these are in df allpoll_LCC2) (you may need to source BES02_tidy.R to 
# generate this df if you get a not found error)
allpoll_comb <- allpoll_rich %>% 
  left_join(allpoll_LCC2, by=c("Site_code", "Depth", "Age_BP")) 

# Create a new column with abbreviated region names for plotting
allpoll_comb <- allpoll_comb %>%
  mutate(Region_abb = abbreviate(Region))

with(allpoll_comb, cor(Affinity, richness))

allpoll_corr <- allpoll_comb %>%
  group_by(Region_abb, Site_code) %>%
  summarise(r=cor(Affinity, richness), .groups="drop") 

allpoll_corr

# Is there regional variation in this relationship?
ggplot(allpoll_corr, aes(Region_abb, r)) +
   geom_boxplot()

# Which are the correlations significant?
# use cor.test instead of cor
if (0) {
  allpoll_corr2 <- allpoll_comb %>%
    group_by(Region_abb, Site_code) %>%
    summarise(r=cor.test(Affinity, richness), .group="drop") 
}

# Can only use summarise on functions that return a vector or 
# data frame, but cor.test returns a list 
# Therefore use map on a nested table
# An easy way to nest the data is to first use group_by, 
# then call nest()

allpoll_corr2 <- allpoll_comb %>%
  group_by(Region_abb, Site_code) %>%
  nest() %>%
  mutate(data, r=map(data, ~cor.test(.$Affinity, .$richness))) 

allpoll_corr2

################################
#
# Back to Slides
#
################################

with(allpoll_comb, cor.test(Affinity, richness))

with(allpoll_comb, tidy(cor.test(Affinity, richness)))

allpoll_corr2 <- allpoll_comb %>%
  group_by(Region_abb, Site_code) %>%
  nest() %>%
  mutate(data, r = map(data, ~tidy(cor.test(.$Affinity, .$richness)))) %>%
  unnest(r)

allpoll_corr2

# We really should plot these data to visualise the relationship 
# between richness and affinity
# allpoll_corr2 still has the original Affinity and richness data nested by 
# site. unnest these data and plot original affinity vs. richness, 
# grouped by region, but summarised for each site by a smooth

allpoll_corr3 <- allpoll_corr2 %>% unnest(data)

allpoll_corr3

ggplot(allpoll_corr3, aes(Affinity, richness, group=Site_code, 
                       col=p.value<0.05)) +
  geom_point(col="lightgrey") +
  geom_smooth(se=FALSE) +
  facet_wrap(~Region) +
  theme(legend.position="top")

# OR code by +ve / -ve relationship

ggplot(allpoll_corr3, aes(Affinity, richness, group=Site_code, 
                       col=estimate>0, linetype=p.value > 0.05)) +
  geom_point(col="lightgrey") +
  geom_smooth(se=FALSE) +
  facet_wrap(~Region) +
  theme(legend.position="top") + 
  scale_colour_discrete(name="r > 0")

# See Gavin's talk on Friday to fit GAMs to data like these
# Now you can fit models to multiple sites.


# How about trends in richness through time?

ggplot(allpoll_corr3, aes(Age_BP, richness, group=Site_code)) +
  geom_point(col="lightgrey") +
  geom_smooth(se=FALSE) +
  facet_wrap(~Region) +
  theme(legend.position="top") + 
  scale_colour_discrete(name="r > 0")

# One more thing...
# Sometimes we want to iterate over a list and call a function for its side-effects,
# that is, the function doesn't return new data but does something else, 
# like save a file or plot a figure.  
# Rather than use map to map our function fun to a list-column we use walk - 
# walk works in the same way as map but is just used to call fun for its side effects.

# We created a pollen diagram in Part 1.  What about creating diagrams for 
# all sites?  We could iterate round a loop but now we know better than that.
# We use poll_nested - this has original nested % data in wide format and 
# converts to % in long format using our existing fun_long_tidy function
# Then we create a new function (fun1) that take the % long data,
# removes taxa < 2% max, and transforms to wide format
# We then create a second function (fun2) that takes this wide data and creates 
# a stratigraphic diagram (we could combine fun1 and fun2 in one function, 
# but it is likely we could reuse fun1 elsewhere so we separate them)

fun_wide_5pc <- function(data, cut=0) {
  fun_max <- function(x, cut) max(x) > cut
  data %>% group_by(VarName) %>%
     filter(fun_max(Percent, cut)) %>% 
     pivot_wider(id_cols=c(Depth, Age_BP), names_from=VarName,
              values_from="Percent",
              values_fill=0) 
}

fun_strat_plot <- function(x, title="") {
  strat.plot(x[, -c(1:2)], x[, 1:2], scale.percent=TRUE, 
             y.rev=TRUE, title=title, yTop=0.7, cex.title=1.2,
             cex.xlabel=1)
}

# We then use mutate / map to apply our transformation functions, 
# and then use walk2 to call fun2 to plot the diagrams.
# walk2 is a variant of walk that takes 2 arguments and passes them 
# to the mapped function.  Here the first argument to walk is .x and 
# we supply the list-column of data, the second, .y is the Site_code 
# that we use to add a title to the plot.
# We surround the call with pdf() and dev.off() to send all graphic 
# output within the call to a pdf file.
if (0) {

  pdf("Pollen_plots.pdf", width=11, height=9, paper="a4r")

  allpoll_nested %>%
    mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
    mutate(data, data=map2(data, .y=2, .f=~fun_wide_5pc(.x, .y))) %>%
    walk2(.x=.$data, .y=.$Site_code, .f=~fun_strat_plot(.x, .y))
  
  dev.off()
}

# We could do something similar to save our cleaned percent-transformed data 
# to an Excel file.

if (0) {

  fun_save_excel <- function(x, y, wb) {
    addWorksheet(wb, sheetName=y)
    writeData(wb, y, x)
    addStyle(wb, y, style1, 2:(nrow(x)+1), 1:2, gridExpand = TRUE)
    addStyle(wb, y, style2, 2:(nrow(x)+1), 3:ncol(x), gridExpand = TRUE)
  }
  
  style1 <- openxlsx::createStyle(numFmt="0")
  style2 <- openxlsx::createStyle(numFmt="0.000")
  wb <- openxlsx::createWorkbook()
  
  allpoll_nested %>%
    mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
    mutate(data, data=map(data, ~fun_wide_5pc(.x))) %>%
    walk2(.x=.$data, .y=.$Site_code, .f=~fun_save_excel(.x, .y, wb))
  
  openxlsx::saveWorkbook(wb, "Woodbridge_Percent.xlsx", overwrite=TRUE)

# In this case there is a quicker way to export a list of dfs to a single Excel 
# file using function export in package rio.  Disadvantage of this approach 
# is that we cannot apply custom formatting to each sheet as we can with 
# write.xlsx in the package openxlsx

allpoll_PC_nested <- allpoll_nested %>%
    mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
    mutate(data, data=map(data, ~fun1(.x))) 

rio::export(allpoll_PC_nested$data, "Test.xlsx")

}

# Ever wondered where these sites are?
# A quick map using leaflet
if (0) {

library(leaflet)
  lab_opt <- labelOptions(noHide=TRUE, textsize="10px",
             direction="top", textOnly=TRUE,
             style=list(color="darkblue"))

  leaflet(data=Site_info) %>%
    addTiles() %>%
    addLabelOnlyMarkers(lng=~Longitude, lat=~Latitude, 
                      label=~Site_code, 
                      labelOptions=lab_opt) %>%
    addProviderTiles(providers$CartoDB.Positron)
}


# Your turn!!!

##########################
#
#  Back to Slides
#
#########################

# Let's put what we have learned together in a simple case study

# Task 1
# Using the UK lakes dataset, calculate the correlation between 
# log10-transformed Chlorophyll-a (Chla) and log10 total phosphorus (TP) for 
# the following lake depth classes: 0-2m 2-5m 5-10m and > 10 and fit regressions 
# of Chla ~ TP separately for each depth class.

UKLakes <- read_excel("UKLakes.xlsx")

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 

# Task 1: Reproduce the above scatter plot grouped by depth class
# classes = 0-2, 2-5, 5-10 and >10m
# Hint: you will need to create a new variable DepthClass using mutate
# and use cut to classify each observation into a depth class
# Then use facet_wrap on DepthClass

# Task 2: Use group_by and summarise to calculate the correlation 
# between log10-transformed Chla and TP for each depth class
# Hint: you can either use mutate to create new log10 variables or
# include log10 transformation in the call to cor

# Task 3: Use nest and mutate / map to calculate correlation & significance 
# of each depth class using cor.test

# Task 4: Use nest and mutate / map to fit a linear regression model 
# of log10 transformed Chla on TP for each Depth class


# Task 5: Back to our pollen data and examine rates of palynological 
# change through the Holocene.
# We can use the function fc_estimate_RoC in package RRatepol to estimate 
# rates of change.  Ondrej Mottl will talk about this package on Thursday.
# The code for running RRatepol is at the end of this script.  It takes 
# several minutes to run on a fast computer, substantially longer on a slow PC. 
# To save time I have pre-computer the results - load these results from 
# an excel file:

allpoll_roc <- read_excel("Woodbridge_ROC.xlsx")

# Task 5A: Plot the ROC data against age, split by region
# Hint age is in column Age, rate of change in column ROC
# There are some extreme values - you can use ylim to limit the range 
# of the y-axis, but ylim removes data, rather than truncates or clips the plot
# This is fine with point data but for lines it is better to use 
# coord_cartesian(ylim=c(0, 2)) to clip so no data is removed. 

# Task 5B: Aggregate the data by Region and 200-year Age slice and 
# calculate the mean and sd of each slice

###############################################################################

# Code for applying Rate of Change analysis using package RRatepol described in: 
# Ondřej Mottl et al. 2021. Rate-of-change analysis in palaeoecology revisited: 
# a new approach. https://doi.org/10.1101/2020.12.16.422943
# This function uses mutate() with across() to calculate column-wise sums to avoid
# converting to long and back.  It will take several minutes to run, hence 
# we use the pre-calculated data.

if (0) {

fun_ROC <- function(x) {
  x <- x %>% select(!contains(non_pollen)) %>%
     rename(age="Cal. yr. BP", depth="Depth (cm)") %>%
     rowwise() %>%
     mutate(total=sum(across(-c(depth, age)))) %>% 
     mutate(across(-total, ~.x/total*100)) %>%
     select(-total) %>%
     ungroup() %>% 
     filter(age <= 9000) %>%
     arrange(age) %>%
     mutate(sample_id=as.character(1:n())) %>%
     relocate(sample_id, .before=depth)
     roc <- fc_estimate_RoC(x[, -c(2:3)], x[, 1:3], 
                        Working_Units="bins", 
                        bin_size=200, treads=T, 
                        rand=1000) # reduce for testing to reduce time
  roc
}

allpoll_roc <- allpoll_nested %>% 
  mutate(data, data=map(data, ~fun_ROC(.x))) %>%
  unnest(data) %>%
  left_join(select(Site_info, Site_code, Region, Region2), 
            by="Site_code") 

# openxlsx::write.xlsx(poll_roc, "Woodbridge_ROC.xlsx")

}
