## Plot pollen data for Hockham Mere

library(readxl)
library(tidyverse)

poll <- read_excel("Woodbridge_et_al_2014_Data.xlsx",sheet="SIONASCA")
poll %>% print(n=4)

### It is possible to do MOST of what we need with the data in wide format.

### But it is not tidy, and tidyverse works best with tidy data (long format).

####   1. Code for tidy data is more intuitive 
####   2. I don't know a tidy easy way to group into LCC classes in wide format


## Plot pollen data for Hockham Mere

# create character vector of non-pollen variables to remove
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

poll_long <- poll %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  pivot_longer(cols=-c("Depth", "Age_BP"),                   
               names_to="VarName", values_to="Count") %>%    
  filter(!(VarName %in% non_pollen)) 

poll_long %>% print(n=6)

## Transform to percentages

poll_pc_long <- poll_long %>% 
  group_by(Depth) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  select(-Count) %>%
  ungroup() 

poll_pc_long

## Plot the diagram

# But we need it in wide format to plot!

poll_pc_wide <- poll_pc_long %>% 
  pivot_wider(names_from=VarName, values_from=Percent, values_fill=0)

# Use function strat.plot in my package `rioja`

library(rioja)

strat.plot(poll_pc_wide[, -c(1:2)], yvar=poll_pc_wide$Age_BP, y.rev=TRUE, 
  scale.percent=TRUE, ylabel="Age BP", plot.poly=TRUE, col.poly="darkgreen", 
  plot.bar=TRUE, col.poly.line=NA, srt.xlabel=45, y.tks=seq(1000, 11500, by=1000), 
  cex.xlabel=1, yTop=0.75, ylabPos=3, xLeft=0.09)

poll_pc_sub <- poll_pc_long %>%
  group_by(VarName) %>%
  filter(max(Percent) > 5) %>%
  pivot_wider(id_cols=c(Depth, Age_BP), names_from=VarName, values_from=Percent, 
              values_fill=0)

strat.plot(poll_pc_sub[, -c(1:2)], yvar=poll_pc_sub$Age_BP, y.rev=TRUE, 
  scale.percent=TRUE, ylabel="Age BP", plot.poly=TRUE, col.poly="darkgreen", 
  plot.bar=TRUE, col.poly.line=NA, srt.xlabel=45, y.tks=seq(1000, 11500, by=1000), 
  cex.xlabel=1, yTop=0.75, ylabPos=3, xLeft=0.09)

## Group taxa into landscape classes

# Join pollen % data to file of land classes 

poll_LCC <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_List")
poll_LCC %>% print(n=4)

# Names of the classes

LCC <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_Classes")
LCC %>% print(n=4)

## Group taxa into landscape classes II

# Join the tables

poll_pc_LCC <- poll_pc_long %>%
  left_join(poll_LCC, by="VarName") %>%
  left_join(LCC, by="LCC_ID") 
poll_pc_LCC %>% print(n=6)

## Group taxa into landscape classes III

# Sum percentages across land classes

poll_pc_LCC2 <- poll_pc_LCC %>%
  group_by(Depth, Age_BP, LCC_name) %>%
  summarise(Percent=sum(Percent)) %>%
  pivot_wider(id_cols=c(Depth, Age_BP), names_from=LCC_name, values_from=Percent, 
              values_fill=0)

strat.plot(poll_pc_LCC2[, -c(1:2)], yvar=poll_pc_LCC2$Age_BP, y.rev=TRUE, 
  scale.percent=TRUE, ylabel="Age BP", plot.poly=TRUE, col.poly="darkgreen", 
  plot.bar=TRUE, col.poly.line=NA, srt.xlabel=45, y.tks=seq(1000, 11500, by=1000), 
  cex.xlabel=1, yTop=0.75, ylabPos=3, xLeft=0.09)

## Reproducing Woodbridge *et al*. 2014

# Aggregating data for all 41 sites

# import all sheets from excel file in once go, data from each sheet is 
# an element in a list
library(rio)
allpoll <- import_list("Woodbridge_et_al_2014_Data.xlsx")
# bind_rows binds all the elements of the list together in a single table, 
# merging columns with the same name
allpoll2 <- bind_rows(allpoll, .id="Site")

# create character vector of non-pollen variables to remove
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

# remove convert to long format, remove non-pollen variables
allpoll_long <- allpoll2 %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  pivot_longer(cols=-c(Site, Depth, Age_BP),                   
               names_to="VarName", values_to="Count") %>%    
  filter(!(VarName %in% non_pollen)) %>%
  na.omit()

## Reproducing Woodbridge *et al*. 2014 II

# convert to percentage
allpoll_pc <- allpoll_long %>% 
  group_by(Site, Depth) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  select(-Count) %>%
  ungroup() 

# join pollen taxa with land cover classes and calculate % of each cover class
allpoll_pc_LCC <- allpoll_pc %>% 
  left_join(poll_LCC, by="VarName") %>%
  group_by(Site, Depth, Age_BP, LCC_ID) %>%
  summarise(Percent=sum(Percent), .groups="drop") 

# Use slice_max to select rows in each level with the highest % land cover class 
allpoll_pc_LCC2 <- allpoll_pc_LCC %>%
  group_by(Site, Depth, Age_BP) %>%
  slice_max(Percent, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  left_join(LCC, by="LCC_ID")

# filter to remove levels older than 9K, create new variable Age2 that is age bins 
# of 200 year intervals, calculate the number of levels in each interval in each 
# land cover class, convert these numbers to percentage of each cover class in 
# each 200 year interval
allpoll_agg <- allpoll_pc_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, 
                        labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC_name) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(id_cols=c(Age2), names_from=LCC_name, values_from=PC, 
              values_fill=0) %>%
  pivot_longer(cols=-c(Age2), names_to="LCC_name", values_to="PC")

## Reproducing Woodbridge *et al*. 2014 IV

# plot the data
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

unique(allpoll_agg$LCC_name)

# This does look exactly like the digram in woodbridge because they use 
# a slightly more complicated algorithm to allocate the levels in the NA column 
# to one of the other classes depending on the percent value.  

# But I think you get the idea how we can do a pretty sophisticated data 
# aggregation that would take weeks of manual site-by-site calculations in a 
# few minutes one we learn a few tidyverse functions.
