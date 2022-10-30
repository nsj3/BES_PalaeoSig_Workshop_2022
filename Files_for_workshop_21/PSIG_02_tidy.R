##################################################################
#
# Script to apply LCC to all sites using Tidyverse
#
##################################################################

suppressWarnings(library(tidyverse))
suppressWarnings(library(readxl))
suppressWarnings(library(rio))
source("PSIG_Functions.R")

# Create functions for our tidy code to convert pollen counts to LCC classes

# Test our tidy versions of the LCC functions on a single site 
poll <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")

poll_long_trans <- fun_long_tidy(poll, non_pollen=non_pollen)

poll_LCC2 <- fun_LCC2_tidy(poll_long_trans, LCC_taxon_list=LCC_taxon_list, 
                  LCC_taxon_classes=LCC_taxon_classes,
                  LCC2_assem_classes=LCC2_assem_classes)

ggplot(poll_LCC2, aes(x=Age_BP, y=Affinity, col=LCC2_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top") +
  guides(col=guide_legend(nrow=2))

# Now apply to all sites

##################################################################
# Back to Slides
##################################################################


# import all worksheets in file as a named list of data frames
# using function import_list in package rio
allpoll_list <- import_list("Woodbridge_et_al_2014_Data.xlsx")

# Convert to a tibble with 2 cols - Site_code and the list-column of data frames
# That is, the data frame for each site is nested within the tibble
allpoll_nested <- tibble(Site_code=names(allpoll_list), data=allpoll_list) 

allpoll_nested

allpoll_LCC2 <- allpoll_nested %>% 
  mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
  mutate(data, data=map(data, ~fun_LCC2_tidy(.x, LCC_taxon_list=LCC_taxon_list,
                                         LCC_taxon_classes=LCC_taxon_classes,
                                         LCC2_assem_classes=LCC2_assem_classes))) %>%
  unnest(data)


tmp0 <- allpoll_nested %>% 
  mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
  unnest(data)

stop()
library(data.table)

tmp <- lapply(1:500, function(x) tmp0)

tmp2 <- bind_rows(tmp)

system.time(
for (i in 1:10) {
tmp3 <- tmp2 %>% group_by(Site_code, Depth) %>%
  mutate(Perc = Count/sum(Count) * 100)
#  mutate(SS = sum(Count))
}
)

dt1 <- data.table(tmp2)
system.time(
for (i in 1:10) {
#tmp4 <- dt1[, Count/sum(Count)]
tmp4 <- dt1[, Perc := Count/sum(Count) * 100, by=list(Site_code, Depth)]
}
)

# Remove levels with Age_BP > 9000
# Create new variable, Age2, specifying the 200-year interval for each level 
# using mutate() with function cut_width()
# Then group by Age2 and LCC2 anc count the number of levels in each LCC2 at 
# each age slice using summarise(N=n())
# Then Normalise the count to 100%
# Our long data only contains positive counts for LCC2.  That is, LCC2 with zero 
# count for a particular Age slice are implied but not present in the data. 
# This will create problems when we plot the data as ggplot will "join the gaps"
# with zero counts.  A simple solution is to convert the data to wide format, and 
# filling implied zeros with explicit zero.  Then convert back to long.
# During this process the LCC2 codes ar convert to character so we turn them back 
# to integer and link to the lookup table of LCC2 class names
allpoll_agg <- allpoll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC2_ID) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(id_cols=c(Age2), names_from=LCC2_ID, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Age2), names_to="LCC2_ID", values_to="PC") %>%
  mutate(LCC2_ID=as.integer(LCC2_ID)) %>%
  left_join(LCC2_assem_classes, by="LCC2_ID")

# We then plot the data!
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

# The order of th LCC groups is slightly different from Woodbridge et al 2014, 
# so convert to factor and specify the order of the factor levels 
LCC2_order=c("Coniferous woodland",
          "Deciduous woodland",
          "Wet/fen woodland",
          "Heath",
          "Pasture/meadow",
          "Semi-open (arboreal)",
          "Semi-open (heath)",
          "Semi-open (pasture)",
          "Semi-open (arable)",
          "Arable",
          "Unused")

allpoll_agg$LCC2_name <- factor(allpoll_agg$LCC2_name, levels=LCC2_order)
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

# modify ggplot theme
my_theme <- theme_bw() + theme(
         panel.grid=element_blank(),
         strip.background=element_blank(),
         axis.line.x=element_line(),
         axis.line.y=element_line(),
         strip.text.x=element_text(angle=90, hjust=0, vjust=1)
)

ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP") + my_theme


# add geom_blank to force x-axis limits to minimum of 10
ggplot(allpoll_agg, aes(x=PC, y=Age2)) +
  geom_blank(aes(x=10)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP") + my_theme

# Now your turn!

# Task1: Produce a plot of land cover classes but showing England and 
# Scotland separately.
# Hint: Copy the aggregation code on line 72, join with the Site_info table, 
# and add Region2 to the grouping, then facet the plot by Region2.

# Task2: Create a figure showing the number of core levels per timeslice. 
# Hint, use summarise(N=n()) to count the number of rows

# Task3: Create a plot showing trends in Affinity scores for 0-9000 years 
# for for each site grouped by Region.

# Task4: Look at the elm decline.  Create a plot the %elm (Ulmus) against time
# grouped by region.
# Hint: Modify the code on line 43 so it only applies the function fun_long_tidy.
# This will create a long table of pollen counts and percentages.
# Filter this to extract Ulmus, and plot.

# Task5: Calculate the mean and sd affinity score for each time slice for 
# 0-9000 years and create plot showing overall mean affinity vs. Age with 
# error-bars showing standard deviation.

# Task6: Repeat the above plot, split by Region.
