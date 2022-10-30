library(tidyverse)
library(readxl)
library(rio)
theme_set(theme_bw(base_size=12))

# import supporting lookup tables
LCC_taxon_list <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_List")
LCC_taxon_classes <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_Classes")
LCC2_assem_classes <- read_excel("LCC_Info.xlsx", sheet="LCC2_Assemblage_Classes")

# import all worksheets in file as a named list of data frames
poll_list <- import_list("Woodbridge_et_al_2014_Data.xlsx")

red <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")

tmp <- calc_long_trans(red)
tmp2 <- calc_LCC2(tmp, LCC_taxon_list=LCC_taxon_list, 
                  LCC_taxon_classes=LCC_taxon_classes,
                  LCC2_assem_classes=LCC2_assem_classes)

ggplot(red_long_LCC2, aes(x=Age_BP, y=Affinity, col=LCC2_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top")
p1

# Convert to a tibble with 2 cols - Site_code and the list of data frames
# That is, the data frame for each site is nested within the tibble
# unnest this to convert to a regular tibble with cols for each (taxon)
poll_wide <- tibble(Site_code=names(poll_list), data=poll_list) %>%
   tidyr::unnest(cols=-Site_code)

# character vector of columns we want to exclude
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

# Rename age and depth columns
# Convert to long format, using Site_code, Depth and Age_BP as our 
# indicator columns and VarName and Count as the key / variable columns
# Then filter to remove non-pollen variables
poll_long <- poll_wide %>% 
  pivot_longer(cols=-c("Site_code", "Depth (cm)", "Cal. yr. BP"), 
              names_to="VarName", values_to="Count") %>%
  filter(!is.na(Count) & !(VarName %in% non_pollen)) %>%
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`)


# Transform counts to percentage and sqrt percent
# Group by Site_code and Depth 
# Mutate creates a new variable Percent which is the count for each row, 
# divided by the sum of count for each depth (ie. each core level) * 100
# Then ungroup to tidy up
poll_trans <- poll_long %>% 
  group_by(Site_code, Depth) %>%
  mutate(Percent = Count / sum(Count) * 100, 
         SQRT_PC = sqrt(Percent)) %>%
  ungroup()

# Group taxa into LCC classes and sum sqrt_percent data over each class
# Merge the transformed data with the list of taxon codes to link to 
# the LCC class for each taxon using left_join
# Then group by Site_code Depth, Age and LCC class, and sum 
# Percent and SQRT_PC columns over each group
# Normalise the sqrt_percent column to 100% 
# Calculate the sum of tree (A = LCC classes 1-3) 
# and open (C = LCC classes) and  Affinity score (A-C)
# Then select the row with the maximum value of normalised 
# sqrt percent for each level using slice_max
# Finally ungroup to tidy up

poll_LCC <- poll_trans %>% left_join(LCC_taxon_list, by="VarName") %>%
  left_join(LCC_taxon_classes, by="LCC_ID") %>%
  group_by(Site_code, Depth, Age_BP, LCC_ID, LCC_group) %>%
  summarise(Percent=sum(Percent), 
            SQRT_PC=sum(SQRT_PC), 
            .groups="drop") %>%
  group_by(Site_code, Depth, Age_BP) %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100,
         A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]),
         Affinity=A-C) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  ungroup()

# Now calculate land cover class for each level (LCC2)
# if Affinity > 20, LCC2 = max class max from 1-3 (woodland types)
# if Affinity < 20, LCC2 = max class  class from 5-7 (open types)
# if Affinity -20 to +20, if LCC = woodland type, then LCC2 = semi-open arboreal (8)
# if Affinity -20 to +20, if LCC = open type, then LCC2 = semi-open type (9-11)
# Do this using the case_when function, which is dplyr's 
# vectorised / nested version of ifelse
# Finally join with lookup table to add LCC2 class names

poll_LCC2 <- poll_LCC %>%
  ungroup() %>%
  mutate(LCC2_ID=case_when(
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
               TRUE ~ LCC_ID
               )
         )   

# Remove levels with Age_BP > 9000
# Create new variable, Age2, specifying the 200-year interval for each level 
# using mutate() with function cut_width()
# Then group by Age2 and LCC2 anc count the number of levels in each LCC2 at 
# each age slice using summarise(n())
# Then Normalise the count to 100%
# Our long data only contains positive counts for LCC2.  That is, LCC2 with zero 
# count for a particular Age slice are implied but not present in the data. 
# This will create problems when we plot the data as ggplot will "join the gaps"
# with zero counts.  A simple solution is to convert the data to wide format, and 
# filling implied zeros with explicit zero.  Then convert back to long.
# During this process the LCC2 codes ar convert to character so we turn them back 
# to integer and link to the lookup table of LCC2 class names

poll_agg <- poll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC2_ID) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(names_from=LCC2_ID, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Age2, N), names_to="LCC2_ID", values_to="PC") %>%
  mutate(LCC2_ID=as.integer(LCC2_ID)) %>%
  left_join(LCC2_assem_classes, by="LCC2_ID")

# We then plot the data!

ggplot(poll_agg, aes(y=Age2, x=PC)) +
  geom_col(orientation="y", width=200, position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

stop()

# the order is slighty different from Woodbridge et al 2014, so reorder the class names

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

poll_agg$LCC2_name <- factor(poll_agg$LCC2_name, levels=LCC2_order)
ggplot(poll_agg, aes(y=Age2, x=PC)) +
  geom_col(orientation="y", width=200, position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

