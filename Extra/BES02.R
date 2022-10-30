library(tidyverse)
library(readxl)
library(rio)
theme_set(theme_bw(base_size=12))

LCC_taxon_list <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_List")
LCC_taxon_classes <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_Classes")
LCC2_assem_classes <- read_excel("LCC_Info.xlsx", sheet="LCC2_Assemblage_Classes")

d_nested <- import_list("Woodbridge_et_al_2014_Data.xlsx")

d <- tibble(Site_code=names(d_nested), data=d_nested) %>%
   tidyr::unnest(cols=-Site_code)

non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

d2 <- d %>% pivot_longer(cols=-c("Site_code", "Depth (cm)", "Cal. yr. BP"), 
                          names_to="VarName", values_to="Count") %>%
  filter(!is.na(Count) & !(VarName %in% non_pollen)) %>%
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`)

d3 <- d2 %>% 
  group_by(Site_code, Depth) %>%
  mutate(Percent = Count / sum(Count) * 100, 
         SQRT_PC = sqrt(Percent)) %>%
  ungroup()

d4 <- d3 %>% left_join(LCC_taxon_list, by="VarName") %>%
  left_join(LCC_taxon_classes, by="LCC_ID") %>%
  group_by(Site_code, Depth, Age_BP, LCC_ID, LCC_group) %>%
  summarise(Percent=sum(Percent), 
            SQRT_PC=sum(SQRT_PC), 
            .groups="drop") %>%
  group_by(Site_code, Depth, Age_BP) %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100,
         A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"])
         ) %>%
  ungroup()

d5 <- d4 %>%
  group_by(Site_code, Depth, Age_BP) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(Affinity=A-C, 
         LCC2_ID=case_when(
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
               TRUE ~ LCC_ID
               )
         )   

d6 <- d5 %>% filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC2_ID) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(names_from=LCC2_ID, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Age2, N), names_to="LCC2_ID", values_to="PC") %>%
  mutate(LCC2_ID=as.integer(LCC2_ID)) %>%
  left_join(LCC2_assem_classes, by="LCC2_ID")

ggplot(d6, aes(y=Age2, x=PC)) +
  geom_col(orientation="y", width=200, position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

LCC_order=c("Coniferous woodland",
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

d6$LCC2_name <- factor(d6$LCC2_name, levels=LCC_order)
ggplot(d6, aes(y=Age2, x=PC)) +
  geom_col(orientation="y", width=200, position=position_nudge(y=-100), fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")

