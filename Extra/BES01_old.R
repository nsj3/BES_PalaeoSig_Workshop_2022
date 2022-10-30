library(rio)
library(dplyr)
library(purrr)
library(readxl)
library(tidyr)
library(fuzzyjoin)
library(ggplot2)
library(pammtools)
library(forcats)
theme_set(theme_bw(base_size=12))

lst <- import_list("Test.xlsx")
LCC <- read_excel("LCC_Taxa_List.xlsx")
SiteInfo <- read_excel("SiteInfo.xlsx")

non_pollen <- c("Sample",
                "Depth [cm]",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

d <- tibble(Site_code=names(lst), data=lst) 

d2 <- d %>% tidyr::unnest(cols=-Site_code)

d3 <- d2 %>% pivot_longer(cols=-c("Site_code", "Depth (cm)", "Cal. yr. BP"), 
                          names_to="P_Var", values_to="Count") %>%
  filter(!is.na(Count) & !(P_Var %in% non_pollen)) %>%
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`)

d4 <- d3 %>% 
  group_by(Site_code, Depth) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  mutate(SQRT_PC = sqrt(Percent)) %>%
  ungroup()

d5 <- d4 %>% left_join(LCC, by=c("P_Var"="VarName")) %>%
  group_by(Site_code, Depth, Age_BP, LCC_ID, LCC_group) %>%
  summarise(SQRT_PC=sum(SQRT_PC), .groups="drop") %>%
  group_by(Site_code, Depth, Age_BP) %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100,
         A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
         Affinity=A-C, 
         LCC_ID2=case_when(
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
               TRUE ~ LCC_ID
               )
         ) %>%
  ungroup()

wc <- read_excel("Woodbridge_Classes.xlsx") %>%
  filter(Site_code %in% unique(d6$Site_code))

tmp <- bind_cols(d6, wc[, -3])

d6 <- d5 %>%
  group_by(Site_code, Depth, Age_BP) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  ungroup()

d7A <- tmp %>% filter(Age_BP < 9200) %>%
  mutate(Age2=cut_width(Age_BP2, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, Class) %>%
  summarise(PC=sum(Norm_SQRT_PC), .groups="drop_last") %>%
  mutate(PC=PC/sum(PC)*100) %>%
  pivot_wider(names_from=Class, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-Age2, names_to="Class", values_to="PC") %>%
  mutate(Class=as.integer(Class)) %>%
  left_join(LCC_assemblage, by=c("Class"="LCC_ID"))

d7A$LCC_name <- factor(d7A$LCC_name, levels=LCC_order)
ggplot(d7A, aes(Age2, ymin=0, ymax=PC)) +
  geom_stepribbon() +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  coord_flip() + 
  labs(y="%")

d8A <- tmp %>% filter(Age_BP < 9000) %>%
  mutate(Age2=cut_width(Age_BP2, boundary=-100, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, Class) %>%
  summarise(PC=n(), .groups="drop_last") %>%
  pivot_wider(names_from=Class, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-Age2, names_to="Class", values_to="PC") %>%
  mutate(PC=PC/sum(PC)*100, Class= as.integer(Class)) %>%
  left_join(LCC_assemblage, by=c("Class"="LCC_ID"))

d8A$LCC_name <- factor(d8A$LCC_name, levels=LCC_order)
ggplot(d8A, aes(Age2, ymin=0, ymax=PC)) +
  geom_stepribbon() +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  coord_flip() + 
  labs(y="%")



d7 <- d6 %>% filter(Age_BP < 9200) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC_ID2) %>%
  summarise(PC=sum(Norm_SQRT_PC), .groups="drop_last") %>%
  mutate(PC=PC/sum(PC)*100) %>%
  pivot_wider(names_from=LCC_ID2, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-Age2, names_to="LCC_ID2", values_to="PC") %>%
  mutate(LCC_ID2= as.integer(LCC_ID2)) %>%
  left_join(LCC_assemblage, by=c("LCC_ID2"="LCC_ID"))

d8 <- d6 %>% filter(Age_BP < 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=-100, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2, LCC_ID2) %>%
  summarise(PC=n(), .groups="drop_last") %>%
  mutate(PC=PC/sum(PC)*100, LCC_ID2= as.integer(LCC_ID2)) %>%
  left_join(LCC_assemblage, by=c("LCC_ID2"="LCC_ID"))
 
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

d7$LCC_name <- factor(d7$LCC_name, levels=LCC_order)
d7A$LCC_name <- factor(d7A$LCC_name, levels=LCC_order)
d8$LCC_name <- factor(d8$LCC_name, levels=LCC_order)
#d7$LCC_name <- sjmisc::word_wrap(d7$LCC_name, wrap=8)

ggplot(d7A, aes(Age2, ymin=0, ymax=PC)) +
  geom_stepribbon() +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  coord_flip() + 
  labs(y="%")

ggplot(d8, aes(Age2, ymin=0, ymax=PC)) +
  geom_stepribbon() +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  coord_flip() + 
  labs(y="%")

SiteInfo %>%
  distinct(Region)

SiteInfo <- SiteInfo %>% mutate(Region2=fct_collapse(Region,
                                        Scotland=c("East Scotland", 
                                                   "North Scotland", 
                                                   "South Scotland"), 
                                        England=c("Southeast England",
                                                  "Central England",
                                                  "Northeast England",
                                                  "Southwest England")
                                        )
                   )

d7 <- d6A %>% 
  left_join(SiteInfo, by="Site_code") %>%
  group_by(Region2, Age2, LCC_ID2) %>%
  summarise(PC=sum(Norm_SQRT_PC), .groups="drop_last") %>%
  mutate(PC=PC / sum(PC) * 100) %>%
  pivot_wider(names_from=LCC_ID2, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Region2, Age2), names_to="LCC_ID2", values_to="PC") %>%
  mutate(LCC_ID2= as.integer(LCC_ID2)) %>%
  left_join(LCC_assemblage, by=c("LCC_ID2"="LCC_ID"))

my_theme <- theme_bw  + theme(
         panel.grid=element_blank(),
         strip.background=element_blank(),
         panel.border=element_blank(),
         axis.line.x=element_line(),
         axis.line.y=element_line(),
         strip.text.x=element_text(angle=90, hjust=0, vjust=1)
)

ggplot(d7, aes(Age2, ymin=0, ymax=PC)) +
  geom_stepribbon() +
  facet_grid(Region2 ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() + 
  labs(y="%", x="Years BP") + my_theme 
