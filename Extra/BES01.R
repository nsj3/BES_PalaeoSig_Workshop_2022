library(tidyverse)
library(readxl)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(princurve)
library(rioja)
theme_set(theme_bw(base_size=12))

LCC_taxon_list <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_List")
LCC_taxon_classes <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_Classes")
LCC2_assem_classes <- read_excel("LCC_Info.xlsx", sheet="LCC2_Assemblage_Classes")

red <- read_excel("Woodbridge_et_al_2014_Data.xlsx", sheet="REDMERE")

non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

red_long <- red %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%
  pivot_longer(cols=-c("Depth", "Age_BP"), 
               names_to="VarName", values_to="Count") %>%
  filter(!(VarName %in% non_pollen))
  
red_long %>%
  distinct(VarName)

red_long_trans <- red_long %>% 
  group_by(Depth) %>%
  mutate(Percent = Count / sum(Count) * 100, 
         SQRT_PC = sqrt(Percent)) %>%
  ungroup() 

red_long_LCC <- red_long_trans %>% 
  left_join(LCC_taxon_list, by="VarName") %>%
  group_by(Depth, Age_BP, LCC_ID) %>%
  summarise(Percent=sum(Percent), 
            SQRT_PC=sum(SQRT_PC), 
            .groups="drop_last") %>%
  mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100) %>%
  ungroup() %>%
  left_join(LCC_taxon_classes, by="LCC_ID") 

ggplot(red_long_LCC, aes(y=Age_BP, x=Percent)) +
  geom_line(orientation="y") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 16000, by=1000)) +
  labs(x="%", y="Years BP")

red_long_LCC2 <- red_long_LCC %>%
  group_by(Depth, Age_BP) %>%
  mutate(A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
         C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
         Affinity=A-C, 
         LCC2_ID=case_when(
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
               Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
               TRUE ~ LCC_ID
               )
         ) %>%
  slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  left_join(LCC2_assem_classes, by="LCC2_ID")

p1 <- ggplot(red_long_LCC2, aes(x=Age_BP, y=Affinity, col=LCC2_name)) +
  geom_line(col="grey") +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP") + 
  scale_colour_discrete(name="LCC") +
  theme(legend.position="top")
p1

red_wide <- red3 %>% pivot_wider(id_cols=c(Depth, Age_BP), 
                                 names_from=VarName, values_from=Count, values_fill=0)

rare <- data.frame(Sum=rowSums(red_wide[, -c(1:4)]), 
rare=vegan::rarefy(red_wide[, -c(1:4)], 200))

red6 <- bind_cols(red5, rare)

p2 <- ggplot(red6, aes(x=Age_BP, y=rare)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP")  

p1 + p2 + plot_layout(ncol=1)

red_wide2 <- red3 %>% pivot_wider(id_cols=c(Depth, Age_BP), 
                                 names_from=VarName, values_from=Percent, values_fill=0)

pc1 <- princurve::principal_curve(x=sqrt(as.matrix(red_wide2[, -c(1:4)])), 
                                  trace=FALSE, penalty=1.4)

red6$PrCurve <- pc1$lambda

p3 <- ggplot(red6, aes(x=Age_BP, y=PrCurve)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP")  

p3 <- ggplot(red6, aes(x=Age_BP, y=PrCurve)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP")  

p1 + p2 + p3 + plot_layout(ncol=1)

red_wide2 <- red3 %>% pivot_wider(id_cols=c(Depth, Age_BP), 
                                 names_from=VarName, values_from=Percent, values_fill=0)
red_wide2 %>% {
  strat.plot(.[, -c(1:2)], yvar=.$Age_BP, 
             scale.percent=TRUE, y.rev=TRUE)
}

red3A <- red3 %>% group_by(VarName) %>%
  mutate(max = max(Percent)) %>%
  filter(max > 1) %>% 
  select(-max) %>%
  pivot_wider(id_cols=c(Depth, Age_BP), 
                                 names_from=VarName, values_from=Percent, values_fill=0)
  
red3A %>% {
  strat.plot(.[, -c(1:2)], yvar=.$Age_BP, 
             scale.percent=TRUE, y.rev=TRUE)
}

# save as excel file
