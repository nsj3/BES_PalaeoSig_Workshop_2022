library(riojaExtra)

#d <- read.CEP("Perc146.cep")
#d <- round(d, 3)
#d <- d / rowSums(d) * 100
#openxlsx::write.xlsx(d, "UKLakes_sp.xlsx", row.names=TRUE)

chem <- read_excel("UKLakes.xlsx")
diat <- read_excel("UKLakes.xlsx", sheet = "Diatoms")
taxon_list <- read_excel("UKLakes.xlsx", sheet = "Lifeform")

diat_l <- diat %>% 
  pivot_longer(cols=-SiteCode, names_to="Taxon", values_to="Percent")


codes <- diat_l %>% distinct(Taxon) %>%
        filter(str_sub(Taxon, 1, 2) %in% 
           c("AS", "AU", "CY", "CC", "ST", "TH") | Taxon == "FR008A")

plank <- diat_l %>% distinct(Taxon) %>%
     mutate(Planktic = case_when(
                Taxon %in% codes$Taxon ~ TRUE,
                TRUE ~ FALSE))

openxlsx::write.xlsx(plank, "tmp.xlsx")

n <- apply(diat[, -1]>0, 1, sum)

nTaxa <- data.frame(SiteCode=diat[, 1], N=n)

chem2 <- merge(chem, nTaxa, by="SiteCode", all.x=TRUE) %>%
  mutate(DepthClass = cut(MaxDepth, breaks=c(0, 2, 5, 10, 100))) %>%
  mutate(TPClass = cut_number(TP, 4)) #, breaks=c(0, 30, 100, 200, 2000))) 
  

ggplot(chem2, aes(TP, N)) +
  geom_point() +
  scale_x_log10() + 
  geom_smooth(method="lm") +
  facet_wrap(~DepthClass)

ggplot(chem2, aes(MaxDepth, N)) +
  geom_point() +
  scale_x_log10() + 
  geom_smooth(method="lm") +
  facet_wrap(~TPClass)



dd <- diat_l %>% 
  left_join(taxon_list, by="Taxon") %>%
  group_by(SiteCode, Lifeform) %>%
  summarise(Percent=sum(Percent)) %>%
  pivot_wider(id_cols=SiteCode, names_from=Lifeform, values_from=Percent) %>%
  right_join(chem2) 

ggplot(dd, aes(MaxDepth, Planktic)) +
  geom_point() +
  scale_x_log10() + 
  geom_smooth() 

ggplot(dd, aes(MaxDepth, Planktic)) +
  geom_point() +
  scale_x_log10() + 
  geom_smooth() +
  facet_wrap(~TPClass)

ggplot(dd, aes(TP, Planktic)) +
  geom_point() +
  scale_x_log10() + 
  geom_smooth(method="lm") +
  facet_wrap(~DepthClass)



+
#  facet_wrap(~DepthClass)

