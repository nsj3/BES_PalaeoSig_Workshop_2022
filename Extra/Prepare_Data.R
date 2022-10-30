library(rio)
library(dplyr)
library(purrr)
library(fuzzyjoin)
library(openxlsx)

SiteInfo <- read_excel("SiteInfo.xlsx")

SiteInfo <- SiteInfo %>% mutate(Region2=fct_collapse(Region,
                                        Scotland=c("East Scotland", 
                                                   "North Scotland", 
                                                   "South Scotland"), 
                                        England=c("Southeast England",
                                                  "Central England",
                                                  "Northeast England",
                                                  "Southwest England")))
write.xlsx(SiteInfo, "Site_Info.xlsx")

paths <- Sys.glob("Data/*.xls")
lst <- import_list(paths, setclass="tbl_df", skip=1, sheet="Raw", 
                   rbind=F, rbind_label="Dataset")
names(lst) <- gsub(" PFTs", "", names(lst))
names(lst) <- gsub(" short chron", "", names(lst))
names(lst)
depth <- c("Depth (cm)",
          "Depth [cm]")

d <- tibble(Site_code=names(lst), data=lst) %>%
   mutate(data = map(data, ~rename_with(.x, 
              function(x) c("Depth (cm)"), starts_with("Depth")))) %>%
   mutate(data = map(data, ~rename_with(.x, 
              function(x) c("Sum"), starts_with("...")))) 
   mutate()

SiteInfo <- read_excel("SiteInfo.xlsx")

d <- filter(d, Site_code %in% SiteInfo$`Site_code`)

export(d$data, "Woodbridge_et_al_2014_Data.xlsx")

LCC_taxa <- data.frame(LCC_ID=1:7, LCC_name=c(
  "Coniferous woodland",
  "Deciduous woodland",
  "Wet/fen woodland",
  "Unused",
  "Heath",
  "Pasture/meadow",
  "Arable indicators"),
  LCC_group=c("A", "A", "A", "B", "C", "C", "C")
  )

LCC_assemblage <- data.frame(LCC_ID=c(1:11), 
                            LCC_name=c(
                              "Coniferous woodland",
                              "Deciduous woodland",
                              "Wet/fen woodland",
                              "Unused",
                              "Heath",
                              "Pasture/meadow",
                              "Arable",
                              "Semi-open (arboreal)",
                              "Semi-open (heath)",
                              "Semi-open (pasture)",
                              "Semi-open (arable)"
                           )
                  )

lst <- import_list(paths, setclass="tbl_df", sheet="Raw", 
                   rbind=F, rbind_label="Dataset")
names(lst) <- gsub(" PFTs", "", names(lst))
names(lst) <- gsub(" short chron", "", names(lst))

fun <- function(x) {
  data.frame(LCC_ID=colnames(x), VarName=as.character(x[1, ])) %>%
    filter(regexpr("^[1-7]", LCC_ID) > 0) %>%
    mutate(LCC_ID=as.numeric(substring(LCC_ID, 1, 1)))
}
  
dX <- tibble(Site_code=names(lst), data=lst) %>%
   filter(Site_code %in% SiteInfo$`Site_code`) %>%
   mutate(data = map(data, ~fun(.x))) %>%
   tidyr::unnest(cols=-Site_code) %>% 
   select(-Site_code) %>%
#   group_by(VarName, LCC_ID) %>%
#   count() %>%
   group_by(VarName) %>%
   summarise(LCC_ID=as.integer(round(mean(LCC_ID), 0))) %>%
   left_join(LCC_taxa, by="LCC_ID")


p_vars <- read_excel("EMPD_pft_P_Vars.xlsx") %>%
  select(VarName, GroupId) %>%
  rename("VarName2" = "VarName")

w <- c(d=0.05, i=0.2, s=.5, t=1)
nms3 <- stringdist_left_join(dX, p_vars, by=c("VarName"="VarName2"), 
                             weight=w, max_dist=5, distance_col="Dist") %>%
  group_by(VarName, LCC_ID, LCC_name) %>%
  slice_min(Dist, n=1, with_ties=FALSE) %>%
  rename("P_Group"="GroupId") %>%
  ungroup()

pg <- read_excel("P_GROUPS.xlsx")
#nms3 %>% distinct(P_Group) %>% left_join(pg, by=c("P_Group"="GroupId"))

#export(dX, "LCC_Taxa_List.xlsx")
wb <- openxlsx::createWorkbook()
addWorksheet(wb, "LCC_Taxon_Lookup")
addWorksheet(wb, "LCC_Taxon_Classes")
addWorksheet(wb, "LCC_Assemblage_Classes")
addWorksheet(wb, "P_Groups")
writeData(wb, "LCC_Taxon_Lookup", nms3)
writeData(wb, "LCC_Taxon_Classes", LCC_taxa)
writeData(wb, "LCC_Assemblage_Classes", LCC_assemblage)
writeData(wb, "P_Groups", pg)
#saveWorkbook(wb, "LCC_Info.xlsx", overwrite = TRUE)

rm(lst, wb, fun, dX, d, paths, depth, LCC_assemblage, LCC_taxa, SiteInfo, pg, nms3)

paths <- Sys.glob("Data/*.xls")
lst <- import_list(paths, setclass="tbl_df", skip=0, sheet="Affinity scores graph", 
                 rbind=F, rbind_label="Dataset")
f <- function(x) {
  grep("Cal. yr. BP", colnames(x))
}

tt <- tibble(nms=names(lst), data=lst) 
ttt <- tt %>% 
  mutate(data = map(data, ~f(.x))) %>%
  unnest(data)

lst <- import_list(paths, setclass="tbl_df", skip=1, sheet="Affinity scores graph", 
                 rbind=F, rbind_label="Dataset")
names(lst) <- gsub(" PFTs", "", names(lst))
names(lst) <- gsub(" short chron", "", names(lst))

res <- NULL
for (i in 1:length(lst)) {
  xx <- which(colnames(lst[[i]])=="Class")
  x <- data.frame(lst[[i]][, c(ttt$data[i], xx)])
  names(x) <- c("Age_BP2", "Class")
  x$Site_code <- names(lst)[i]
  x$Age_BP2 <- round(x$Age_BP2+0.01, 0)
  res <- bind_rows(res, x)
}  
  
write.xlsx(res, "Woodbridge_Classes.xlsx")

rm(tt, ttt, res, lst, paths, f, i, x, xx)
