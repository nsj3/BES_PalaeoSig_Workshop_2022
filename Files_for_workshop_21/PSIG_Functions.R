# Change ggplot default theme
theme_set(theme_bw(base_size=12))

# import supporting lookup tables
LCC_taxon_list <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_List")
LCC_taxon_classes <- read_excel("LCC_Info.xlsx", sheet="LCC_Taxon_Classes")
LCC2_assem_classes <- read_excel("LCC_Info.xlsx", sheet="LCC2_Assemblage_Classes")
Site_info <- read_excel("LCC_Info.xlsx", sheet="Site_Info")

# create character vector of non-pollen variables to remove
non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")

##################################################################
#
# Functions for calculating LCC classes using Base R and Tidyverse
#
##################################################################

fun_LCC2_base <- function(x, LCC_taxon_list, LCC_taxon_classes) {
   non_pollen <- c("Sample",
                "Radiocarbon years B.P.",
                "EPD default [yrs.BP.]",
                "EPD [yrs.BP.]",
                "Fossilva [yrs.BP.]", 
                "Sum")
   del <- colnames(x) %in% non_pollen
   x <- x[, !del]
   # split data into depth/age & species parts
   depth_age <- subset(x, select=c(`Depth (cm)`, `Cal. yr. BP`))
   spec <- subset(x, select=-c(`Depth (cm)`, `Cal. yr. BP`))
   #rename depth / age columns
   colnames(depth_age) <- c("Depth", "Age_BP")
   # transform counts to percentages
   spec_pc <- spec / rowSums(spec) * 100
   # transform to sqrt
   spec_sqrt <- sqrt(spec_pc)

   # aggregate columns into LCC groups
   # create vector of LCC groups corresponding to columns names
   # first, look up index of colnames in LCC taxon list
   sel <- match(colnames(spec_sqrt), LCC_taxon_list$VarName)

   # extract LCC class for each taxon
   LCC_code <- LCC_taxon_list$LCC_ID[sel]

   # calculate column sums
   # no easy way to calc grouped columns sums directly
   # so we transpose data, calculate grouped rowsums then transpose 
   # result back
   spec_LCC <- t (rowsum(t(spec_sqrt), group=LCC_code))
   spec_LCC <- data.frame(spec_LCC, check.names=FALSE)
   # normalise data
   spec_LCC <- spec_LCC / rowSums(spec_LCC) * 100

   # Calculate the sum of tree (A = LCC classes 1-3) 
   # and open (C = LCC classes) and  Affinity score (A-C)
   spec_LCC$A <- rowSums(spec_LCC[, c("1", "2", "3")])
   spec_LCC$C <- rowSums(spec_LCC[, c("5", "6", "7")])
   spec_LCC$Affinity <- spec_LCC$A - spec_LCC$C
   # Create land cover class for each level in the core (LCC2)
   # if Affinity > 20, LCC2 = max class max from 1-3 (woodland types)
   # if Affinity < 20, LCC2 = max class  class from 5-7 (open types)
   # if Affinity -20 to +20, if LCC = woodland type, then LCC2 = semi-open arboreal (8)
   # if Affinity -20 to +20, if LCC = open type, then LCC2 = semi-open type (9-11)

   # This is horrible to code in base R.  
   # Could do in a loop using switch on each row, or vectorise it, 
   # or use nested ifelse() but this is really horrible... 
   # Help is on hand with a vectorised nested if function (nif) in package kit:
   A_nms <- c("1", "2", "3")
   C_nms <- c("5", "6", "7")
   AC_nms <- c(A_nms, C_nms)

   LCC2_ID <- kit::nif(
   spec_LCC$Affinity > 20, as.integer(A_nms[apply(spec_LCC[, A_nms], 1, which.max)]),
   spec_LCC$Affinity < -20, as.integer(C_nms[apply(spec_LCC[, C_nms], 1, which.max)]),
   default = ifelse(as.integer(AC_nms[apply(spec_LCC[, AC_nms], 1, which.max)]) < 4, 
                   8L,
                   as.integer(AC_nms[apply(spec_LCC[, AC_nms], 1, which.max)]) + 4L)
   )

   spec_LCC$LCC2_ID <- LCC2_ID
   # merge with depths
   LCC2 <- cbind(depth_age, spec_LCC)
   
   # merge with LCC2 lookup table to add LCC2 names for plotting
   LCC2 <- merge(LCC2, LCC2_assem_classes, by="LCC2_ID", sort=FALSE)
   # merge re-sorts data by LCC2_ID!!! Why?
   # sort data by Depth
   LCC2 <- LCC2[order(LCC2$Depth), ]
   LCC2
}

fun_long_tidy <- function(x, non_pollen) {
  x_long <- x %>% 
    rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
    pivot_longer(cols=-c("Depth", "Age_BP"),                   
                 names_to="VarName", values_to="Count") %>%    
    filter(!(VarName %in% non_pollen))                         

  x_long_trans <- x_long %>% 
     group_by(Depth) %>%
     mutate(Percent = Count / sum(Count) * 100, 
            SQRT_PC = sqrt(Percent)) %>%
     ungroup() 
  x_long_trans
}


fun_LCC2_tidy <- function(x, LCC_taxon_list, LCC_taxon_classes, 
                      LCC2_assem_classes) {
  x_long_LCC <- x %>% 
    left_join(LCC_taxon_list, by="VarName") %>%
    group_by(Depth, Age_BP, LCC_ID) %>%
    summarise(Percent=sum(Percent), 
              SQRT_PC=sum(SQRT_PC), 
              .groups="drop_last") %>%
    mutate(Norm_SQRT_PC=SQRT_PC/sum(SQRT_PC)*100) %>%
    ungroup() %>%
    left_join(LCC_taxon_classes, by="LCC_ID") 

  x_long_LCC2 <- x_long_LCC %>%
    group_by(Depth, Age_BP) %>%
    mutate(A=sum(Norm_SQRT_PC[LCC_group=="A"]), 
           C=sum(Norm_SQRT_PC[LCC_group=="C"]), 
           Affinity=A-C) %>%
    slice_max(Norm_SQRT_PC, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(LCC2_ID=case_when(
                 Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 1:3 ~ 8,
                 Affinity >= -20 & Affinity <= 20 & LCC_ID %in% 5:7 ~ LCC_ID + 4,
                 TRUE ~ LCC_ID
                 )
           ) %>%
    left_join(LCC2_assem_classes, by="LCC2_ID")
    x_long_LCC2
}



