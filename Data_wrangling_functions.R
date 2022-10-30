# Change ggplot default theme
theme_set(theme_bw(base_size=12))

##################################################################
#
# Functions for calculating LCC classes using Base R and Tidyverse
#
##################################################################

fun_LCC_base <- function(x, LCC_lookup, non_pollen) {
   del <- colnames(x) %in% non_pollen
   x <- x[, !del]
# split data into depth/age & species parts
   depth_age <- subset(x, select=c(`Depth (cm)`, `Cal. yr. BP`))
   poll_count <- subset(x, select=-c(`Depth (cm)`, `Cal. yr. BP`))
#rename depth / age columns
   colnames(depth_age) <- c("Depth", "Age_BP")
# transform counts to percentages
   poll_pc <- poll_count / rowSums(poll_count) * 100
# transform to sqrt
   poll_sqrt <- sqrt(poll_pc)
# aggregate columns into LCC groups
# First create a vector of LCC groups corresponding to taxon names
   sel <- match(colnames(poll_sqrt), lcc_lookup$VarName)
# Then extract LCC class for each taxon
   taxa_lcc <- lcc_lookup[sel, ]
# calculate column sums
   poll_lcc <- t (rowsum(t(poll_sqrt), group=taxa_lcc$LCC_name))
   poll_lcc <- data.frame(poll_lcc, check.names=FALSE)
# we have sums of sqrt data, so normalise data
   poll_norm <- poll_lcc / rowSums(poll_lcc) * 100
   dominant_class <- apply(poll_norm, 1, which.max)
   poll_norm$lcc_class <- colnames(poll_norm)[dominant_class]
   sum_arboreal <- rowSums(poll_norm[, c("Coniferous woodland", 
                                     "Deciduous woodland", 
                                     "Wet/fen woodland")])  # arboreal classes
   sum_open <- rowSums(poll_norm[, c("Heath", 
                                     "Pasture/meadow", 
                                     "Arable indicators")])  # open classes
   poll_norm$Affinity <- sum_arboreal - sum_open
   result <- cbind(depth_age, poll_norm)
}

fun_long_tidy <- function(x, non_pollen) {
  x_long <- x %>% 
    rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
    pivot_longer(cols=-c("Depth", "Age_BP"),                   
                 names_to="VarName", values_to="Count") %>%    
    filter(!(VarName %in% non_pollen)) 
  x_long
}


