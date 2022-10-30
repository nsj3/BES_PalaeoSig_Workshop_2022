df <- read_excel("extra/PPTX_figures.xlsx", sheet=6)

df %>% mutate(MaxDepth=max(Depth1, Depth2, Depth3))

df %>% 
  rowwise() %>%
  mutate(MaxDepth=max(Depth1, Depth2, Depth3))

df %>% mutate(MaxDepth = pmax(Depth1, Depth2, Depth3))

del <- colnames(red) %in% non_pollen
red <- red[, !del]
sel <- colnames(red) %in% c("Depth (cm)", "Cal. yr. BP")
depth_age <- red[, sel]
spec <- red[, !sel]
colnames(depth_age) <- c("Depth", "Age_BP")
spec_pc <- spec / rowSums(spec) * 100
spec_sqrt <- sqrt(spec_pc)
sel <- match(colnames(spec_sqrt), LCC_taxon_list$VarName)
LCC_code <- LCC_taxon_list$LCC_ID[sel]
spec_LCC <- t (rowsum(t(spec_sqrt), group=LCC_code))
spec_LCC <- data.frame(spec_LCC, check.names=FALSE)
spec_LCC <- spec_LCC / rowSums(spec_LCC) * 100
spec_LCC$A <- rowSums(spec_LCC[, c("1", "2", "3")])
spec_LCC$C <- rowSums(spec_LCC[, c("5", "6", "7")])
spec_LCC$Affinity <- spec_LCC$A - spec_LCC$C
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
LCC2 <- cbind(depth_age, spec_LCC)
LCC2 <- merge(LCC2, LCC2_assem_classes, by="LCC2_ID", all.x=TRUE, sort=FALSE)
LCC2 <- LCC2[order(LCC2$Depth), ]
LCC_names <- LCC_taxon_classes$LCC_name[as.integer(AC_nms)]


red_long <- red %>% 
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  pivot_longer(cols=-c("Depth", "Age_BP"),                   
               names_to="VarName", values_to="Count") %>%    
  filter(Count > 0 & !(VarName %in% non_pollen))                         

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


df <- read_excel("extra/PPTX_Figures.xlsx")

df %>% mutate(across(-(1:2), ~replace_na(.x, 0)))


replace(df, is.na(d), 0)

df[is.na(df)] <- 0

df <- read_excel("extra/PPTX_Figures.xlsx", sheet=3)
df %>% tidyr::separate(Site_Depth, into=c("Site_code", "Depth" )) 


df <- read_excel("extra/UKLakes.xlsx")

df %>% 
  rowwise %>% 
  mutate(Sum_cat = sum(Na+K+Mg+Ca, na.rm=TRUE))



iris <- tibble(iris)

iris

n_iris <- iris %>% nest(data=-Species)
n_iris

n_iris$data %>% map(summary) 

n_iris %>% mutate(summary = map(data, summary))

n_iris %>% mutate(model = map(data, 
              ~lm(Sepal.Width~Sepal.Length, data=.)))


mods <- 



library(tidyr)
library(openxlsx)

tmp <- read_excel("PPTX_figures.xlsx", sheet=1)

tmp2 <- tmp %>% pivot_longer(cols=-c(`Depth (cm)`, `Cal. yr. BP`), 
                             names_to="Taxon", values_to="Count")
write.xlsx(tmp2, "tmp.xlsx")

tmp3 <- tmp2 %>% pivot_wider(id_cols=c(`Depth (cm)`, `Cal. yr. BP`),
                    names_from=Taxon, values_from=Count)

df <- tmp2[1:9, ]

df %>% 
  group_by(`Depth (cm)`) %>%
  summarise(Total=sum(Count))

df %>% 
  group_by(`Depth (cm)`) %>%
  mutate(Percent=Count/sum(Count)*100)

df %>% 
  mutate(SQRT_Count=sqrt(Count))

x <- 1:5
for (i in x) {
  print(i)
}

my_fun <- function(a, b) {
  c <- a + b
  return(c)
}
my_fun(2, 3)


with(poll_comb, glance(cor.test(Affinity, richness)))
