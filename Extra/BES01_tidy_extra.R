SiteInfo <- SiteInfo %>% mutate(Region2=fct_collapse(Region,
                                        Scotland=c("East Scotland", 
                                                   "North Scotland", 
                                                   "South Scotland"), 
                                        England=c("Southeast England",
                                                  "Central England",
                                                  "Northeast England",
                                                  "Southwest England")))







red_wide <- red_long_trans %>% 
  group_by(VarName) %>%
  mutate (max=max(Percent)) %>%
  filter(max > 2) %>%
  pivot_wider(id_cols=c(Depth, Age_BP), 
  names_from=VarName, values_from=Percent, values_fill=0)

# Then plot:
strat.plot(red_wide[,-c(1:2)], yvar=red_wide$Age_BP, scale.percent=TRUE, 
           y.rev=TRUE, yTop=0.7, cex.title=1.2, cex.xlabel=0.7, 
           plot.poly=TRUE, col.poly="darkgreen", ylabel="Age (yr BP)", 
           srt.xlabel=45, xLeft=0.1, xRight=0.95, ylabPos=2.8) 




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

tmp <- tibble(sample_id=as.character(1:nrow(red_wide2)), red_wide2) %>% 
  rename(age="Age_BP", depth="Depth")
roc <- fc_estimate_RoC(tmp[, -c(2:3)], tmp[, 1:3], Working_Units="bins", bin_size=20)
roc

roc2 <- fc_detect_peak_points(roc, method="GAM_deriv")
roc2
fc_plot_RoC_sequence(roc2)

p4 <- ggplot(roc, aes(x=Age, y=ROC)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) +
  labs(x="Years BP")

p1 + p2 + p3 + p4 + plot_layout(ncol=1)

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
  

# save as excel file
