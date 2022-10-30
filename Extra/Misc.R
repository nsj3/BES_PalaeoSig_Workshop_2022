ggplot(d7, aes(ymin=Age2, ymax=Age2+200, xmin=0, xmax=PC)) +
  geom_rect(fill="black") +
  facet_grid(. ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP")



nms <- d3 %>% distinct(P_Var) %>% 
  data.frame()

nms

#w <- c(d=0.1, i=0.1, s=1, t=1)
#nms3 <- stringdist_left_join(nms, p_vars2, by=c("P_Var"="VarName"), 
#                             weight=w, max_dist=5, distance_col="Dist")

#nms4 <- nms3 %>% group_by(P_Var) %>% 
#  slice_min(Dist, n=1, with_ties=FALSE) %>%
#  select(P_Var, VarName, PFT, GroupId) %>%
#  ungroup()
# nms4 %>% distinct(GroupId)

#keep_groups <- c("TRSH", "HERB", "DWAR")

#keep_vars <- nms4 %>%
#  filter(GroupId %in% keep_groups)

#d4 <- d3 %>% filter(P_Var %in% keep_vars$P_Var) %>%
#  group_by(Dataset, Depth) %>%
#  mutate(Percent=Count / sum(Count) * 100) %>%
#  ungroup()
