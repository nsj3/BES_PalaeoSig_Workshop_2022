wc <- read_excel("Woodbridge_Classes.xlsx") 

%>%
  filter(Site_code %in% unique(d6$Site_code))

wc3 <- wc %>% 
  filter(Age_BP2 <= 9200) %>%
  mutate(Age2=cut_width(Age_BP2, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Site_code, Age2, Class) %>%
  count() %>%
  filter(Class==6 & Age2 > 8000)
wc3

wc2 <- wc %>% 
  filter(Age_BP2 <= 9000) %>%
  mutate(Age2=cut_width(Age_BP2, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Site_code, Age2, Class) %>%
  count() %>%
  slice_max(n, n=1, with_ties=FALSE) %>%
  group_by(Age2, Class) %>%
  summarise(PC=sum(n), .groups="drop_last") %>%
  mutate(PC=PC/sum(PC)*100) %>%
  pivot_wider(names_from=Class, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-Age2, names_to="Class", values_to="PC") %>%
  mutate(Class=as.integer(Class)) %>%
  left_join(LCC_assemblage, by=c("Class"="LCC_ID")) %>%
  ungroup()

wc2$LCC_name <- factor(wc2$LCC_name, levels=LCC_order)
ggplot(wc2, aes(x=Age2, xmin=Age2, xmax=Age2+200, ymin=0, ymax=PC)) +
  geom_rect(data=wc2, fill="black") +
  facet_grid(. ~ LCC_name, scales="free_x", space="free_x") +
  scale_y_continuous(breaks=seq(0, 100, by=20)) +
  scale_x_reverse(breaks=seq(0, 10000, by=1000)) +
  coord_flip() + 
  labs(y="%")
