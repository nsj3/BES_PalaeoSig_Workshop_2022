##################################################################
#
# Solutions for PSIG02_tidy.R
#
##################################################################

# Task1: Produce a plot of land cover classes (as above, line 144) but 
# showing England and Scotland separately.
# Hint: Copy the aggregation code on line 78 and add Region2 to the grouping, 
# then facet the plot by region.

allpoll_agg2 <- allpoll_LCC2 %>% 
  left_join(Site_info, by="Site_code") %>%
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Region2, Age2, LCC2_ID) %>%
  summarise(N=n(), .groups="drop_last") %>%
  mutate(PC=N/sum(N)*100) %>%
  pivot_wider(id_cols=c(Region2, Age2), 
              names_from=LCC2_ID, values_from=PC, values_fill=0) %>%
  pivot_longer(cols=-c(Region2, Age2), names_to="LCC2_ID", values_to="PC") %>%
  mutate(LCC2_ID=as.integer(LCC2_ID)) %>%
  left_join(LCC2_assem_classes, by="LCC2_ID") 

ggplot(allpoll_agg2, aes(x=PC, y=Age2)) +
  geom_col(orientation="y", width=200, 
           position=position_nudge(y=-100), fill="black") +
  facet_grid(Region2 ~ LCC2_name, scales="free_x", space="free_x") +
  scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 10000, by=1000)) +
  labs(x="%", y="Years BP") 


# Task2: Create a figure of the number of levels per timeslice for 0-9000 years 

tmp <- allpoll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2) %>%
  summarise(N=n(), .groups="drop_last") 

ggplot(tmp, aes(Age2, N)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=1000)) 


# Task3: Create a plot showing trends in Affinity scores for 0-9000 years 
# for for each site grouped by Region.  Note here we just pipe the transformed data 
# directly to ggplot

allpoll_LCC2 %>%
  left_join(Site_info, by="Site_code") %>%
  filter(Age_BP <= 9000) %>%
  ggplot(aes(Age_BP, Affinity, col=Site_code)) +
  geom_line() +
  facet_wrap(~Region) +
  theme(legend.position="none")
  

# Task4: Let's look at the elm decline.  Create a plot of % elm (Ulmus) against 
# time for the last 9k yr, grouped by region.
# Hint: Modify the code on line 43 so it only applies the function fun_long_tidy.
# This will create a long table of pollen counts and percentages.
# Filter this by Age_BP < 9000 and by VarName to extract Ulmus, 
# join with Site_info to add Region column and plot.

allpoll_Ulmus <- allpoll_nested %>% 
  mutate(data, data=map(data, ~fun_long_tidy(.x, non_pollen=non_pollen))) %>%
  unnest(data) %>%
  filter(Age_BP<9000 & VarName=="Ulmus") %>%
  left_join(Site_info, by="Site_code")

ggplot(allpoll_Ulmus, aes(Age_BP, Percent, col=Site_code))  +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 10000, by=2000)) +
  facet_wrap(~Region, scales="free_y") +
  theme(legend.position="none") 

# Task5: Calculate the mean and sd affinity score for each time slice for 
# 0-9000 years and create plot showing overall mean affinity vs. Age with 
# error-bars showing standard deviation.

allpoll_aff <- allpoll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2) %>%
  summarise(Aff_mean=mean(Affinity), Aff_sd=sd(Affinity)) 

# a more concise way of writing the above summarise operation 
# when summarising many columns is:
allpoll_aff <- allpoll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, 
                        width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Age2) %>%
  summarise(across(.cols=Affinity, 
                   .fns=list(mean=mean, sd=sd), 
                   .names="Aff_{.fn}")) 

ggplot(allpoll_aff, aes(Age2, y=Aff_mean, ymin=Aff_mean-Aff_sd, 
                        ymax=Aff_mean+Aff_sd)) +
  geom_line() +
  geom_errorbar()

# Task6: Repeat the above plot, split by Region.

allpoll_aff <- allpoll_LCC2 %>% 
  filter(Age_BP <= 9000) %>%
  left_join(Site_info, by="Site_code") %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Region, Age2) %>%
  summarise(Aff_mean=mean(Affinity), Aff_sd=sd(Affinity), .groups="drop") 

ggplot(allpoll_aff, aes(Age2, y=Aff_mean, ymin=Aff_mean-Aff_sd, 
                     ymax=Aff_mean+Aff_sd)) +
  geom_line() +
  geom_errorbar() +
  facet_wrap(~Region)

