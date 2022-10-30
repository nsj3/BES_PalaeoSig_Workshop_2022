##################################################################
#
# Solutions for PSIG02_tidy.R
#
##################################################################

UKLakes <- read_excel("UKLakes.xlsx")

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 

# Task 1: Reproduce the above scatter plot grouped by depth class
# classes = 0-2, 2-5, 5-10 and >10m
# Hint: you will need a create a new variable DepthClass using mutate
# and use cut to classify each observation into a depth class
# Then use facet_wrap on your DepthClass

UKLakes <- UKLakes %>%
  mutate(DepthClass = cut(MaxDepth, breaks=c(0, 2, 5, 10, 100)))

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~DepthClass)

# Add smooth

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~DepthClass)

# Task 2: Use group_by and summarise to calculate the correlation 
# between log10-transformed Chla and TP for each depth class
# Hint: you can either use mutate to create new log10 variables or
# include log10 transformation in the call to cor

# I prefer to calculate new log10 variables as we will use the data 
# in a number of different tasks
UKLakes <- UKLakes %>%
  mutate(TP_log=log10(TP), Chla_log=log10(Chla))

# or if you have many variables to transform:
UKLakes <- UKLakes %>%
  mutate(across(c(TP, Chla), .fns=list(log10=log10), .names="{.col}_log"))

UKLakes %>% 
  group_by(DepthClass) %>%
  summarise(r=cor(TP_log, Chla_log))

# Task 3: Use nest and mutate / map to calculate correlation & significance 
# of each depth class using cor.test

UKLakes %>% 
  group_by(DepthClass) %>%
  nest() %>%
  mutate(data, r=map(data, ~tidy(cor.test(.$TP_log, .$Chla_log)))) %>%
  unnest(r)

# Task 4: Use nest and mutate / map to fit a linear regression model 
# of log10 transformed Chla on TP for each Depth class

UKLakes %>% 
  group_by(DepthClass) %>%
  nest() %>%
  mutate(data, reg=map(data, ~tidy(lm(Chla_log ~ TP_log, data=.)))) %>% 
  unnest(reg)

# or ...
UKLakes %>% 
  group_by(DepthClass) %>%
  nest() %>%
  mutate(data, reg=map(data, ~glance(lm(Chla_log ~ TP_log, data=.)))) %>% 
  unnest(reg)

# Task 5: Examine rates of palynological change through the Holocene
# We can use the function fc_estimate_RoC in package RRatepol to estimate 
# rates of change.  The code for this is at the end of this script.  It takes 
# 30second on a fast machine, substantially longer on a slow PC. 
# To save time load the results from an excel file:

poll_roc <- read_excel("Woodbridge_ROC.xlsx")

# Task 5A: Plot the ROC data against age, split by region
# Hint age is in column Age, rate of change in column ROC
# There are some extreme values - you can use ylim to limit the range 
# of the y-axis, but ylim removes data, rather than truncates or clips the plot
# This is fine with point data but for lines it is better to use 
# coord_cartesian(ylim=c(0, 2)) to clip so no data is removed. 

ggplot(poll_roc, aes(Age, y=ROC, col=Site_code)) +
  geom_line() +
  facet_wrap(~Region) +
  coord_cartesian(ylim=c(0, 2)) + 
  theme(legend.position="none") 

# Task 5C: Aggregate the data by Region and 200-year Age slice and 
# calculate the mean and sd of each slice

poll_roc2 <- poll_roc %>%
  mutate(Age2=cut_width(Age, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Region, Age2) %>%
  summarise(ROC_mean = mean(ROC), ROC_sd=sd(ROC), .groups="drop")

# Task 5B Plot the mean + sd values
ggplot(poll_roc2, aes(Age2, y=ROC_mean, 
                      ymin=ROC_mean-ROC_sd, 
                      ymax=ROC_mean+ROC_sd)) +
  geom_line() +
  geom_errorbar() +
  facet_wrap(~Region) +
  coord_cartesian(ylim=c(0, 2))

