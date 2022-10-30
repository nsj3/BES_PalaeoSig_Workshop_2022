library(tidyverse)
library(readxl)
library(broom)

theme_set(theme_bw(base_size=12))

UKLakes <- read_excel("UKLakes.xlsx")

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() 
  
ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm", se=FALSE) 

UKLakes <- UKLakes %>% mutate(DepthClass=
                                cut(MaxDepth, breaks=c(0, 2, 5, 10, 100), 
                                    include.lowest=TRUE), 
                              DepthClass = factor(DepthClass))

ggplot(UKLakes, aes(TP, Chla)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~DepthClass)

UKLakes <- UKLakes %>% mutate(TP_log=log10(TP), Chla_log=log10(Chla))

library(magrittr)
UKLakes %$% cor(Chla_log, TP_log)

UKLakes %>% 
  group_by(DepthClass) %>%
  summarise(r=cor(Chla_log, TP_log))

UKLakes %>% 
  nest(data=-DepthClass) %>%
  mutate(cor = map(data, ~cor.test(.x$TP, .x$Chla))) 

UKLakes %>% 
  nest(data=-DepthClass) %>%
  mutate(cor = map(data, ~tidy(cor.test(.x$TP_log, .x$Chla_log)))) %>%
  unnest(cor)

UKLakes %>% 
  nest(data=-Dataset) %>%
  mutate(reg = map(data, ~ glance(lm(Chla_log ~ TP_log, data=.x)))) %>%
  unnest(reg)

UKLakes %>% 
  nest(data=-Dataset) %>%
  mutate(reg = map(data, ~ tidy(lm(Chla_log ~ TP_log, data=.x)))) %>%
  unnest(reg)

UKLakes %>% 
  nest(data=-Dataset) %>%
  mutate(reg = map(data, ~ augment(lm(Chla_log ~ TP_log, data=.x)))) %>%
  unnest(reg, names_repair="unique") 

