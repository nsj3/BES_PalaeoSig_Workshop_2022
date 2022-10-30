##################################################################
#
# Solutions for Wrangling palaeodata with the Tidyverse
#
##################################################################
#
# Base R:
#
##################################################################
#
# Task 1
# poll_pc contains % data.  We can use this to plot a full diagram:

riojaPlot(poll_pc, depth_age, 
       scale.percent=TRUE,
       yvar.name="Age_BP")

# Use whatever method you want, plot the stratigraphic diagram with taxa
# less that 5% max abundance omitted
# Hint: in base R I would use apply or sapply to calculate the max of 
# each column (ie. taxon) and then use a condition to either create a
# logical vector to index the required columns or a numeric vector of 
# column numbers

sp_max <- apply(poll_pc, 2, max)
#or 
sp_max <- sapply(poll_pc, max)

# Then use the max values as a condition to subset columns
tmp1 <- poll_pc[, sp_max > 5]

riojaPlot(tmp1, depth_age, 
       scale.percent=TRUE,
       yvar.name="Age_BP")

############################################################

# Tidyverse

############################################################
#
# Task 2. Using tidy methods and function riojaPlot 
# reproduce the stratigraphic diagram of pollen types omitting taxa 
# with a maximum abundance of less than 5%
# Hint: polllong_sqrt contains the % data (and sqrt data).  Start with 
# this and calculate the max abundance of each taxon and use 
# this to filter out those < 5%.

tmp2 <- polllong_sqrt %>%
  group_by(VarName) %>%
  mutate(max = max(Percent)) %>%
  ungroup() %>%
  filter(max > 5) %>% 
  pivot_wider(id_cols=c(Depth, Age_BP), 
      names_from=VarName, values_from=Percent, values_fill=0)

riojaPlot(tmp2[, -c(1:2)], tmp2[, c(1:2)], 
       scale.percent=TRUE,
       yvar.name="Age_BP")

# A more concise and marginally faster but slightly less 
# intuitive way is to create a function that calculates max value 
# of a vector and compares to a cut value, (5%), then call 
# that function in filter to select rows where max > 5.

fun <- function(x, cut) { max(x) > cut }
tmp3 <- polllong_sqrt %>%
  group_by(VarName) %>%
  filter(fun(Percent, 5)) %>% 
  ungroup() %>%
  pivot_wider(id_cols=c(Depth, Age_BP), 
      names_from=VarName, values_from=Percent, values_fill=0)

riojaPlot(tmp3[, -c(1:2)], tmp3[, c(1:2)], 
       scale.percent=TRUE,
       yvar.name="Age_BP")

# There is a way to accomplish the above without converting to long format
# But it is less intuitive and I don't know a simple way of summing across 
# LCC classes with the data in wide format. We have to use long format 
# form summing LCC classes so it is more efficient to use it for all 
# transformations.
# But if you are curious, here is how to convert the original table of pollen counts 
# to % and remove taxa < 5% max without converting to long format
fun_max <- function(x, cut) { max(x) > cut }
tmp4 <- polldata %>% 
  select(!contains(non_pollen)) %>%   # select all columns except those in non_pollen
  rename("Depth"=`Depth (cm)`, "Age_BP"=`Cal. yr. BP`) %>%   
  rowwise(Depth, Age_BP) %>%  # group by rows, eith Depth and Age_BP as unique row identifiers
  mutate(total=sum(across())) %>%  # use across to calc rowsums
  mutate(across(-total, ~.x/total*100)) %>%  #divide by row totals to calc %
  ungroup() %>%
  select(-total) %>%  # remove column total
  select(where(~fun_max(., 5)) | contains(c("Depth", "Age_BP"))) # select cols with max >5 or Age and Depth

riojaPlot(tmp4[, -c(1:2)], tmp4[, c(1:2)], 
       scale.percent=TRUE,
       yvar.name="Age_BP")


############################################################

# Modelling with the Tidyverse

############################################################
#
# Task 3
# How does richness vary with Affinity (degree of Openness of landscape)?
#
# Join the table of palynolical richness (rich3) to the Affinity scores (allpoll_lcc) and 
# plot richness against affinity.  Is there a relationship for the whole dataset?
# Is the relationship different for different sites.  Is there a regional pattern in the relationship?

# join richness and affinity tables
# remember that -ve values of Affinity indicate more open landscape (less wooded)
# here we multiple affinity by -1 to reverse, so high values indicate more open landscape

rich4 <- rich3 %>%
  left_join(allpoll_lcc, by = c("Site", "Depth", "Age_BP")) %>%
  mutate(Affinity = -Affinity)

# plot relationship for the whole dataset

ggplot(rich4, aes(Affinity, richness)) +
  geom_point() +
  geom_smooth(method="lm")

# Remember that -ve values of Affinity indicate more open landscape (less wooded)
# clear relationship with increasing richness in more open landscapes
# Is it the same for all sites?

ggplot(rich4, aes(Affinity, richness, col=Site)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  theme(legend.position="none")

# Clearly not, but a bit of a mess. Facet by region

ggplot(rich4, aes(Affinity, richness, col=Site)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~Region) +
  theme(legend.position="none")

# Summarise relations with linear regression.  
# Slope give change in richness per unit change in affinity
regs <- rich4 %>%
  group_by(Site, Region, Latitude) %>%
  summarise(broom::tidy(lm(richness ~ Affinity)), .groups="drop") %>%
  filter(term != "(Intercept)") 

regs %>%
    mutate(Region=fct_reorder(Region, Latitude, mean)) %>%
    ggplot(aes(Region, estimate)) +
      geom_boxplot() +
      coord_flip()

# Sites in the S show greater increase in richness for the same degree of "openness" of lansdscape
