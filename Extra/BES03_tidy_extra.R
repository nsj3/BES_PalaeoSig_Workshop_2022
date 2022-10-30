stop()

tmp <- rename(red_wide, depth="Depth", age="Age_BP")
tmp$sample_id <- as.character(1:nrow(tmp))
tmp <- tmp %>% relocate(sample_id, .before=depth)
library(RRatepol)
system.time(
  roc2 <- fc_estimate_RoC(tmp[, -c(2:3)], tmp[, 1:3], 
                  Working_Units="bins", treads=F, 
                  bin_size=200, 
                  rand=100)

    )



SiteInfo <- SiteInfo %>% mutate(Region2=fct_collapse(Region,
                                        Scotland=c("East Scotland", 
                                                   "North Scotland", 
                                                   "South Scotland"), 
                                        England=c("Southeast England",
                                                  "Central England",
                                                  "Northeast England",
                                                  "Southwest England")),
                                col=palette()[as.integer(factor(Region))])

library(leaflet)

leaflet(data=SiteInfo) %>%
  addTiles() %>%
  addLabelOnlyMarkers(lng=~Longitude, lat=~Latitude, 
                      label=~Site_code, 
                      labelOptions=labelOptions(noHide=TRUE, textsize="10px",
                                                direction="top", textOnly=TRUE)) %>%
  addProviderTiles(providers$CartoDB.Positron)

leaflet() %>%
  addTiles() %>%
  addLabelOnlyMarkers(data=SiteInfo, lng=~Longitude, lat=~Latitude, 
                      label=purrr::map(glue::glue("<span style='color:{SiteInfo$col}'>{as.character(SiteInfo$Site_code)}<span>"), htmltools::HTML), labelOptions=labelOptions(noHide=TRUE, textsize="10px",
                                                direction="top", textOnly=TRUE,
                                                style=list('color'=~col))) %>%
  addProviderTiles(providers$CartoDB.Positron)


fun_PrCurve <- function(x) {
   x <- x %>% filter(Age_BP <= 9000)
   prc <- principal_curve(sqrt(as.matrix(x[, -c(1:2)])), 
                        trace=FALSE, penalty=1.4)
   prc2 <- prc$lambda
   if (cor(prc2, x$Age_BP) < 0)
     prc2 <- -prc2
   prc2 <- scale(prc2, T, T)
   tibble(x[, 1:2], PrCurve=prc2)
}

poll_PrCurve <- tmp %>%
  nest(data=-Site_code) %>%
  mutate(data, data=map(data, ~pivot_wider(.x, id_cols=c(Depth, Age_BP), 
                                           names_from=VarName, 
                                           values_from="Percent", 
                                           values_fill=0))) %>%
  mutate(data, data=map(data, ~fun_PrCurve(.x)))

poll_PrCurve <- poll_PrCurve %>% unnest(data) %>%
  left_join(SiteInfo2, by="Site_code")

ggplot(poll_PrCurve, aes(x=Age_BP, y=PrCurve, group=Site_code)) +
  geom_line() +
  facet_wrap(~ Region)

 d2A_PrCurve <- d2_PrCurve %>% 
  mutate(Age2=cut_width(Age_BP, boundary=0, width=400, labels=FALSE) * 400 - 400) %>%
  group_by(Site_code, Age2) %>%
  summarise(PrCurve=mean(PrCurve), .groups="drop_last") %>%
  left_join(SiteInfo2, by="Site_code") %>%
  mutate(N=n()) %>%
  filter(N > 11) 

d3_PrCurve <- d2A_PrCurve  %>%
#  filter(Region2 == "England") %>%
  pivot_wider(id_cols=Age2, names_from=Site_code, values_from=PrCurve) %>%
  arrange(Age2)

cc <- 1 - abs(cor(d3_PrCurve[, -c(1)], use="pairwise.complete.obs"))
cc2 <- as.dist(cc)
hc <- hclust(cc2)
plot(hc, hang=-1)
gr1 <- cutree(hc, k=7)
gr <- data.frame(Site_code=names(gr1), Group=gr1)

dd <- d2A_PrCurve %>% left_join(gr)
ggplot(dd, aes(x=Age2, y=PrCurve, group=Site_code)) +
  geom_line() +
  facet_wrap(~ Group)

d4A_wide <- d4  %>%
  mutate(Age2=cut_width(Age_BP, boundary=0, width=400, labels=FALSE) * 400 - 400) %>%
  group_by(Site_code, Age2) %>%
  summarise(Ulmus=mean(Ulmus), .groups="drop_last") 

d4_wide <- d4A_wide %>%
  pivot_wider(id_cols=Age2, names_from=Site_code, values_from=Ulmus) %>%
  arrange(Age2)

cc <- 1 - abs(cor(d4_wide[, -c(1)], use="pairwise.complete.obs"))
cc[is.na(cc)] <- 0
cc2 <- as.dist(cc)
hc <- hclust(cc2)
plot(hc, hang=-1)
gr1 <- cutree(hc, k=7)
gr <- data.frame(Site_code=names(gr1), Group=gr1)

dd <- d4A_wide %>% left_join(gr)
ggplot(dd, aes(x=Age2, y=Ulmus, group=Site_code)) +
  geom_line() +
  facet_wrap(~ Group)

fun2 <- function(x) {
  if (!("Ulmus" %in% colnames(x)))
    x$Ulmus <- 0
  select(Depth, Age_BP, Ulmus)
}

d2 <- d %>% 
  mutate(data, data=map(data, ~pivot_wider(.x, names_from=VarName, 
                                      values_from="Percent", 
                                      values_fill=0))) %>%
  unnest(data)

d4 <- d2 %>% filter(Age_BP <= 9000) %>%
  left_join(SiteInfo2, by="Site_code") 

ggplot(d4, aes(x=Age_BP, y=Ulmus, col=Site_code)) +
  geom_line() +
  facet_wrap(~ Region, scales="free_y") +
  theme(legend.position="none")

fun3 <- function(x, title="") {
  strat.plot(x[, -c(1:2)], x[, 2], scale.percent=TRUE, 
             y.rev=TRUE, title=title, yTop=0.7, cex.title=1.2,
             cex.xlabel=1)
}

fun5 <- function(data) {
  data %>% group_by(VarName) %>%
  mutate(max = max(Percent)) %>%
  filter(max > 2) %>% 
  select(-max) %>%
  pivot_wider(names_from=VarName,
              values_from="Percent", 
              values_fill=0)
}

pdf("Pollen_plots.pdf", width=11, height=9, paper="a4r")
d %>%
  mutate(data, data=map(data, ~fun5(.x))) %>%
  walk2(.x=.$data, .y=.$Site_code, .f=~fun3(.x, .y))
dev.off()


fun_ROC <- function(x) {
   x <- x %>% filter(Age_BP <= 9000)
   x2 <- tibble(sample_id=as.character(1:nrow(x)), x) %>%
     rename(age="Age_BP", depth="Depth")
   roc <- fc_estimate_RoC(x2[, -c(2:3)], x2[, 1:3], 
                          Working_Units="bins", bin_size=20)
   if (nrow(x) > 25)
      fc_detect_peak_points(roc, method="GAM_deriv")
   else
      fc_detect_peak_points(roc, method="trend_non_linear")
}

d_ROC <- d %>%
  mutate(data, data=map(data, ~pivot_wider(.x, names_from=VarName, 
                                      values_from="Percent", 
                                      values_fill=0))) %>%
  mutate(data=map(data, ~fun_ROC(.x))) %>%
  unnest(data)

d_ROC <- d_ROC %>% 
  left_join(SiteInfo2, by="Site_code")

ggplot(d_ROC, aes(Age, ROC, group=Site_code, col=Peak)) +
  geom_line() +
#  geom_point() +
  facet_wrap(Region~., scales="free_y") +
  theme(legend.position="none")

ggplot(d_ROC, aes(x=Age, y=ROC, group=Site_code)) +
  geom_line() +
  facet_wrap(~ Region, scales="free_y")

d2A_ROC <- d_ROC %>% 
  mutate(Age2=cut_width(Age, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Site_code, Age2) %>%
  summarise(nROC=sum(Peak), .groups="keep") %>%
  summarise(nROC2=sum(nROC>0), .groups="drop") %>%
  left_join(SiteInfo2, by="Site_code") %>%
  group_by(Region, Age2) %>%
  summarise(nROC2=sum(nROC2) / n(), .groups="drop")

ggplot(d2A_ROC, aes(Age2, nROC2)) +
  geom_col(width=200) +
  facet_wrap(~Region)

d2A_ROC <- d_ROC %>% 
  mutate(Age2=cut_width(Age, boundary=0, width=200, labels=FALSE) * 200 - 200) %>%
  group_by(Site_code, Age2) %>%
  summarise(ROC=mean(ROC), .groups="drop_last") %>%
  left_join(SiteInfo2, by="Site_code") %>%
  mutate(N=n()) %>%
  filter(N > 20) 

tmp2 <- d2A_ROC  %>%
#  filter(Region2 == "England") %>%
  pivot_wider(id_cols=Age2, names_from=Site_code, values_from=ROC) %>%
  arrange(Age2)

cc <- 1 - abs(cor(tmp2[, -c(1)], use="pairwise.complete.obs"))
cc2 <- as.dist(cc)
#cc2 <- replace_na(cc2, 10)
hc <- hclust(cc2)
plot(hc, hang=-1)
gr1 <- cutree(hc, k=8)
gr <- data.frame(Site_code=names(gr1), Group=gr1)

dd <- d2A_ROC %>% left_join(gr)
ggplot(dd, aes(x=Age2, y=ROC, group=Site_code)) +
  geom_line() +
  facet_wrap(~ Group, scales="free_y")


system.time(source("BES02_base.R"))
system.time(source("BES02_tidy.R"))
