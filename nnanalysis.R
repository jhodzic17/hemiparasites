## Paired Plot Analysis
## Jasna Hodzic
## 17.11.20

# run originally with R 4.0.3. in RStudio 1.3.1093

#library(sf)
library(dplyr)
#library(leaflet)
#library(htmltools)
#library(nngeo)
library(ggplot2)
library(rcompanion)
library(tidyverse)
library(lme4)
library(lmerTest)

###  - can scroll down to analysis, first ~ 70 lines of code are creating the nn in GIS

# read in shapefile
setwd("~/NPS")
shp = st_read('GIS/nn1.shp') ## shape file with ecoregion information
# shp = shp %>% slice(1:1000)

# find the nn distance (this takes a few minutes)

nn_id = st_nn(shp, shp, k = 2, sparse = T, returnDist = T, parallel = parallel::detectCores()-1)
 nn_tib = tibble(id = as.integer(), nn_id = as.integer(), nn_dist_m = as.numeric())

# Start building lookup table with nn data
for(i in seq(1,length(shp$SITEPLOT))){
  id_tmp = nn_id$nn[[i]][[2]]
  dist_tmp = nn_id$dist[[i]][2]
  
  nn_tib = nn_tib %>%
    bind_rows(c(id = i, nn_id = id_tmp, nn_dist_m = dist_tmp))
}

# Add SITEPLOT for joining purposes
nn_tib = nn_tib %>% bind_cols(SITEPLOT = shp$SITEPLOT)

# Add nn plot name
nn_plot_name = c()
for(i in nn_tib$nn_id){
  nn_plot_name = c(nn_plot_name, nn_tib$SITEPLOT[i])
}

nn_tib = nn_tib %>% 
  bind_cols(nn_plot_name = nn_plot_name)

# visualize chunks of the data and verify
# smaller radius, stroke, and weight makes a faster drawing speed
leaflet() %>%
  addTiles() %>%
  addCircleMarkers(data = shp %>% filter(SITEPLOT %in% c('AGFO 18', 'AGFO 20')), lng = ~LONG_, lat = ~LAT, radius = 10, stroke = 1, weight = 2, popup = ~htmlEscape(SITEPLOT)) 
  
# if everything looks good, join the nn data to the shapefile and export
shp = shp %>%
  inner_join(nn_tib %>% select(SITEPLOT, nn_dist_m, nn_plot_name), by = c('SITEPLOT' = 'SITEPLOT')) %>%
  relocate(any_of(c('nn_dist_m', 'nn_plot_name')), .after = SITEPLOT)

# Write out CSV with paired plots for reading 
write.csv(shp, "PARKHEMIFINAL/gisfinalpairs.csv", row.names = F)

# you can also write this out as a shapefile with the sf library
# st_write(obj = shp, dsn = 'pts_4326_nn.shp')

############### Analysis ###################
#
## csv with pairs
setwd("~/NPS")
data.n <- read.csv("PARKHEMIFINAL/gisfinalpairs.csv" , header=TRUE)
data.eco <- read.csv("PARKHEMIFINAL/parkhemi.ecoall.020421.csv", header = TRUE)
area <- read.csv("plots_areasampled.csv", header = TRUE)

pairs <- data.n %>% dplyr::select(SITEPLOT , geometry, Shape_Leng, nn_dist_m, nn_plot_name)

data.nn <- merge(data.eco, pairs)
data.nn <- merge(data.nn, area)
head(data.nn)
##testing to replace "NA" in herbeven with 0



data.nn<- data.nn%>% mutate(From = nn_plot_name,
                      To = SITEPLOT)
                 

data.nn$unique.pair <- with(data.nn, paste(From, To, sep = "_"))


pairs <- data.nn %>% dplyr::select(To, From, nn_dist_m, geometry, unique.pair)

data.JH <- data.nn %>% dplyr::select(! c(To, From, nn_dist_m, geometry, unique.pair))




## create subset of plots that contain hemiparasites
IN.hemi <- pairs %>%
  merge(y = data.JH, by.x = "To", by.y= "SITEPLOT")%>%
  filter(HEMIAB > 0)



##  create subset of plots that do not contain hemiparasites
NEAR.hemi <- pairs %>%
  merge(y = data.JH, by.x = "From", by.y = "SITEPLOT") %>%
  filter(HEMIAB == 0)


NEAR.hemi. <- NEAR.hemi %>% filter(!is.na(Area_sampled))

master.JH <- IN.hemi %>%
  merge(y = NEAR.hemi, by = c("To", "From", "nn_dist_m", "unique.pair"),
        suffixes = c(".IN", ".NEAR"))


# columns that end in ".IN" are data from the plot with hemiparasites;
# columns that end in ".NEAR" are data from the nearby plot without hemiparasites

## 38 taken out due to lat/long (same codes)

master.JH <- master.JH %>% filter(!nn_dist_m==0) ## final is 1244 observations

## average distance


master.JH$park.same <- as.factor(ifelse(master.JH$PARK.IN == master.JH$PARK.NEAR, "Y","N"))
master.JH$L4.same <- as.factor(ifelse(master.JH$L4_KEY.IN == master.JH$L4_KEY.NEAR, "Y","N"))
master.JH$L3.same <- as.factor(ifelse(master.JH$L3_KEY.IN == master.JH$L3_KEY.NEAR, "Y","N"))
master.JH$Area.same <- as.factor(ifelse(master.JH$Area_sampled.IN == master.JH$Area_sampled.NEAR, "Y","N"))

summary(master.JH$park.same)
summary(master.JH$Area.same)
master.JH[master.JH$park.same == "N",] ## consider dropping?
L4<-master.JH[master.JH$L4.same == "N",]
master.JH[master.JH$L3.same == "N",] ## 5
area.n <- master.JH[master.JH$Area.same == "N",]
area.n <- area.n %>% na.omit()
### drop pairs:
 
master.JH <- master.JH %>% filter(L4.same == "Y" & park.same == "Y") ## 1195 pairs

## what about restricting to the same area


master.JH <- master.JH %>% filter(!is.na(NOHEMIEVEN.IN), !is.na(NOHEMIEVEN.NEAR), !is.na(HERBEVEN.IN), !is.na(HERBEVEN.NEAR), 
        !is.na(ALLEVEN.NEAR), !is.na(ALLEVEN.IN)) ## 11

head(master.JH)
## dropping 57 pairs which have NA's in herbeven or NOHEMIEVEN; i.e. all due to richness = 1.0

write.csv(master.JH , "masterpairs.csv")

############# READ IN HERE IF LINES 1 - 141 already ceated

library(dplyr)
#library(leaflet)
#library(htmltools)
#library(nngeo)
library(ggplot2)
library(rcompanion)
library(tidyverse)
library(lme4)
library(lmerTest)

###  

# read in shapefile
setwd("~/NPS")


master.JH <- read.csv("masterpairs.csv", header = TRUE)

head(master.JH)
## read it in here

## JDB - can use this to calculate differences between the two
## recal NEAR is w/o hemiparasites 

## richness/evenness with hemiparaistes
#master.JH$HERBRICH.diff <- with(master.JH, NOHEMIRICH.IN - NOHEMIRICH.NEAR) #
#master.JH$HERBEVEN.diff <- with(master.JH, NOHEMIEVEN.IN - NOHEMIEVEN.NEAR)
#master.JH$HERBEVEN.0.diff <- with(master.JH, NOHEMIEVEN.0.IN - NOHEMIEVEN.0.NEAR)  ## makign the "NA's (see below) = 0"

## richness/evennes with hemiparasites
master.JH$HERBHEMI.diff <- with(master.JH, HERBEVEN.IN - HERBEVEN.NEAR)
master.JH$HERBRICH.diff <- with(master.JH, HERBRICH.IN - HERBRICH.NEAR)
master.JH$TREEEVEN.diff <- with(master.JH, TREEEVEN.IN - TREEEVEN.NEAR)
master.JH$TREERICH.diff <- with(master.JH, TREERICH.IN - TREERICH.NEAR)

## all of hte metrics:
master.JH$ALLRICH.diff <- with(master.JH, ALLRICH.IN - ALLRICH.NEAR)
master.JH$ALLEVEN.diff <- with(master.JH, ALLEVEN.IN - ALLEVEN.NEAR)


# positive values = greater abundance where hemiparasites present

## T- tests


library(rcompanion)

## herbaceous with hemiparasites

herbhemieven.t <- with(master.JH, t.test(HERBHEMI.diff)); herbhemieven.t 
herbhemirich.t <- with(master.JH, t.test(HERBRICH.diff)); herbhemirich.t 

## tree
treeeven.t <- with(master.JH, t.test(TREEEVEN.diff)); treeeven.t
treerich.t <- with(master.JH, t.test(TREERICH.diff)); treerich.t


## all of it
alleven.t <- with(master.JH, t.test(ALLEVEN.diff)); alleven.t #
allrich.t <- with(master.JH, t.test(ALLRICH.diff)) ; allrich.t

## are the response variables normally distributed
plotNormalHistogram(master.JH$HERBRICH.diff)
plotNormalHistogram(master.JH$HERBEVEN.diff)
plotNormalHistogram(master.JH$ALLRICH.diff)
plotNormalHistogram(master.JH$ALLEVEN.diff)

head(master.JH)

## with means in dashed blue: 


##  can also pivot longer and use GLM

## 



head(master.JH)
master.JH <- master.JH %>% select(unique.pair, ALLRICH.IN:NOHEMIAB.IN, ALLRICH.NEAR:NOHEMIAB.NEAR , Area.IN, Area.NEAR)

summary(master.JH$Area.IN)
summary(master.JH$Area.NEAR)
head(master.JH)

master.JH <- master.JH%>% filter(!is.na(Area.IN), !is.na(Area.NEAR))
master.long <- master.JH %>% 
  pivot_longer(cols = -unique.pair, 
               names_to = c(".value", "in_near"), 
               names_sep = "\\.", values_drop_na = TRUE)

head(master.long)

master.long$HEMI.status = ifelse(master.long$HEMIAB >0, "Present", "Absent")


master.long <- master.long %>% mutate(ALLEVEN.logit = car::logit(ALLEVEN),
                       HERBEVEN.logit = car::logit(HERBEVEN),
                       TREEEVEN.logit = car:: logit(TREEEVEN),
                       ALLRICH.log = log(ALLRICH),
                       HERBRICH.log = log(HERBRICH),
                       TREERICH.log = log(TREERICH))

master.long$HERBRICH.log[is.finite(master.long$HERBRICH.log) == F] <- 0

master.long$TREERICH.log[is.finite(master.long$TREERICH.log) == F] <-0


#### MODELS 

alleven.mod <- lmer(ALLEVEN.logit ~ Area + HEMI.status + (1|unique.pair), master.long ); summary(alleven.mod)

plot(alleven.mod)

alleven.mod.n <- lmer(ALLEVEN.logit ~ Area + (1|unique.pair), master.long )
tab_model(alleven.mod)
deltaAIC = AIC(alleven.mod.n) - AIC(alleven.mod); deltaAIC ## 



alleven.mod <- lmer(ALLEVEN.logit ~ HEMI.status + (1|unique.pair), master.long ); summary(alleven.mod) ## for the graph. without area

ALLEVEN.pred <- data.frame(HEMI.status = master.long$HEMI.status)

ALLEVEN.pred <- master.long %>%
  mutate(ALLEVEN.logit = predict(alleven.mod, newdata = ALLEVEN.pred, re.form = NA),
         ALLEVEN = inv.logit(ALLEVEN.logit),
         HEMISTATUS = HEMI.status)

alleven.box <- ggboxplot(master.long, x = "HEMI.status", y = "TREEEVEN")+theme_bw(); alleven.box


alleven.dens <- ggplot(master.JH) +
  xlim(c(-0.7,0.7))+
  geom_density(aes(x = ALLEVEN.diff), fill = "#00AFBB", color = "#00AFBB", alpha = 0.3)+
  geom_vline(xintercept= inv.logit(.110), color = "black", linetype = "dashed")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # panel.border = element_blank(),
        axis.line.x = element_blank(),
        text = element_text(size = 10, family = "Helvetica"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        strip.background = element_rect(color  = "white",fill = "white"), ## takes away the facet borders 
        legend.position = "")+
  geom_vline(xintercept = 0, color = "black", linetype = "solid")+
  labs(x = " ", y = "All Plants", title = "Evenness");alleven.dens




all.even <- ggplot(data = master.long, aes(x = HEMI.status, y = ALLEVEN)) +
  geom_point(alpha = 0.3, color = "#00AFBB") +
  geom_line(data = ALLEVEN.pred, color = "black", size = 1)+
 # scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
  #  scale_x_log10(labels = comma)+ ## what the eff how do I make this not the scientific notation
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5))+
  labs(x = expression(""),
       y = expression("All Plants"),
       subtitle = "Evenness"); all.even
#ggsave("abundance.herb.tiff", path = "Figures/")

herbeven.mod <- lmer(HERBEVEN.logit ~ Area + HEMI.status + (1|unique.pair), master.long); summary(herbeven.mod)
herbeven.mod.n <- lmer(HERBEVEN.logit ~ Area + (1|unique.pair), master.long )
deltaAIC = AIC(herbeven.mod.n) - AIC(herbeven.mod); deltaAIC ## qhoa
tab_model(herbeven.mod)

treeeven.mod <- lmer(TREEEVEN.logit ~ Area + HEMI.status + (1|unique.pair), master.long ); summary(treeeven.mod)
treeeven.mod.n <- lmer(TREEEVEN.logit  ~ Area  + (1|unique.pair), master.long )
deltaAIC = AIC(treeeven.mod.n) - AIC(treeeven.mod); deltaAIC ## qhoa
tab_model(treeeven.mod)

allrich.mod <- lmer(ALLRICH.log ~ Area + HEMI.status + (1|unique.pair), master.long); summary(allrich.mod)
allrich.mod.n <- lmer(ALLRICH.log ~ Area + (1|unique.pair), master.long); summary(allrich.mod.n)
deltaAIC = AIC(allrich.mod.n ) - AIC(allrich.mod); deltaAIC

allrich.t <- with(master.JH, t.test(allrich.diff)) ## area lowe
tab_model(allrich.mod)

herbrich.mod <- lmer(HERBRICH.log ~ Area + HEMI.status + (1|unique.pair), master.long); summary(herbrich.mod)
herbrich.mod.n <- lmer(HERBRICH.log ~ Area + (1|unique.pair), master.long); summary(herbrich.mod.n)
deltaAIC = AIC(herbrich.mod.n) - AIC(herbrich.mod); deltaAIC
tab_model(herbrich.mod)
herbrich.t <- with(master.JH, t.test(HERBRICH.diff)) ## area lower


treerich.mod <- lmer(TREERICH.log ~ Area + HEMI.status + (1|unique.pair), master.long); summary(treerich.mod)
treerich.mod.n <- lmer(TREERICH.log ~ Area + (1|unique.pair), master.long); summary(treerich.mod.n)
deltaAIC = AIC(treerich.mod.n ) - AIC(treerich.mod); deltaAIC
tab_model(treerich.mod)
treerich.t <- with(master.JH, t.test(treerich.diff)) ## area lowe


## GRAPHING BOX PLOTS PAPER FIGURE


ggplot(master.long, aes(x = HEMI, y = ALLEVEN, group = unique.pair))+ ## forgot how to link the plots together...
  geom_point() +
  geom_line()+
  theme_bw()

rich <- ggplot(master.long, aes(x = HEMI.status, y = ALLRICH)) +
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#00c1cf", "#00AFBB"))+## how do you change hte mena line to work?
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "")+
  labs(y = expression(""),
       x = expression ("")); rich



tree.rich <- ggplot(master.long, aes(x = HEMI.status, y = TREERICH)) +
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#9ccdc1", "#69b3a2"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        legend.position = "")+
  labs(y = expression(""),
       x = expression ("")); tree.rich

even<-ggplot(master.long, aes(x = HEMI.status, y = ALLEVEN)) +
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#00c1cf", "#00AFBB"))+
  geom_boxplot(aes(color = HEMI.status))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "")+
  labs(y = expression("All Plants"),
       x = expression (""));even

herb.rich<-ggplot(master.long, aes(x = HEMI.status, y = HERBRICH)) +
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#fd8453", "#FC4E07"))+
  geom_boxplot(aes(color = HEMI.status))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "")+
  labs(y = expression(""),
       x = expression (""));herb.rich


herb.even<-ggplot(master.long, aes(x = HEMI.status, y = HERBEVEN)) +
  geom_boxplot(aes(color = HEMI.status))+
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#fd8453", "#FC4E07"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "")+
  labs(y = expression("Herbaceous Plants"),
       x = expression (""));herb.even

tree.even <- ggplot(master.long, aes(x = HEMI.status, y = TREEEVEN)) +
  geom_boxplot(aes(color = HEMI.status))+
  theme_bw()+
  geom_boxplot(aes(color = HEMI.status))+
  scale_color_manual(values = c( "#9ccdc1", "#69b3a2"))+
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        legend.position = "")+
  labs(y = expression("Woody Plants"),
       x = expression ("")); tree.even 

richness <- ggarrange(rich, herb.rich, tree.rich, nrow = 3, ncol= 1)
richness <- annotate_figure(richness, top = text_grob("Richness"))
evenness <- ggarrange(even, herb.even, tree.even, nrow = 3, ncol = 1)
evenness <- annotate_figure(evenness, top = text_grob("Evenness")); evenness

all <- ggarrange(evenness, richness); all
all <- annotate_figure(all, bottom = text_grob("Hemiparasite Presence")); all

all + ggsave("Figures/paper/boxplots.tiff", dpi = 300, height = 6, width = 6)






