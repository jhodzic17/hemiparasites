
## Abunaunce analysis
## Jasna Hodzic & Ian Pearse & Jon Bakker

#packages


library(lme4)
library(ggeffects)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lmerTest)
library(piecewiseSEM)
library(lmtest)
library(redres)
library(rcompanion)
library(gtools)
library(car)
library(scales)
require(scales)
## to download  redres ( to see diagnostic graphs in a Shiny app )
## devtools::install_github("goodekat/redres")

#working directory
setwd("~/NPS")


data <- read.csv("PARKHEMIFINAL/parkhemi.eco.csv", header = TRUE); summary(data)
area <- read.csv("plots_areasampled.csv", header = TRUE)  ## plot size from Eve

area <- area %>% select(SITEPLOT, Area_sampled)
data <- merge(data, area)


head(data)

############ Analysis ################


## transforming variables
data <-data %>% mutate(NOHEMIEVEN.logit = NOHEMIEVEN,
                             TREEEVEN.logit = TREEEVEN,
                             ALLEVEN.logit =ALLEVEN,
                             HEMIAB.log = log(HEMIAB),
                             ALLEVEN.logit = car::logit(ALLEVEN),
                             HERBEVEN.logit = car::logit(HERBEVEN),
                             NOHEMIEVEN.logit = car::logit(NOHEMIEVEN.logit),
                       TREEEVEN.logit = car::logit(TREEEVEN),
                       ALLRICH.log = log(ALLRICH),
                       HERBRICH.log = log(HERBRICH),
                       TREERICH.log = log(TREERICH))

data$HERBRICH.log[is.finite(data$HERBRICH.log)==F]<-0
data$TREERICH.log[is.finite(data$TREERICH.log) == F] <-0
#data$HEMIAB.log[is.finite(data$HEMIAB.log == F)] <- 0
data.ab<-data %>% filter(HEMIAB > 0) 
summary(data.ab)

data.ab <- data %>% filter(HEMIAB > 0 & Area_sampled > 398 & Area_sampled < 401) ## plots with hemiparasites and ~ 400m2
summary(data.ab)


############### Analysis #################

### Is hemiparasite abundance correlated with an increase in richness or evenness?

### ALLEVEN 

model.alleven <-lmer(ALLEVEN.logit ~ HEMIAB.log  + (1|PARK) + (1|L4_KEY) + (1|YEAR), data.ab)
model.alleven.n <-lmer(ALLEVEN.logit ~ 1 + (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab) ## null without HEMIAB

deltaAIC = AIC(model.alleven) - AIC(model.alleven.n); deltaAIC
summary(model.alleven) 
confit()
library(performance)
r2_nakagawa(model.alleven)
##
model.alleven <-lmer(ALLEVEN.logit ~ HEMIAB.log  + (1|PARK) + (1|L4_KEY) + (1|YEAR), data.ab) ## taking area out for the graph

ALLEVEN.pred <- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                            to = max(data.ab$HEMIAB.log),
                                            by = 0.05))
ALLEVEN.pred <- ALLEVEN.pred %>%
  mutate(ALLEVEN.logit = predict(model.alleven, newdata = ALLEVEN.pred, re.form = NA),
         ALLEVEN = inv.logit(ALLEVEN.logit),
         HEMIAB = exp(HEMIAB.log))

library(ggpubr)

all.even <- ggplot(data = data.ab, aes(x = HEMIAB, y = ALLEVEN)) +
  geom_point(alpha = 0.3, color = "#00AFBB") +
  geom_line(data = ALLEVEN.pred, color = "black", size = 1)+
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
       # axis.line.x = element_blank(),
      #  axis.line.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5))+
  labs(x = expression(""),
       y = expression("All Plants"),
       subtitle = "Evenness"); all.even
#ggsave("abundance.herb.tiff", path = "Figures/")


### HERBEVEN (

model.herb <- lmer(HERBEVEN.logit ~  HEMIAB.log + (1|PARK) + (1|L4_KEY) +  (1|YEAR) ,data.ab) 

model.herb.n <- lmer(HERBEVEN.logit ~  (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab) ## null wthout hemi abundaunce

summary(model.herb)
AIC(model.herb, model.herb.n) ## 
deltaAIC = AIC(model.herb) - AIC(model.herb.n); deltaAIC
r2_nakagawa(model.herb)

## without area_sampled; for the graphs below

model.herb <- lmer(HERBEVEN.logit ~  HEMIAB.log + (1|PARK) + (1|L4_KEY) +  (1|YEAR) ,data.ab) ## without area

## graph
HERBEVEN.pred <- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                             to = max(data.ab$HEMIAB.log),
                                             by = 0.05))
HERBEVEN.pred <- HERBEVEN.pred %>%
  mutate(HERBEVEN.logit = predict(model.herb, newdata = HERBEVEN.pred, re.form = NA),
         HERBEVEN = inv.logit(HERBEVEN.logit),
         HEMIAB = exp(HEMIAB.log))




herb.even <- ggplot(data = data.ab, aes(x = HEMIAB, y = HERBEVEN)) +
  geom_point(color = "#FC4E07", alpha = 0.3) +
 # geom_density( fill="#69b3a2", color="#69b3a2", alpha=0.4)+
 geom_line(data = HERBEVEN.pred, color = "black", size = 1)+
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
 ## what the eff how do I make this not the scientific notation
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
       # axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
      #  axis.line.y = element_blank(),
        text = element_text(size = 10))+
  labs(x = expression(""),
       y = expression("Herbaceous Plants")); herb.even
#ggsave("abundance.herb.tiff", path = "Figures/")



############  Woody Plants
summary(data.ab$HERBEVEN.logit)
summary(data.ab$TREEEVEN)


data.ab$TREEEVEN_LOGIT <- car::logit(data.ab$TREEEVEN, percents = FALSE, adjust = 0.025)

model.tree <-lmer(TREEEVEN.logit ~  HEMIAB.log + (1|PARK) + (1 | L4_KEY) + (1|YEAR), data.ab)
model.tree.n <- lmer(TREEEVEN.logit ~  (1|PARK) + (1 | L4_KEY) + (1|YEAR), data.ab) ## again park doesn't do anything
summary(model.tree)
 
deltaAIC = AIC(model.tree) - AIC(model.tree.n); deltaAIC ## AIC lowered by taking out hemiparasites

launch_redres(mod1hemi.tree) ## assumptions
r2_nakagawa(model.tree)
##
model.tree <-lmer(TREEEVEN.logit ~  HEMIAB.log + (1|PARK) + (1 | L4_KEY) + (1|YEAR), data.ab) ## for graphing; removing area


HEMIEVEN.pred.tr <- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                             to = max(data.ab$HEMIAB.log),
                                             by = 0.05))
HEMIEVEN.pred.tr <- HEMIEVEN.pred.tr %>%
  mutate(TREEEVEN.logit = predict(model.tree, newdata = HEMIEVEN.pred.tr, re.form = NA),
         TREEEVEN = inv.logit(TREEEVEN.logit),
         HEMIAB = exp(HEMIAB.log))



tree.even <- ggplot(data = data.ab, aes(x = HEMIAB, y = TREEEVEN)) +
  geom_point(color = "#69b3a2", alpha = 0.3) +
  #geom_line(data = HEMIEVEN.pred.tr, color = "black", size = 1)+
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
#  scale_x_log10(labels = comma)+ ## what the eff how do I make this not the scientific notation
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = expression(""),
       y = expression("Woody Plants")); tree.even



s


############### Richness ~ Abundance #############3

##all plants
model.allrich  <-lmer(ALLRICH.log ~  HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)
model.allrich.n <-lmer(ALLRICH.log ~ (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)
deltaAIC = AIC(model.allrich) - AIC(model.allrich.n); deltaAIC ## null preferred 
summary(model.allrich)

deltaAIC =  AIC(model.allrich) - AIC(model.allrich.n); deltaAIC
r2_nakagawa(model.allrich)


anova(mod1rich.all, mod1rich.n)
summary(mod1rich.all) 
tab_model(mod1rich.all)
launch_redres(mod1rich.all) ## nice
anova(mod1rich.all) 
## faint signal; i.e. vines, etc. doing something (those features not in herb or tree. Makes sense for richness)
## though graph shows richness decreasing with abundance; opposite of presence test



## graph  

model.allrich  <-lmer(ALLRICH.log ~  HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)

HEMIRICH.pred.all <- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                                 to = max(data.ab$HEMIAB.log),
                                                 by = 0.05))

HEMIRICH.pred.all<- HEMIRICH.pred.all %>%
  mutate(ALLRICH.log = predict(model.allrich, newdata = HEMIRICH.pred.all, re.form = NA),
         ALLRICH = exp(ALLRICH.log),
         HEMIAB = exp(HEMIAB.log)) ## n ot backtransforemd for visual


rich.all<- ggplot( data.ab, aes(x = HEMIAB, y = ALLRICH)) +
  geom_point( color = "#00AFBB", alpha = 0.3) +
#  geom_line(data = HEMIRICH.pred.all, color = "red", size = 1)+
  scale_x_log10(labels = comma)+
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
      #  axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
       # axis.line.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5))+
  labs(x = expression(),
       y = expression(),
       subtitle = "Richness");rich.all
# ggsave("richall.norm.tiff", path = "Figures/", height = 5, width = 7)


## herb
model.herbrich <- lmer(HERBRICH.log ~ HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)
model.herbrich.n <- lmer(HERBRICH.log ~ (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)
deltaAIC = AIC(model.herbrich) - AIC(model.herbrich.n); deltaAIC
summary(model.herbrich)
r2_nakagawa(model.herbrich)

HEMIRICH.pred <- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                             to = max(data.ab$HEMIAB.log),
                                             by = 0.05))

## graph
model.herbrich <- lmer(HERBRICH.log ~ HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR) , data.ab)

HEMIRICH.pred <- HEMIRICH.pred%>%
  mutate(HERBRICH.log = predict(model.herbrich, newdata = HEMIRICH.pred, re.form = NA),
         HERBRICH = exp(HERBRICH.log),
         HEMIAB = exp(HEMIAB.log)) ## n ot backtransforemd for visual



rich.herb<- ggplot( data.ab, aes(x = HEMIAB, y = HERBRICH)) +
  geom_point(color = "#FC4E07", alpha= 0.3) +
  #geom_line(data = HEMIRICH.pred, color = "black", size = 1)+
  scale_x_log10(labels = comma)+
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank())+
        #axis.line.x = element_blank(),
       # axis.line.y = element_blank())+
  labs(x = expression(),
       y = expression());rich.herb
# ggsave("richall.norm.tiff", path = "Figures/", height = 5, width = 7)



summary(data.ab$NOHEMIRICH)


## transformations


data.ab <- data.ab %>% mutate(NOHEMIRICH.log = log(NOHEMIRICH),
                   TREERICH.log = log(TREERICH),
                   ALLRICH.log = log(ALLRICH),
                   HERBRICH.log = log(HERBRICH))



## trees
data.ab$TREERICH.log[is.finite(data.ab$TREERICH.log) == F] <-0

model.rich.tree <- lmer(TREERICH.log ~  HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR), data.ab)
model.rich.tree.n <- lmer(TREERICH.log ~ (1|PARK) + (1|L4_KEY) + (1|YEAR), data.ab)
summary(model.rich.tree) ## rando effects nto accounting for much variation.. 
deltaAIC = AIC(model.rich.tree) - AIC(model.rich.tree.n);deltaAIC

r2_nakagawa(model.rich.tree)
TREERICH.pred<- data.frame(HEMIAB.log = seq(from = min(data.ab$HEMIAB.log), 
                                                 to = max(data.ab$HEMIAB.log),
                                                 by = 0.05))


## graph
model.rich.tree <- lmer(TREERICH.log ~ HEMIAB.log + (1|PARK) + (1|L4_KEY) + (1|YEAR), data.ab)

TREERICH.pred<- TREERICH.pred %>%
  mutate(TREERICH.log = predict(model.rich.tree, newdata = TREERICH.pred, re.form = NA),
         TREERICH = exp(TREERICH.log),
         HEMIAB = exp(HEMIAB.log)) ## n ot backtransforemd for visual


rich.tree<- ggplot( data.ab, aes(x = HEMIAB, y = TREERICH)) +
  geom_point(alpha = 0.3,  color="#69b3a2") +
  scale_x_log10(labels = scales::number_format(accuracy = 0.01))+
  #geom_line(data = TREERICH.pred, color = "red", size = 1)+
 # scale_x_log10(labels = comma)+
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        plot.subtitle = element_text(hjust = 0.5))+
  labs(x = expression(),
       y = expression());rich.tree
# ggsave("richall.norm.tiff", path = "Figures/", height = 5, width = 7)



library(ggpubr)
graph.ab <- ggarrange(all.even, rich.all,  herb.even,rich.herb, tree.even, rich.tree, ncol = 2, nrow = 3, align = "v"); graph.ab ## do the same here as you did with the presence graph 
graph.ab<- annotate_figure(graph.ab, bottom = text_grob("Hemiparasite Abundance"))
graph.ab + ggsave("Figures/paper/abundancemod400.tiff", height = 7, width = 9, dpi = 300)

