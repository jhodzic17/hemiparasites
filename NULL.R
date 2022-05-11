###analysis of relationship between hemiparasite abundance / presence and the evenness of plant communities in the US National Parks using NPS inventory and monitoring data.

###R libraries and functions

library(vegan)
library (dplyr)
library(lme4)
library(magrittr)
library(car)


###combine data from each National Park into one big file with all I&M data from NPS
# setwd("C:/Users/ipearse/Desktop/hemiparasites/analysis/Formatted Cover Data Per Site Mar23/park data")
# files<-list.files()
# read.csv("clean_ZION_cover.csv")->aaa
# parkall<-aaa[0,]
# for(i in 1:length(files)){
	# parki<-read.csv(files[i])
	# parkall<-rbind(parkall,parki)
# }


##set working directory and load data



setwd("~/NPS")
area <- read.csv("plots_areasampled.csv", header = TRUE) %>%
  dplyr::select(SITEPLOT, Area_sampled) 
park <- read.csv("otherdatasheets/park.eco.csv", header = TRUE)
hemi<-read.csv("otherdatasheets/hemiparasite_USDAplantcodes.csv")
park <- merge(area, park)


###create a list of species that are not hemiparasites, but fall within hemiparasite global ranks to draw communities from
spdat_ab<-tapply(park$Pct_Cov, park$Species, sum)
spdat_ab<-sort(spdat_ab, decreasing=T)
absran<-range(which(names(spdat_ab)%in%hemi$Accepted.Symbol))
spoth<-spdat_ab[absran[1]:absran[2]]
spoth<-spoth[which(!(names(spoth)%in%hemi$Accepted.Symbol))]

####make random assemblages of non-hemiparasites with the same number of species as there are hemiparasites in NPS dataset drawn from those species that fall within hemiparasite rank abundance.
n_hemi<-sum(unique(hemi$Accepted.Symbol)%in%names(spdat_ab))  ##there are 97 hemiparasite 'taxa' (including non-specified genera)
rass<-sample(spoth, n_hemi)


###make a reduced dataframe where each row is a plot within site (i.e one sampling place within each park)
park$SITEPLOT<-paste(park$Site, park$Plot)
park <- park %>% filter(!Growth_habit %in% c("Nonvascular", "Lichenous"))

unsite<-unique(park$SITEPLOT)
parkhemi<-data.frame(SITEPLOT=unsite, PARK=rep(NA,length(unsite)), PLOT=rep(NA,length(unsite)),LONG=rep(NA,length(unsite)),LAT=rep(NA,length(unsite)),YEAR=rep(NA,length(unsite)), ECO3=rep(NA,length(unsite)),ECO4=rep(NA,length(unsite)),
                     Area_sampled = rep(NA, length(unsite)), ##all the useful info from 'park'
	ALLRICH=rep(NA,length(unsite)),ALLEVEN=rep(NA,length(unsite)), ALLSHANON=rep(NA,length(unsite)),ALLAB=rep(NA,length(unsite)),ALLINVAB=rep(0,length(unsite)),ALLNATRICH=rep(0,length(unsite)),ALLNATAB=rep(0,length(unsite)),
	TREERICH = rep(0, length(unsite)), TREEEVEN = rep(NA, length(unsite)), TREEAB = rep(0, length(unsite)), TREESHANON = rep(NA, length(unsite)),
	HERBRICH=rep(0,length(unsite)),HERBEVEN=rep(NA,length(unsite)),HERBSHANON=rep(NA,length(unsite)),HERBAB=rep(0,length(unsite)),HERBNATRICH=rep(0,length(unsite)),HERBNATAB=rep(0,length(unsite)),HERBINVAB=rep(0,length(unsite)),
	HEMIAB=rep(0,length(unsite)),HEMIRICH=rep(0,length(unsite)),Comandra=rep(0,length(unsite)),Castilleja=rep(0,length(unsite)),Dasistoma=rep(0,length(unsite)),Orthocarpus=rep(0,length(unsite)),Pedicularis=rep(0,length(unsite)))
hnam<-c(as.character(hemi$Accepted.Symbol), as.character(hemi$Synonym.Symbol))
park$Species<-as.character(park$Species)
park$Species[which(is.na(park$Species))]<-"other"
for(i in 1:nrow(parkhemi)){
	flush.console()
	print(paste(i, "of",nrow(parkhemi)))
	parki<-park[which(park$SITEPLOT==parkhemi$SITEPLOT[i]),]
	parkhemi$PARK[i]<-parki$Site[1]
	parkhemi$PLOT[i]<-parki$Plot[1]
	parkhemi$LONG[i]<-parki$Long[1]
	parkhemi$LAT[i]<-parki$Lat[1]
	parkhemi$YEAR[i]<-parki$Year[1]
	parkhemi$ECO3[i]<-parki$L3_KEY[1]
	parkhemi$ECO4[i]<-parki$L4_KEY[1]
	parkhemi$Area_sampled[i] <- parki$Area_sampled[1]
	parkvec<-tapply(parki$Pct_Cov, as.character(parki$Species),sum)
	parkhemi$ALLRICH[i]<-specnumber(parkvec)
	parkhemi$ALLSHANON[i]<-diversity(parkvec,index='shannon')
	parkhemi$ALLEVEN[i]<-parkhemi$ALLSHANON[i]/log(parkhemi$ALLRICH[i])
	parkhemi$ALLAB[i]<-sum(parkvec)
	if("I" %in% parki$Exotic) parkhemi$ALLINVAB[i]<-sum(parki$Pct_Cov[which(parki$Exotic=="I")])
	if("N" %in% parki$Exotic) parkhemi$ALLNATAB[i]<-sum(parki$Pct_Cov[which(parki$Exotic=="N")])
	if("N" %in% parki$Exotic) parkhemi$ALLNATRICH[i]<-length(parki$Pct_Cov[which(parki$Exotic=="N")])
	parkiherb<-parki[which(parki$Growth_habit %in% c("Graminoid","Forb/herb","Subshrub" , "Vine")),]  ###changed this line to be based not on stratum, but plant type
	if (nrow(parkiherb)>0){ 
		herbvec<-tapply(parkiherb$Pct_Cov, as.character(parkiherb$Species),sum)
		parkhemi$HERBRICH[i]<-specnumber(herbvec)
		parkhemi$HERBSHANON[i]<-diversity(herbvec,index='shannon')
		parkhemi$HERBEVEN[i]<-parkhemi$HERBSHANON[i]/log(parkhemi$HERBRICH[i])
		parkhemi$HERBAB[i]<-sum(herbvec)
		if("I" %in% parkiherb$Exotic) parkhemi$HERBINVAB[i]<-sum(parkiherb$Pct_Cov[which(parkiherb$Exotic=="I")])
		if("N" %in% parkiherb$Exotic) parkhemi$HERBNATAB[i]<-sum(parkiherb$Pct_Cov[which(parkiherb$Exotic=="N")])
		if("N" %in% parkiherb$Exotic) parkhemi$HERBNATRICH[i]<-length(parkiherb$Pct_Cov[which(parkiherb$Exotic=="N")])
	}
	parkitree<-parki[which(parki$Growth_habit %in% c("Shrub", "Tree")),] 
	if(nrow(parkitree) >0){
	  treevec<-tapply(parkitree$Pct_Cov, as.character(parkitree$Species),sum)
	  parkhemi$TREERICH[i]<-specnumber(treevec)
	  parkhemi$TREESHANON[i] <- diversity(treevec, index = 'shannon')
	  parkhemi$TREEEVEN[i]<-parkhemi$TREESHANON[i]/log(parkhemi$TREERICH[i])
	  parkhemi$TREEAB[i]<-sum(treevec)
	}
	phemi<-parki[which(as.character(parki$Species)%in%hnam),]
	if(nrow(phemi)>0){
		hemivec<-tapply(phemi$Pct_Cov, as.character(phemi$Species),sum)
		parkhemi$HEMIAB[i]<-sum(hemivec)
		parkhemi$HEMIRICH[i]<-length(hemivec)
	}
} 

head(parkhemi)
## transformation
parkhemi$ALLEVEN_LOGIT<-car::logit(parkhemi$ALLEVEN, percents = FALSE, adjust = 0.025)
parkhemi$HERBEVEN_LOGIT<-car::logit(parkhemi$HERBEVEN, percents = FALSE, adjust = 0.025)
parkhemi$TREEEVEN_LOGIT <- car::logit(parkhemi$TREEEVEN, percents = FALSE, adjust = 0.025)

parkhemi$ALLRICH_LOG<-log(parkhemi$ALLRICH)
parkhemi$HERBRICH_LOG<-log(parkhemi$HERBRICH)
parkhemi$HERBRICH_LOG[is.finite(parkhemi$HERBRICH_LOG)==F]<-0
parkhemi$TREERICH_LOG <- log(parkhemi$TREERICH)
parkhemi$TREERICH_LOG[is.finite(parkhemi$TREERICH_LOG)==F]<-0
#parkhemi$HEMIAB2<-parkhemi$HEMIAB
parkhemi$HEMIAB.log<-log(parkhemi$HEMIAB)



 write.csv(parkhemi, "parkhemirandom.csv")


################# ABUNDANCE ANALYSIS
setwd("~/NPS")
parkhemi.original <- read.csv("parkhemirandom.csv"); 
parkhemi <- parkhemi.original

library(vegan)
library (dplyr)
library(lme4)
library(tidyverse) ##JDB - so that tidyr functions also available
library(car)


## let us restrict it to 400m2 , rounding up 399

parkhemi <- parkhemi %>% filter(Area_sampled > 398 & Area_sampled < 401) ## allowing these plots in

hemiab <- parkhemi[parkhemi$HEMIAB >0,]
summary(hemiab$Area_sampled)
plot.number<-hemiab %>% group_by(Area_sampled) %>% count()

# let's select relevant columns so this runs faster, dont need Area_sampled anymore
parkhemi <- parkhemi %>%
  select(c(HEMIAB.log, ALLEVEN_LOGIT, HERBEVEN_LOGIT, TREEEVEN_LOGIT, ALLRICH_LOG,
           HERBRICH_LOG, TREERICH_LOG, ECO4, YEAR, PARK, SITEPLOT))

park <- read.csv("otherdatasheets/park.eco.csv", header = TRUE)
hemi<-read.csv("otherdatasheets/hemiparasite_USDAplantcodes.csv")
sumpark <- park %>%
  filter(SITEPLOT %in% parkhemi$SITEPLOT)
##JDB - filter to 400m^2 plots

#### 
###create a list of species that are not hemiparasites, but fall within hemiparasite global ranks to draw communities from
spdat_ab<-tapply(park$Pct_Cov, park$Species, sum)
spdat_ab<-sort(spdat_ab, decreasing=T)
absran<-range(which(names(spdat_ab)%in%hemi$Accepted.Symbol))
spoth<-spdat_ab[absran[1]:absran[2]]
spoth<-spoth[which(!(names(spoth)%in%hemi$Accepted.Symbol))]


####make random assemblages of non-hemiparasites with the same number of species as there are hemiparasites in NPS dataset drawn from those species that fall within hemiparasite rank abundance.
n_hemi <- sum(unique(hemi$Accepted.Symbol) %in% names(spdat_ab))  ##there are 97 hemiparasite 'taxa' (including non-specified genera)
##JDB - only 78 taxa once we filter to 400m^2 plots

## how many unique taxa?
park$Species <- as.factor(park$Species)
nlevels(park$Species)
##Create reduced dataframe using random subsets of data to describe the presence and abundance of non-hemiparasites
n_rand<-1000
parkhemi_rand <- c()
for(j in 1:n_rand) {
	print(j); set.seed(j)
	rass <- sample(spoth, n_hemi)  #grab random sample of non-hemiparasites
	##JDB - drop non-selected species and calculated abundance of randomly selected non-hemiparasites
	rass_cov <- park %>%
	  filter(Species %in% names(rass)) %>%
	  group_by(SITEPLOT) %>%
	  summarize(Ab = sum(Pct_Cov))
	rass_cov <- rass_cov %>% mutate(Log_Ab = log(Ab), Ab = NULL, RandomSubset = paste('Ran_Ab',j, sep=""))
	parkhemi_rand <- rbind(parkhemi_rand, rass_cov)
 }
parkhemi_rand <- parkhemi_rand %>%
  pivot_wider(names_from = RandomSubset, values_from = Log_Ab)
##JDB - 11490 SITEPLOTs here but 12498 meet size criteria
##JDB - therefore, some never contained any of the randomly selected species

parkhemi_full <- merge(x = parkhemi[ , c("SITEPLOT", "PARK", "YEAR", "ECO4", "ALLEVEN_LOGIT",
                                         "ALLRICH_LOG", "HERBEVEN_LOGIT", "HERBRICH_LOG",
                                         "TREEEVEN_LOGIT", "TREERICH_LOG", "HEMIAB.log")],
                       y = parkhemi_rand, by = "SITEPLOT", all.x = TRUE)
##JDB - calling columns by name


############# ABUNDANCE ANALYSIS #################
models <- c("ALLEVEN_LOGIT", "ALLRICH_LOG",
            "HERBEVEN_LOGIT", "HERBRICH_LOG",
            "TREEEVEN_LOGIT", "TREERICH_LOG")
model_results <- c()
for(j in 1:n_rand) {
	print(j)
  random.column <- paste("Ran_Ab", j, sep = "")
	parkhemi2 <- parkhemi_full[which(is.finite(parkhemi_full[ , random.column]) == T) , ] ##filtering to null plant > 0 abundance

	for(i in 1:length(models)) {
	  response <- models[i]
	  modform <- formula(paste(response,"~",random.column,"+(1|PARK)+(1|ECO4)+(1|YEAR)"))
	  aaa <- lmer(modform, data = parkhemi2)
	  model_results_i <- data.frame(RandomSubset = j,
	                                NumberPlots = nrow(parkhemi2),
	                                Response = response,
	                                Coefficient = fixef(aaa)[2],
	                                Singular = isSingular(aaa))
	  model_results <- rbind(model_results, model_results_i)
	}
}

##JDB - repeat same analysis using actual hemiparasite abundance values and add to model_results
parkhemi2 <- parkhemi_full[which(is.finite(parkhemi_full$HEMIAB.log) == T) , ] ##filtering to null plant > 0 abundance
for(i in 1:length(models)) {
  response <- models[i]
  modform <- formula(paste(response,"~HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR)"))
  aaa <- lmer(modform, data = parkhemi2)
  model_results_i <- data.frame(RandomSubset = "Actual",
                                NumberPlots = nrow(parkhemi2),
                                Response = response,
                                Coefficient = fixef(aaa)[2],
                                Singular = isSingular(aaa)) ## still would like to know the number of observations of each null plant ; i.e. how many plots have null plant ab >0 to generate the null model coeff
  model_results <- rbind(model_results, model_results_i)
}
rm(parkhemi2)

write.csv(model_results, "null_abunduance_results.csv")

model_results<-read.csv("null_abunduance_results.csv")
##JDB - number of random subsamples with singular solutions (high for trees!)
model_results %>% filter(Singular == TRUE) %>% group_by(Response) %>% count(Singular)
model_results %>% group_by(Response, Singular) %>% count(Response)
##JDB - calculate p-values for each response
model_results %>% 
  #filter(Singular == FALSE) # %>% ##JDB - uncomment this line to drop singular models (but p-values minimally affected)
  group_by(Response) %>%
  summarize(P = sum(Coefficient >= Coefficient[RandomSubset == "Actual"]) / length(Coefficient))
##JDB - herb evenness strongly significant; all evenness marginal; others not significant


# ## FIGURES 

##JDB - create function once here and then call as needed 
##JDB - allows easy changes to cascade through all figures
generic_graphing_function <- function(data, color, ylab, textual_note) {
  require(ggplot2)
  graphic <- ggplot(data = data, aes(x = Coefficient)) +
    geom_density(fill = color, color = color, alpha = 0.3) +
    xlim(c(-0.12, 0.2)) +
    geom_vline(aes(xintercept = Coefficient[RandomSubset == "Actual"]),
               color="black", linetype="dashed", size=1) +
    geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text = element_text(size = 10),
          # panel.border = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    labs(y = ylab,
         x = expression(""),
         subtitle = textual_note)
  graphic
}

alleven <- model_results %>%
  filter(Response == "ALLEVEN_LOGIT") %>%
  generic_graphing_function(data = .,
                            color = "#00AFBB",
                            ylab = expression("All Plants"),
                            textual_note = expression("p = 0.070"))

herbeven <- model_results %>%
  filter(Response == "HERBEVEN_LOGIT") %>%
  generic_graphing_function(data = .,
                            color = "#FC4E07",
                            ylab = expression("Herbaceous Plants"),
                            textual_note = expression("p < 0.001"))

treeeven <- model_results %>%
  filter(Response == "TREEEVEN_LOGIT") %>%
  generic_graphing_function(data = .,
                            color = "#69b3a2",
                            ylab = expression("Woody Plants"),
                            textual_note = expression("p = 0.914"))

allrich <- model_results %>%
  filter(Response == "ALLRICH_LOG") %>%
  generic_graphing_function(data = .,
                            color = "#00AFBB",
                            ylab = expression(""),
                            textual_note = expression("p = 0.453"))

herbrich <- model_results %>%
  filter(Response == "HERBRICH_LOG") %>%
  generic_graphing_function(data = .,
                            color = "#FC4E07",
                            ylab = expression(""),
                            textual_note = expression("p = 0.560"))

treerich <- model_results %>%
  filter(Response == "TREERICH_LOG") %>%
  generic_graphing_function(data = .,
                            color = "#69b3a2",
                            ylab = expression(""),
                            textual_note = expression("p = 0.196"))


## gather graphs together into one figure
require(ggpubr)
even <- ggarrange(alleven, herbeven, treeeven, nrow = 3, ncol = 1)
even <- annotate_figure(even, top = text_grob("Evenness")); even

rich <- ggarrange(allrich, herbrich, treerich, nrow = 3 , ncol = 1)
rich <- annotate_figure(rich, top = text_grob("Richness")); rich

graph.null <- ggarrange(even, rich); graph.null
graph.null <- annotate_figure(graph.null,    bottom = text_grob("Model Coefficients"), top = text_grob("Abundance",size = 14, vjust = 0.50))

graph.null + ggsave("Figures/null.ab.tiff", height = 7, width = 7, dpi = 300)
ggsave("null.png", height = 7, width = 7, dpi = 300) ##JDB - smaller file for sharing


##################################PAIRED PLOTS PRESENCE ANALYSIS############################################################3
#setwd("~/NPS")
parkhemi <- read.csv("parkhemirandom.csv")
 
parkhemi <- parkhemi %>% filter(Area_sampled > 398 & Area_sampled < 401) ## allowing these plots in
 
## let's select relevant columns so this runs faster, dont need Area_sampled anymore
 
parkhemi <- parkhemi %>% select(c(HEMIAB.log, HEMIAB, LAT, LONG, ALLEVEN_LOGIT, HERBEVEN_LOGIT, TREEEVEN_LOGIT, ALLRICH_LOG, HERBRICH_LOG, TREERICH_LOG, ECO4, YEAR, PARK, SITEPLOT))
 
 ###Add to reduced dataframe using random subsets of data to describe the presence and abundance of non-hemiparasites
 n_rand<-1000
 for(j in 1:n_rand){
   print(j); set.seed(j)
   rass<-sample(spoth, n_hemi)  #grab random sample of non-hemiparasites
   coln<-rep(NA,nrow(parkhemi))
   parkhemi<-cbind(parkhemi, coln)
   colnames(parkhemi)[(14+j)]<-paste('Ran_Ab',j, sep="") ## changing here to reflect the transformed variables
   parktr<-park
   parktr$Pct_Cov[which(!(parktr$Species%in%names(rass)))]<-0
   rass_cov<-tapply(parktr$Pct_Cov, parktr$SITEPLOT, sum)
   parkhemi[,(14+j)]<-log(rass_cov[parkhemi$SITEPLOT])## changing here to reflect the transformed variables
   #parkhemi[,(29+j)][which(is.finite(parkhemi[,(29+j)])==F)]<-0
 }
 


 ###make global plot pairs

 
 parkhemi_d <- parkhemi
 colnames(parkhemi)
 parkhemi_d<-parkhemi_d[,c(1, 2, 3, 4, 11, 13, 14)] ## all relevant columns for distane - plot, hemiab, lat, long, eco4
 colnames(parkhemi_d)
 parkhemi_d$pair1<-NA
 parkhemi_d$dist1<-NA
 
 for(i in 1:nrow(parkhemi)){
   print(i);set.seed(i)
   ploti<-parkhemi_d[i,]
   p_red<-parkhemi_d[which(parkhemi_d$ECO4==ploti$ECO4),]
   if(nrow(p_red)>1){
     p_red<-parkhemi_d[which(parkhemi_d$SITEPLOT!=ploti$SITEPLOT),]
     p_red$disti<-sqrt(((p_red$LAT-ploti$LAT)^2)+((p_red$LONG-ploti$LONG)^2))
     plotp<-p_red[which(p_red$disti==min(p_red$disti)),][1,]
     pairnam<-paste("P",i,sep="")
     if(is.na(ploti$pair1)==T & is.na(plotp$pair1)==T){
       parkhemi_d$pair1[i]<-pairnam
       parkhemi_d$dist1[i]<-plotp$disti
       parkhemi_d$pair1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-pairnam
       parkhemi_d$dist1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-plotp$disti
     }
     if(is.na(ploti$pair1)==F & plotp$disti<ploti$dist1){
       parkhemi_d$pair1[which(parkhemi_d$pair1==ploti$pair1)]<-NA
       parkhemi_d$pair1[i]<-pairnam
       parkhemi_d$dist1[i]<-plotp$disti
       parkhemi_d$pair1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-pairnam
       parkhemi_d$dist1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-plotp$disti
     }
     if(is.na(plotp$pair1)==F & plotp$disti<plotp$dist1){
       parkhemi_d$pair1[which(parkhemi_d$pair1==plotp$pair1)]<-NA
       parkhemi_d$pair1[i]<-pairnam
       parkhemi_d$dist1[i]<-plotp$disti
       parkhemi_d$pair1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-pairnam
       parkhemi_d$dist1[which(parkhemi_d$SITEPLOT==plotp$SITEPLOT)]<-plotp$disti
     }
   }
 }
 
 
 
 parkhemi$pair1<-parkhemi_d$pair1
 parkhemi$dist1<-parkhemi_d$dist1
 parkhemi<-parkhemi[which(is.na(parkhemi$pair1)==F),]
 parkhemi<-parkhemi[!parkhemi$dist1 == 0,]
 head(parkhemi)
 
# parkhemi <- read.csv("Figures/parkhemi_null_paired_transformed.final_1.csv") #
 parkhemi$pair1 <- as.factor(parkhemi$pair1)
 summary(parkhemi$pair1)
 
 write.csv(parkhemi, "PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/parkhemi.pairs.minuszero.csv", row.names=F)
 
#parkhemi <- read.csv("PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/parkhemi.pairs.csv", row.names=F)

 ##DO ANALYSIS OF REAL HEMIPARASITES IN PAIRED PLOTS
 ###
 hemidat<-parkhemi[which(parkhemi$HEMIAB>0),]
 nohemidat<-parkhemi[which(parkhemi$HEMIAB==0),]
 parkpair<-parkhemi[which(parkhemi$pair1%in%hemidat$pair1),]
 parkpair<-parkpair[which(parkpair$pair1%in%nohemidat$pair1),]
 parkpair$HEMIPRES<-0
 parkpair$HEMIPRES[which(parkpair$HEMIAB>0)]<-1
 
 
 mod1<-lm(ALLEVEN_LOGIT~HEMIPRES+pair1, data=parkpair); summary(mod1)
 mod2<-lm(ALLRICH_LOG~HEMIPRES+pair1, data=parkpair)
 mod3<-lm(HERBEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
 mod4<-lm(HERBRICH_LOG~HEMIPRES+pair1, data=parkpair)
 mod5<-lm(TREEEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
 mod6<-lm(TREERICH_LOG~HEMIPRES+pair1, data=parkpair)
 
 
 #mod1<-lmer(HERBEVEN~HEMIPRES+Area_sampled+(1|pair1), data=parkpair)
 #mod1n<-lmer(HERBEVEN~Area_sampled+(1|pair1), data=parkpair)
 
 ##DO ANALYSES OF NULL PLANT COMMUNITIES
 n_rand=1000
 null_coeff_tab_pair<-data.frame(ALLEVEN_LOGIT=rep(NA, n_rand), ALLRICH_LOG=rep(NA, n_rand),HERBEVEN_LOGIT=rep(NA, n_rand),HERBRICH_LOG=rep(NA, n_rand),TREEEVEN_LOGIT=rep(NA, n_rand),TREERICH_LOG=rep(NA, n_rand), 
                                 ALLEVEN_SE=rep(NA, n_rand), ALLRICH_SE=rep(NA, n_rand),HERBEVEN_SE=rep(NA, n_rand),HERBRICH_SE=rep(NA, n_rand),TREEEVEN_SE=rep(NA, n_rand),TREERICH_SE=rep(NA, n_rand),
                                 nplot=rep(NA, n_rand))
 
 for(i in 1:n_rand){
   print(i)
   if(i/200 == round(i/200)) write.csv(null_coeff_tab_pair, "null_coeff_pres.csv", row.names=F)
   coln<-paste("Ran_Ab",i,sep="")
   hemidat<-parkhemi[which(is.infinite(parkhemi[,coln])==F),]
   nohemidat<-parkhemi[which(is.infinite(parkhemi[,coln])==T),]
   parkpair<-parkhemi[which(parkhemi$pair1%in%hemidat$pair1),]
   parkpair<-parkpair[which(parkpair$pair1%in%nohemidat$pair1),]
   parkpair$HEMIPRES<-0
   parkpair$HEMIPRES[which(is.infinite(parkpair[,coln])==F)]<-1
   mod1<-lm(ALLEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
   mod2<-lm(ALLRICH_LOG~HEMIPRES+pair1, data=parkpair)
   mod3<-lm(HERBEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
   mod4<-lm(HERBRICH_LOG~HEMIPRES+pair1,  data=parkpair)
   mod5<-lm(TREEEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
   mod6<-lm(TREERICH_LOG~HEMIPRES+pair1, data=parkpair)
   null_coeff_tab_pair$ALLEVEN_LOGIT[i]<-coef(mod1)[2]
   null_coeff_tab_pair$ALLRICH_LOG[i]<-coef(mod2)[2]
   null_coeff_tab_pair$HERBEVEN_LOGIT[i]<-coef(mod3)[2]
   null_coeff_tab_pair$HERBRICH_LOG[i]<-coef(mod4)[2]
   null_coeff_tab_pair$TREEEVEN_LOGIT[i]<-coef(mod5)[2]
   null_coeff_tab_pair$TREERICH_LOG[i]<-coef(mod6)[2]
   null_coeff_tab_pair$ALLEVEN_SE[i]<-summary(mod1)$coefficients[2,2]
   null_coeff_tab_pair$ALLRICH_SE[i]<-summary(mod2)$coefficients[2,2]
   null_coeff_tab_pair$HERBEVEN_SE[i]<-summary(mod3)$coefficients[2,2]
   null_coeff_tab_pair$HERBRICH_SE[i]<-summary(mod4)$coefficients[2,2]
   null_coeff_tab_pair$TREEEVEN_SE[i]<-summary(mod5)$coefficients[2,2]
   null_coeff_tab_pair$TREERICH_SE[i]<-summary(mod6)$coefficients[2,2]
   null_coeff_tab_pair$nplot[i]<-nrow(parkpair)
 }
 
 #write.csv(null_coeff_tab_pair, "null_coeff_tab_pair_transform.final_1.csv", row.names=F)
 
write.csv(null_coeff_tab_pair, "PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/null_pres_minuszero.csv")

null_coeff_tab_pair <- read.csv("PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/null_pres_minuszero.csv") 
summary(null_coeff_tab_pair)
 ## FIGURES 
 
 parkhemi <- read.csv("PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/parkhemi.pairs.minuszero.csv");summary(parkhemi)
 
 hemidat<-parkhemi[which(parkhemi$HEMIAB>0),]
 
 nohemidat<-parkhemi[which(parkhemi$HEMIAB==0),]
 parkpair<-parkhemi[which(parkhemi$pair1%in%hemidat$pair1),]
 parkpair<-parkpair[which(parkpair$pair1%in%nohemidat$pair1),]
 parkpair$HEMIPRES<-0
 parkpair$HEMIPRES[which(parkpair$HEMIAB>0)]<-1

library(ggplot2)
 mod1<-lm(ALLEVEN_LOGIT~HEMIPRES+pair1, data=parkpair); summary(mod1)
 sum(null_coeff_tab_pair$ALLEVEN_LOGIT>coefficients(mod1)[2])/nrow(null_coeff_tab_pair) # p = 0.459
 
 alleven<-ggplot(null_coeff_tab_pair, aes(x=ALLEVEN_LOGIT))+
   geom_density( fill = "#00AFBB", color = "#00AFBB", alpha = 0.3)+
   # xlim(c(-0.05, 0.15))+
   geom_vline(aes(xintercept=mean(coefficients(mod1)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   # annotate("text", x=0.030, y=40, label= "p = 0.081") +
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         text = element_text(size = 10),
         panel.grid.major = element_blank(),
         # panel.border = element_blank(),
         # axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   # plot.title = element_text(hjust = 0.5))+
   labs(y = expression("All Plants"),
        x= expression(""),
        subtitle = "p = 0.459"); alleven
 
 
 
 mod2<-lm(ALLRICH_LOG~HEMIPRES+pair1, data=parkpair)
 sum(null_coeff_tab_pair$ALLRICH_LOG>coefficients(mod2)[2])/nrow(null_coeff_tab_pair) # p = .783
 
 allrich<-ggplot(null_coeff_tab_pair, aes(x=ALLRICH_LOG))+
   geom_density( fill = "#00AFBB", color = "#00AFBB", alpha = 0.3)+
   # xlim(c(-0.05, 0.15))+
   geom_vline(aes(xintercept=mean(coefficients(mod2)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   # annotate("text", x=0.030, y=40, label= "p = 0.081") +
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         text = element_text(size = 10),
         panel.grid.major = element_blank(),
         # panel.border = element_blank(),
         # axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   # plot.title = element_text(hjust = 0.5))+
   labs(y = expression(""),
        x= expression(""),
        subtitle = "p = 0.784"); allrich
 
 
 
 mod3<-lm(HERBEVEN_LOGIT~HEMIPRES+pair1, data=parkpair); summary(mod3)
 
 sum(null_coeff_tab_pair$HERBEVEN_LOGIT>coefficients(mod3)[2])/nrow(null_coeff_tab_pair) # p = .442
 
 herbeven<-ggplot(null_coeff_tab_pair, aes(x=HERBEVEN_LOGIT))+
   geom_density( fill = "#FC4E07", color = "#FC4E07", alpha= 0.3)+
   # ylim(c(0,25))+
   # xlim(c(-0.05,0.15))+
   #  annotate("text", x=0.038, y=25, label= "p < 0.001") +
   geom_vline(aes(xintercept=mean(coefficients(mod3)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         text = element_text(size = 10),
         panel.grid.major = element_blank(),
         #   panel.border = element_blank(),
         axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   labs(y = expression("Herbaceous Plants"),
        x = expression(""),
        subtitle = "p = 0.442"); herbeven
 
 
 mod4<-lm(HERBRICH_LOG~HEMIPRES+pair1,data=parkpair) 
 
 sum(null_coeff_tab_pair$HERBRICH_LOG>coefficients(mod4)[2])/nrow(null_coeff_tab_pair) # p = .861
 
 herbrich<-ggplot(null_coeff_tab_pair, aes(x=HERBRICH_LOG))+
   geom_density( fill = "#FC4E07", color = "#FC4E07", alpha= 0.3)+
   # ylim(c(0,25))+
   # xlim(c(-0.05,0.15))+
   #  annotate("text", x=0.038, y=25, label= "p < 0.001") +
   geom_vline(aes(xintercept=mean(coefficients(mod4)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         text = element_text(size = 10),
         panel.grid.major = element_blank(),
         #   panel.border = element_blank(),
         axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   labs(y = expression(""),
        x = expression(""),
        subtitle = "p = 0.861"); herbrich
 
 
 mod5<-lm(TREEEVEN_LOGIT~HEMIPRES+pair1, data=parkpair)
 sum(null_coeff_tab_pair$TREEEVEN_LOGIT>coefficients(mod5)[2])/nrow(null_coeff_tab_pair) # p = 0.112
 
 treeeven<-ggplot(null_coeff_tab_pair, aes(x=TREEEVEN_LOGIT))+
   # xlim(c(-0.10, 0.05))+
   geom_density( fill="#69b3a2", color="#69b3a2", alpha=0.4)+
   geom_vline(aes(xintercept=mean(coefficients(mod5)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   # annotate("text", x=-0.068, y=25,  label= "p  = 0.062 ") +
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         text = element_text(size = 10),
         # panel.border = element_blank(),
         axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   labs(y = expression("Woody Plants"),
        x = expression(""),
        subtitle = "p = 0.112"); treeeven
 
 
 mod6<-lm(TREERICH_LOG~HEMIPRES+pair1, data=parkpair)
 sum(null_coeff_tab_pair$TREERICH_LOG>coefficients(mod6)[2])/nrow(null_coeff_tab_pair) # p = .316
 
 treerich<-ggplot(null_coeff_tab_pair, aes(x=TREERICH_LOG))+
   # xlim(c(-0.10, 0.05))+
   geom_density( fill="#69b3a2", color="#69b3a2", alpha=0.4)+
   geom_vline(aes(xintercept=mean(coefficients(mod6)[2])),
              color="black", linetype="dashed", size=1)+
   geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
   # annotate("text", x=-0.068, y=25,  label= "p  = 0.062 ") +
   theme_bw()+
   theme(panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         text = element_text(size = 10),
         # panel.border = element_blank(),
         axis.line.x = element_blank(),
         axis.line.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank())+
   labs(y = expression(""),
        x = expression(""),
        subtitle = "p = 0.316"); treerich
 
 ##
 
 library(ggpubr)
 
 even <- ggarrange(alleven, herbeven, treeeven, nrow = 3, ncol = 1)
 even <- annotate_figure(even, top = text_grob("Evenness")); even
 
 rich <- ggarrange(allrich, herbrich, treerich, nrow = 3 , ncol = 1)
 rich <- annotate_figure(rich, top = text_grob("Richness")); rich
 
 graph.null <- ggarrange(even, rich); graph.null
 graph.null <-  annotate_figure(graph.null,    bottom = text_grob("Model Coefficients"), top = text_grob("Presence", size = 14, vjust = 0.5)); graph.null
 
 graph.null + ggsave("Figures/paper/presencenull.tiff", height = 7, width = 7, dpi = 300)
 graph.null
 ##
 
 
