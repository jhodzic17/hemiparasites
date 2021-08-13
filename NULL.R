###analysis of relationship between hemiparasite abundance / presence and the evenness of plant communities in the US National Parks using NPS inventory and monitoring data.

###R libraries and functions
library(vegan)
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

# setwd("C:/Users/ipearse/Desktop/hemiparasites/analysis")
# write.csv(parkall, "allNPSdata.csv",row.names=F)


##set working directory and load data

#setwd("C:/Users/ipearse/Desktop/hemiparasites/analysis")
#park<-read.csv("NPSdata_clean.csv")
#park<-read.csv("park.eco.csv")
#hemi<-read.csv("hemiparasite_USDAplantcodes.csv")

#park <- read.csv("otherdatasheets/park.eco.csv", header = TRUE); ## subset park (as in for create plot level dataframe + ecoregion)
#hemi<-read.csv("otherdatasheets/hemiparasite_USDAplantcodes.csv")

## merge with Eve area
setwd("~/NPS")
area <- read.csv("plots_areasampled.csv", header = TRUE)
park <- read.csv("otherdatasheets/park.eco.csv", header = TRUE)
hemi<-read.csv("otherdatasheets/hemiparasite_USDAplantcodes.csv")
area <- area %>% dplyr::select(SITEPLOT, Area_sampled)
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



parkhemi <- read.csv("parkhemirandom.csv", header = TRUE)

###Add to reduced dataframe using random subsets of data to describe the presence and abundance of non-hemiparasites
n_rand<-1000
for(j in 1:n_rand){
	print(j)
	rass<-sample(spoth, n_hemi)  #grab random sample of non-hemiparasites
	coln<-rep(NA,nrow(parkhemi))
	parkhemi<-cbind(parkhemi, coln)
	colnames(parkhemi)[(29+j)]<-paste('Ran_Ab',j, sep="")
	parktr<-park
	parktr$Pct_Cov[which(!(parktr$Species%in%names(rass)))]<-0
	rass_cov<-tapply(parktr$Pct_Cov, parktr$SITEPLOT, sum)
	parkhemi[,(29+j)]<-log(rass_cov[parkhemi$SITEPLOT])
	#parkhemi[,(29+j)][which(is.finite(parkhemi[,(29+j)])==F)]<-0
}



###abundance, evenness, all plants ##############################################################

#parkhemi$HEMIAB.log[is.finite(parkhemi$HEMIAB.log) == F] <-0

## abundance, hemiparasite ALL Plants
for(i in 1:n_rand){
	#print(i)
	parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
	modform<-formula(paste("ALLEVEN_LOGIT~ Area_sampled +",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR)"))
	aaa<-lmer(modform, data=parkhemi2)
	nullcoef_tab$COEFF_ALLEVEN_AB[i]<-fixef(aaa)[3]
	#nullcoef_tab$WALDT_ALLEVEN_AB[i]<-summary(aaa)$coefficients[2,3]
}
parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(ALLEVEN_LOGIT~ Area_sampled + HEMIAB.log + (1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
summary(aaa)
sum(nullcoef_tab$COEFF_ALLEVEN_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P = 0.108
#sum(nullcoef_tab$WALDT_ALLEVEN_AB>summary(aaa)$coefficients[2,3])/nrow(nullcoef_tab); # p = 0.76

###abundance, richness, all plants
for(i in 1:n_rand){
	#print(i)
	parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
	modform<-formula(paste("ALLRICH_LOG~ Area_sampled + ",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR) "))
	aaa<-lmer(modform, data=parkhemi2)
	nullcoef_tab$COEFF_ALLRICH_AB[i]<-fixef(aaa)[3]
#	nullcoef_tab$WALDT_ALLRICH_AB[i]<-summary(aaa)$coefficients[2,3]
}

parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(ALLRICH_LOG~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2) 
summary(aaa)
sum(nullcoef_tab$COEFF_ALLRICH_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P= 0.516


###abundance, evenness, herbaceaous plants
for(i in 1:n_rand){
	#print(i)
	parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
	modform<-formula(paste("HERBEVEN_LOGIT ~ Area_sampled + ",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR)" ))
	aaa<-lmer(modform, data=parkhemi2)
	nullcoef_tab$COEFF_HERBEVEN_AB[i]<-fixef(aaa)[3]
#	nullcoef_tab$WALDT_HERBEVEN_AB[i]<-summary(aaa)$coefficients[2,3]
}

parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(HERBEVEN_LOGIT~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
summary(aaa)


sum(nullcoef_tab$COEFF_HERBEVEN_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P=0


###abundance, richness, herbaceous plants
for(i in 1:n_rand){
	#print(i)
	parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
	modform<-formula(paste("HERBRICH_LOG~ Area_sampled +",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR)"))
	aaa<-lmer(modform, data=parkhemi2)
	nullcoef_tab$COEFF_HERBRICH_AB[i]<-fixef(aaa)[3]
#	nullcoef_tab$WALDT_HERBRICH_AB[i]<-summary(aaa)$coefficients[2,3]
}
parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(HERBRICH_LOG~ Area_sampled + HEMIAB.log +(1|PARK)+(1|ECO4)+(1|YEAR) , data=parkhemi2)
summary(aaa)
nullcoef_tab <- nullcoef_tab[!is.na(nullcoef_tab$COEFF_HERBRICH_AB),]
sum(nullcoef_tab$COEFF_HERBRICH_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P=0 .584
#sum(nullcoef_tab$WALDT_HERBRICH_AB>summary(aaa)$coefficients[2,3])/nrow(nullcoef_tab) #P=0.2522

## trees
###abundance even

for(i in 1:n_rand){
  #print(i)
  parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
  modform<-formula(paste("TREEEVEN_LOGIT ~ Area_sampled +",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR)"))
  aaa<-lmer(modform, data=parkhemi2)
  nullcoef_tab$COEFF_TREEEVEN_AB[i]<-fixef(aaa)[3]
  #	nullcoef_tab$WALDT_HERBEVEN_AB[i]<-summary(aaa)$coefficients[2,3]
}

parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(TREEEVEN_LOGIT~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR) , data=parkhemi2)
summary(aaa)
sum(nullcoef_tab$COEFF_TREEEVEN_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P=.92



###abundance, richness, tree
for(i in 1:n_rand){
  #print(i)
  parkhemi2<-parkhemi[which(is.finite(parkhemi[,(29+i)])==T),]
  modform<-formula(paste("TREERICH_LOG~ Area_sampled +",colnames(parkhemi)[(29+i)],"+(1|PARK)+(1|ECO4)+(1|YEAR)"))
  aaa<-lmer(modform, data=parkhemi2)
  nullcoef_tab$COEFF_TREERICH_AB[i]<-fixef(aaa)[3]
  #	nullcoef_tab$WALDT_HERBRICH_AB[i]<-summary(aaa)$coefficients[2,3]
}
parkhemi2<-parkhemi[which(is.finite(parkhemi$HEMIAB.log)==T),]
aaa<-lmer(TREERICH_LOG ~ Area_sampled + HEMIAB.log +(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
summary(aaa)
nullcoef_tab <- nullcoef_tab[!is.na(nullcoef_tab$COEFF_TREERICH_AB),]
sum(nullcoef_tab$COEFF_TREERICH_AB>fixef(aaa)[3])/nrow(nullcoef_tab) #P=0 .06



write.csv(nullcoef_tab, "PARKHEMIFINAL/SCRIPT/post edit/random/null model results.csv")



## FIGURES 

nullcoef_tab<-read.csv("PARKHEMIFINAL/SCRIPT/post edit/random/null model results.csv", header = TRUE)

parkhemi <- read.csv("parkhemirandom.csv")
parkhemi2 <-parkhemi[is.finite(parkhemi$HEMIAB.log)==T,]

all.mod<-lmer(ALLEVEN_LOGIT~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
herb.mod <- lmer(HERBEVEN_LOGIT~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR),  data=parkhemi2)
tree.mod <- lmer(TREEEVEN_LOGIT ~ Area_sampled + HEMIAB.log +(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
all.richmod <-lmer(ALLRICH_LOG ~ Area_sampled+ HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR), data=parkhemi2)
herb.richmod <- lmer(HERBRICH_LOG ~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR) ,data=parkhemi2)
tree.richmod <-lmer(TREERICH_LOG ~ Area_sampled + HEMIAB.log+(1|PARK)+(1|ECO4)+(1|YEAR) ,data=parkhemi2)


head(nullcoef_tab)


alleven<-ggplot(nullcoef_tab, aes(x=COEFF_ALLEVEN_AB))+
  geom_density( fill = "#00AFBB", color = "#00AFBB", alpha = 0.3)+
  xlim(c(-0.05, 0.15))+
  geom_vline(aes(xintercept=mean(fixef(all.mod)[3])),
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
       subtitle = "p = 0.11"); alleven

herbeven<-ggplot(nullcoef_tab, aes(x=COEFF_HERBEVEN_AB))+
  geom_density( fill = "#FC4E07", color = "#FC4E07", alpha= 0.3)+
  ylim(c(0,25))+
  xlim(c(-0.05,0.15))+
  #  annotate("text", x=0.038, y=25, label= "p < 0.001") +
  geom_vline(aes(xintercept= 0.15),  color="black", linetype="dashed", size=1)+ ## not fixef(model) becuase that = 0.153 and to keep 0's aligned in graphs, limit needs to be 0.15
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
       subtitle = "p <0.001"); herbeven

treeeven<-ggplot(nullcoef_tab, aes(x=COEFF_TREEEVEN_AB))+
  geom_density( fill="#69b3a2", color="#69b3a2", alpha=0.4)+
  ylim(c(0,25))+
  xlim(c(-0.05,0.15))+
  #  annotate("text", x=0.038, y=25, label= "p =.986")+
  geom_vline(aes(xintercept=mean(fixef(tree.mod)[3])),
             color="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        #  panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y = expression("Woody Plants"),
       x = expression(""),
       subtitle = "p = 0.92"); treeeven

allrich<-ggplot(nullcoef_tab, aes(x=COEFF_ALLRICH_AB))+
  geom_density( fill = "#00AFBB", color = "#00AFBB", alpha = 0.3)+
  xlim(c(-0.10, 0.05))+
  geom_vline(aes(xintercept=mean(fixef(all.richmod)[3])),
             color="black", linetype="dashed", size=1)+ 
  geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
  # annotate("text", x= - 0.055, y=40, label= "p = .465" ) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        #  panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  # plot.subtitle = element_text(hjust = 0.5))+
  labs(y = expression(""),
       x = expression(""),
       subtitle = "p = 0.52"); allrich

herbrich<-ggplot(nullcoef_tab, aes(x=COEFF_HERBRICH_AB))+
  geom_density( fill = "#FC4E07", color = "#FC4E07", alpha= 0.3)+
  xlim(c(-0.10, 0.05))+
  geom_vline(aes(xintercept=mean(fixef(herb.richmod)[3])),
             color="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept = 0 ), color = "black", linetype = "solid", size = 0.5)+
  # annotate("text", x=-0.068, y=25,  label= "p  = .584") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 10),
        #panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y = expression(""),
       x = expression(""),
       subtitle = "p = 0.58"); herbrich


treerich<-ggplot(nullcoef_tab, aes(x=COEFF_TREERICH_AB))+
  xlim(c(-0.10, 0.05))+
  geom_density( fill="#69b3a2", color="#69b3a2", alpha=0.4)+
  geom_vline(aes(xintercept=mean(fixef(tree.richmod)[3])),
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
       subtitle = "p = 0.32"); treerich


even <- ggarrange(alleven, herbeven, treeeven, nrow = 3, ncol = 1)
even <- annotate_figure(even, top = text_grob("Evenness")); even

rich <- ggarrange(allrich, herbrich, treerich, nrow = 3 , ncol = 1)
rich <- annotate_figure(rich, top = text_grob("Richness")); rich

graph.null <- ggarrange(even, rich); graph.null
graph.null <-  annotate_figure(graph.null,    bottom = text_grob("Model Coefficients"))

graph.null + ggsave("Figures/paper/null.tiff", height = 7, width = 7, dpi = 300)
