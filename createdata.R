##Ian Pearse and Jasna Hodzic


## Creating plot-level data from species level sheet

###R libraries and functions
library(vegan)
library(lme4)
library(magrittr)
library(dplyr)

###IP - combine data from each National Park into one big file with all I&M data from NPS


setwd("~/NPS")
park<-read.csv("otherdatasheets/allNPSdata.csv") 
hemi<-read.csv("otherdatasheets/hemiparasite_USDAplantcodes.csv")

n_hemi<-sum(hemi$Accepted.Symbol%in%names(spdat_ab)); n_hemi
## of the list of species that had NA's , I sorted them into herbaceous, grass, nonvascular, aquativ, vine, and tree/shrub: 

herbnam <- c(as.character( c("ASTERA", "AGOSE2", "AGERA2", "AGOSE", "AMARA", "AGROS2" ,"AMBRO" ,"AGROP2", "ANTEN", 
                             "ASTERA", "ARIST", "AQUIL", "ARABI2" ,"ASCLE", "ASTER", "BAHIA", "BIDEN",
                             "BRASSI", "BESSE", "BOUTE" , "CAPSE" ,"CHAMA15" , "HORDE", "DALEA", 
                             "BROMU",  "CASTI2", "DRABA", "CIRSI", "CLEMA", "CHAMA1", "CHENO", "HELIO4",
                             "CORYP", "CYPER", "ELEOC", "EUPHO", "HACKE", "DALEA", "DESCU", "ERYSI",
                             "LAENN", "LEMNA", "MACHA", "MENTZ", "MERTE", "MIRAB", "PACKE", "PECTI",
                             "PHYSA2", "PLANT", "SCHOE6", "SISYM", "TRAGO", "TETRA17", "THELE", 
                             "HELIA3", "LEPID", "MONAR", "MUHLE", "OENOT", "OXALI", "PHLEU", "PHYSA",
                             "CRYPT", "DANTH", "POLYG4", "PORTU", "SOLAN", "SPORO", "SYMPH4", "HYMEN7",
                             "DELPH" ,"EUPHO", "ECHIN3", "ELYMU", "EPILO", "EQUIS" , "ERIGE2", "FESTU",
                             "FRAGA", "GALIU", "GERAN", "HEUCH", "IPOMO2", "LACTU", "LATHY", "LOTUS",
                             "LUPIN", "MEDIC", "MUHLE", "LITHO3", "LOTUS", "PANIC", "PENST", "PHACE", "PHLOX", "POTEN",
                             "RANUN", "ROSA5", "RUBUS", "RUMEX", "SENEC", "PEDIC", "SOLID", "STIPA", "TRIFO",
                             "VICIA", "VIOLA", "VALER", "GILIA", "UNKFORB" , "ERIOG", "CLARK", "SILEN",
                             "ARNIC", "SAXIF", "ZIGAD", "CERAS", "CERAS", "BOARGI", "CLAYT", "STREP", "THELY",
                             "CAMIS", "THALI2", "LESQU", "LINAN", "HIERA", "SILEN", "ALLIU", "ANEMO", "ORTHO4",
                             "CLARK", "PYROL", "SCHIS", "ACHNA", "OXYTR", "CUSCU", "ARENA", "PRENA", "OSMOR",
                             "STEPH", "ARNIC", "MITER", "GENTI", "SPHAG2", "STELL", "GUTIE", "LOMAT", "SAXIF", "AMSIN",
                             "ZIGAD", "AMSIN", "CHAEN", "CREPI", "KRAME", "LINUM", "PLAGI", "PUPAP", "THERM",
                             "BORAGI", "CLAYT","PECTO", "VERON", "FABACE", "FORB", "POLEM" , "TARAX", "BOERH2", "IRIS",
                             "MADIA", "PLATA2", "ANDRO3", "ANGEL", "AUREO", "CORAL", "ESCHS", "HABEN", "HEDEO", "LIATR",
                             "THYSA", "VERAT", "CRYPT3", "FRAXI", "GEUM", "ZIZIA", "ZINNI", "UVULA", "UNFORB" , "TRITE",
                             "TRIOS" , "TOFIE", "THYMO", "THLAS", "DESMO", "DODEC", "GAYOP", "GEUM", "HYPOC", "MIMUL", "MINUA", 
                             "MONAR2", "NAMA4", "SCUTE", "TOWNS", "ACTAE", "APOCY", "BRASS2", "CARDA", "COLLI", "ERITR", "HERAC", 
                             "HETER8", "HYMEN4", "HYPER", "IMPAT", "LESPE", "LOBEL", "MAIAN",  "OROBA", "POLYG", "PSORA2" ,"RORIP",
                             "SEDUM", "SONCHUS", "STREP2", "UNFORB2", "ABRON", "ALYSS", "AMSON", "ASCLE", "ASTRASP1", "ASTRASP2",
                             "BRASS", "BRODI", "CALTH", "CHROI2", "COLLO", "ERODI" , "EUPAT", "EURYB", "FRITI", "GNAPH", "HAPLO11",
                             "HESPE7", "HEXAS", "HYDRO4", "KRIGI", "LANGL", "LAPPU", "LIGUS", "LILIU", "LISTE", "MALVA", "MYOSO",
                             "NAVAR", "NEMOP", "ORTHO", "PERIT", "PETAL", "POBIB4", "PRIMU", "PSATH", "SIDAL", "TOFIE" , "THYMO",
                             "2FORB", "ABUTI", "ACALY" , "AGRIM", "ACLEI", "AGAST", "ALETE", "ALICI", "ALISM", "AMARAN", "AMORP",
                             "APIA" , "APIAXX", "ARALI", "ARCTI", "ARUNC", "ASCLEP", "ASLI14", "ASPER", "ASPLE", "ASTER_sp", "ASTRA",
                             "BORAG", "BERBE", "ORCHID", "ZIZAN", "XYLOR", "LILIUM", "DUDLE", "FABAC", "PLECT" , "VINCA", "VERBE")))
                             ## Including bulbs like Camas and Calochortus, orchid (#VERAT, UVULA, MAIAN POLYG, FRITI,  TOFIELILIU)
## hemis (KRAME, ORTHO, ORTHO4)
## orchids : #HABEN  

grassnam <- c(as.character(c("MELIC", "DICHA2", "POACEA" , "UNGRAS", "CALAM", "LUZUL", "GRASS", "PO3", "POACEA", "POACXX",
                           "UNKGRAM", "NASSE", "AVENA", "XYRIS", "UNGRASS", "HESPE11",  "LOLIU" , "PARNA", "TRISE",
                           "VULPI", "POA", "TRISE", "DESCH", "ENNEA", "ERIOP", "GLYCE", "LEYMU", "POA1", "POA2", "CAREX" ,".grass",
                           "ACORU", "AGVIV", "ALOPE", "ANDRO2", "ARGEM", "ARISA", "ARNOG", "CAREX", "CAREX3", "CLINT",
                           "Legum", "GRASS3", "POACE")))
              
fernnam <- c(as.character(c("ATHYR", "BOTRY", "WOODS", "POLYP", "WOODS")))

hemip <- c(as.character("AGALI" , "CASTI")) ## HEMIPARASITES.. ## CHECK CASTI

woodynam <- c(as.character( c("AMELA", "ABIES", "ACER", "ARTEM", "BRICK","EPHED", "EPTO*", "OPUNT","QUERC",  "ROBIN", "SCLER10",
                             "RIBES", "SPHAE", "SALIX", "SAMBU", "CRATA", "QURU", "CALOC", "POTR5S", "PURSH", "SPIRA", 
                             "CARYA", "SANIC", "SYMPH", "RHODO", "CORDY", "CHRYS9", "GARRY", "MAMMI", "CASTA", "PSORO",
                             "UNKSHB", "ULMUS", "CEANO", "VACCI", "ACARO2", "ARCEU", "ENCEL", "LYCIU", "PHILA", "PHORA",
                             "PINUS", "POBA2S", "POLYS", "PRUNU", "RAMAL2", "YUCCA" , "ATRIP", "BACCH", "BERNA", "NEEA",
                             "TETRA3", "2SHRUB" , "AGAVE", "ALNUS", "ARONI", "BETUL", "GYMNO", "PICEA", "TSUGA", "ERICA2"))) ##including cacti ("MAMMI"")

nonvasnam <- c(as.character (c("MOSS", "SELAG", "LICHEN", "UNKMOSP", "USNEA2", "PELTI2", "PERID", "RHIZO3", "BRYO",
                               "CLADI3", "CLADO3", "UNKMOSP", "POLYT5", "DICRA8", "CRUST", "XANTH7", "LEUCO9", "AULAC2",
                             "CALOP7", "LICHEN3", "BRACH10", "BRYUM2", "TORTU", "THUD",  "GRIMM2",  "LETHA", "MNIUM2",
                             "CANDE2", "FISSI2", "HEDWI", "HELIA", "HELIA2", "nonvasc", "PLEOP2", "ANOMO", "ARCTO", "ARCTO3",
                             "ASPIC2-dark", "BAZZA", "ASPIC2-light", "BRYOR2", "MARCH", "MOSS2", "THUID" , "XANTH",
                               "XANTH2", "ALGAE", "ATHYR", "ATRIC2", "LICHE", "LICHE5", "moss", "MOSS3", ".moss", ".liverwort", ".liver" )))

aquaticnam <- c(as.character(c("JUNCU")))


vinenam <- c(as.character(c("VITIS", "CALYS", "DIOSC", "FUNAS", "LONIC", "SMILA2")))

holoparnam <- c(as.character( c("CUSCU", "ARCEU" , "CALYS" , "OROBA", "PHORA"))) ## Jasna - also put these under whatever USDA said

mycoparnam <- c(as.character(c("CORAL"))) ## corahlizza. Also included in forbs


## assigning them a value in growth habit:
park$Growth_habit <- as.factor(park$Growth_habit)
levels(park$Growth_habit)

park[park$Species %in% herbnam ,  "Growth_habit"] <- "Forb/herb"
park[park$Species %in% grassnam, "Growth_habit"] <- "Graminoid"
park[park$Species %in% woodynam , "Growth_habit"] <- "Shrub"
park[park$Species %in% nonvasnam, "Growth_habit"] <- "Nonvascular"
park[park$Species %in% aquaticnam, "Growth_habit"] <- "Nonvascular" ## Just so I can take itout

park[park$Species %in% vinenam, "Growth_habit"] <- "Vine"
park[park$Species %in% hemip, "Growth_habit"] <-"Forb/herb"
  
  


## checking: 
park %>% filter(Species == "POA") %>% select(Species, Growth_habit)

## did not do anything wit hthe ferns or aquatic species. Fine for them to be NA and excluded, anyways.

## fix the year
park[is.na(park$Year), "Year"] <- "Unknown"
park$Year <-as.factor(park$Year)
levels(park$Year)
## let's make every hemiparsite go to Forb/herb, so we avoid that annoying NOHEMIEVEN stuff

hnam<-c(as.character(hemi$Accepted.Symbol), as.character(hemip))
park[park$Species %in% hnam, "Growth_habit"] <- "Forb/herb"



## cleaning it up
park <- park %>% filter(!Species %in% c("other" , ".liche", ".moss."), !is.na(Species), !Stratum %in% c("a1", "A2", "A3", "FA", "H3"),
                        !Growth_habit %in% c("Nonvascular" , "Lichenous"), !is.na(Lat), !is.na(Long), 
                        !Species %in% c("UNKSP",  "U
                                        NKSP."))

#park.new <- merge(park, genus, by = "Species", all = TRUE)
park.nohemi <- park %>% filter(!Species %in% hnam)
park.nohemi$SITEPLOT <- paste(park.nohemi$Site, park.nohemi$Plot)
park.nohemi$Genus <- "NA"
park.ab <- park %>% filter(Species %in% hnam)
genus <- hemi %>% dplyr::select(Genus, Accepted.Symbol) %>% rename(Species = Accepted.Symbol)

park.ab <- merge(park.ab, genus, by = "Species")
park.ab <- park.ab %>% distinct() ## you are adding some damncolumnshee
park.ab$SITEPLOT <- paste(park.ab$Site, park.ab$Plot)



park<-rbind(park.ab, park.nohemi)

##abundance and frequency






### 02.16.2022 addition 
parkhemi.original <- read.csv("parkhemirandom.csv"); 
parkhemi <- parkhemi.original
parkhemi <- parkhemi %>% filter(Area_sampled > 398 & Area_sampled < 401) ## allowing these plots in
park <- read.csv("otherdatasheets/park.eco.csv")
park <- park %>%
  filter(SITEPLOT %in% parkhemi$SITEPLOT) ## restricting plot

park$Species <- as.factor(park$Species); nlevels(park$Species)
nlevels(park$Species)
park.freq <- park %>% group_by(Species) %>% summarize(freq = length(Species),
                                                      ab = sum(Pct_Cov)) %>% arrange(desc(freq));

park.freq <- park.freq %>%
  mutate(rank.freq = rownames(park.freq)); head(park.freq)

write.csv(park.freq, "PARKHEMIFINAL/manuscript_rep/most_Recent/revisions/allab.csv")

park.ab1 <- park.freq %>% arrange(desc(ab));
park.ab1 <- park.ab1 %>% mutate(rank.ab = rownames(park.ab1)); head(park.ab1)

hemi.freq <- park.ab1 %>% filter(Species %in% hnam)
head(hemi.freq)
## plotting
ggplot(hemi.freq, aes(x= freq, y = ab)) + geom_point() + theme_bw()

freqhead(park)
park$SITEPLOT<-paste(park$Site, park$Plot);park$SITEPLOT

## genus #
park.hemi <- park[park$Species %in% hnam,];head(park.hemi)

genus <- hemi %>% dplyr::select(Genus, Accepted.Symbol) %>% rename(Species = Accepted.Symbol)

## mergeing

hemi.freq. <- merge(hemi.freq,genus) %>% arrange(desc(ab)); head(hemi.freq.) ## adding in the Genus
write.csv(hemi.freq., "otherdatasheets/hemi.freq.new.csv")
####
park.genus <- merge(park.hemi, genus)
park.genus$Genus <- as.factor(park.genus$Genus)
park.genus <- park.genus %>% distinct()
total.ab <- as.numeric(sum(park.genus$Pct_Cov))
park.genus.stat <- park.genus %>% select(Pct_Cov, Genus) %>% group_by(Genus) %>% summarize(tot = sum(Pct_Cov),
                                                                                           percentage = (sum(Pct_Cov)/total.ab) * 100) %>%
  arrange(desc(percentage)); park.genus.stat


## how many distinct hemiparasites
hemi.name <- park[park$Species %in% hnam,]
hemi.name$Species <- as.factor(hemi.name$Species)
hemi.name$Species <- droplevels(hemi.name$Species)
nlevels(hemi.name$Species) ## we are losing some due to lat/long
hemi.list <- as.character(levels(hemi.name$Species))
##
allhemi.stat <- park[park$Species %in% hemi.list,]
allhemi.stat<-allhemi.stat %>% select(Pct_Cov,Species) %>% group_by(Species) %>% summarize(tot = sum(Pct_Cov)) 
allhemi.stat
write.csv(allhemi.stat, "hemi_ab_400.csv") ## ABUNDANCE BY SYMBOL


hemhemi.species <- hemi[hemi$Accepted.Symbol %in% hemi.list,] ## Just you my friend. (non in synonymn symbol) how do we have 203...
hemi.species <- hemi.species %>% select(!Synonym.Symbol) %>% distinct()  ## we have some definite repeats of the same symbol for different species,

## but et's figure out the abundaunce
## just top 3 genera

#hemi.species <- hemi.species %>% filter(Genus %in% c("Castilleja", "Pedicularis", "Krameria" , "Comandra")) ## ok so 132
hemi.ab.list <- as.character(levels(hemi.species$Accepted.Symbol))

hemi.ab.list
park.hemi.ab <- park[park$Species %in% hemi.ab.list,]
head(park.hemi.ab)
park.hemi.ab.stat<- park.hemi.ab%>% select(Pct_Cov,Species) %>% group_by(Species) %>% summarize(tot = sum(Pct_Cov)) %>%
  arrange(desc(tot)); park.hemi.ab.stat

write.csv(park.hemi.ab.stat, "hemi_ab_400.csv")

## to know genus orh emi
## all hemiparasites 
hemi.ab.all <-merge()
hemi.ab.data <- merge(hemi.species, park.hemi.ab.stat, by.y = "Species", by.x = "Accepted.Symbol") 
hemi.ab.data %>% arrange(desc(tot)) ## interesting. ## should we do the analysis?

## KRER seems suspicious ,
## so we could do an analysis of the most common of the three biggest genera - Castilleja miniata, KRameria grayi, and PEdicularis bracteosa...

## let's try for fun, CAMI12

park.all <- park %>% select(Pct_Cov, Species)%>% group_by(Species) %>% summarize(tot = sum(Pct_Cov)) %>%
  arrange(desc(tot)); head(park.all)

## figurng out their abundaunce

head(park.hemi)
###make a reduced datafrom where each row is a plot within site (i.e one sampling place within each park)


unsite<-unique(park$SITEPLOT);unsite
##  "all" factors can never be 0 (if they rae it is NA) . everything ele can be (treeab, etc.) 
## stil will make evenness = NA if richness = 1.0?
parkhemi<-data.frame(SITEPLOT=unsite, PARK=rep(NA,length(unsite)), PLOT=rep(NA,length(unsite)),LONG=rep(NA,length(unsite)),LAT=rep(NA,length(unsite)),YEAR=rep(NA,length(unsite)),
                     Castilleja = rep(0, length(unsite)), Pedicularis = rep (0, length(unsite)), Krameria = rep(0, length(unsite)), Comandra = rep(0, length(unsite)),
                     Aureolaria = rep(0, length(unsite)), Cordylanthus = rep(0, length(unsite)), Orthocarpus = rep(0, length(unsite)), Dasistoma = rep(0, length(unsite)),
                     Geocaulon = rep(0, length(unsite)), Triphysaria = rep(0, length(unsite)), Euphrasia = rep(0, length(unsite)), Seymeria = rep(0,length(unsite)),
                     Castilleja.rich = rep(0, length(unsite)), Pedicularis.rich = rep(0, length(unsite)), Krameria.rich = rep(0, length(unsite)), Comandra.rich = rep(0, length(unsite)),
                     Aureolaria.rich = rep(0, length(unsite)), Cordy.rich = rep(0, length(unsite)), Orthocarpus.rich = rep(0, length(unsite)), Dasi.rich = rep(0, length(unsite)),
                     Geo.rich = rep(0, length(unsite)), Tri.rich = rep(0, length(unsite)), Eu.rich = rep(0, length(unsite)), Sey.rich = rep(0,length(unsite)),##all the useful info from 'park' ## JH - why are some NA
                     ALLRICH=rep(NA,length(unsite)),ALLEVEN=rep(NA,length(unsite)), ALLSHANON=rep(NA,length(unsite)),ALLAB=rep(0,length(unsite)),ALLINVAB=rep(0,length(unsite)),ALLNATRICH=rep(NA,length(unsite)),ALLNATAB=rep(NA,length(unsite)),
                     HERBRICH=rep(0,length(unsite)),HERBEVEN=rep(0,length(unsite)),HERBSHANON=rep(0,length(unsite)),HERBAB=rep(0,length(unsite)),HERBNATRICH=rep(0,length(unsite)),HERBNATAB=rep(0,length(unsite)),HERBINVAB=rep(0,length(unsite)),
                     TREEAB = rep(0, length(unsite)), TREERICH = rep(0, length(unsite)), TREEEVEN = rep(0, length(unsite)), TREESHANON = rep(0, length(unsite)), TREENATRICH = rep(0, length(unsite)), TREENATAB = rep(0, length(unsite)), 
                     TREEINVAB = rep(0, length(unsite)),HEMIAB=rep(0,length(unsite)),HEMIRICH=rep(0,length(unsite)),  NOHEMIRICH = rep(0, length(unsite)), NOHEMIEVEN = rep(0, length(unsite)), NOHEMISH = rep(0, length(unsite)), NOHEMIAB = rep(0,length(unsite)))

## NA's for things I want "NA" of if there is nothing of tht species (i.e. TREEEVEN = NA if there are no trees)
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
 # parkhemi$RATIO[i]<-parki$Ratio[1]
  parkvec<-tapply(parki$Pct_Cov, as.character(parki$Species),sum)
  parkhemi$ALLRICH[i]<-specnumber(parkvec)
  parkhemi$ALLSHANON[i]<-diversity(parkvec,index='shannon')
  parkhemi$ALLEVEN[i]<-parkhemi$ALLSHANON[i]/log(parkhemi$ALLRICH[i])
  parkhemi$ALLAB[i]<-sum(parkvec)
  if("I" %in% parki$Exotic) parkhemi$ALLINVAB[i]<-sum(parki$Pct_Cov[which(parki$Exotic=="I")])
  if("N" %in% parki$Exotic) parkhemi$ALLNATAB[i]<-sum(parki$Pct_Cov[which(parki$Exotic=="N")])
  if("N" %in% parki$Exotic) parkhemi$ALLNATRICH[i]<-length(parki$Pct_Cov[which(parki$Exotic=="N")])
  parkiherb<-parki[which(parki$Growth_habit %in% c("Forb/herb", "Subshrub", "Graminoid" , "Vine")),] ## with Ian's way, we go for Stratum over "Growth_HAbit"
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
    parkhemi$TREESHANON[i]<-diversity(treevec,index='shannon')
    parkhemi$TREEEVEN[i]<-parkhemi$TREESHANON[i]/log(parkhemi$TREERICH[i])
    parkhemi$TREEAB[i]<-sum(treevec)
    if("I" %in% parkitree$Exotic) parkhemi$TREEINVAB[i]<-sum(parkitree$Pct_Cov[which(parkitree$Exotic=="I")])
    if("N" %in% parkitree$Exotic) parkhemi$TREENATAB[i]<-sum(parkitree$Pct_Cov[which(parkitree$Exotic=="N")])
    if("N" %in% parkitree$Exotic) parkhemi$TREENATRICH[i]<-length(parkitree$Pct_Cov[which(parkitree$Exotic=="N")])
  }
  if (nrow(parkiherb)>0){
  phemi<-parki[which(as.character(parki$Species)%in%hnam),]
  if(nrow(phemi)>0){
    hemivec<-tapply(phemi$Pct_Cov, as.character(phemi$Species),sum)
    parkhemi$HEMIAB[i]<-sum(hemivec)
    parkhemi$HEMIRICH[i]<-length(hemivec)
  }
  if(nrow(phemi)>0){
    pCastilleja <- parki[parki$Genus  == "Castilleja",]
    pPedicularis <- parki[parki$Genus == "Pedicularis",]
    pKrameria <- parki[parki$Genus == "Krameria",]
    pComandra <- parki[parki$Genus == "Comandra",]
    pAur <- parki[parki$Genus == "Aureolaria",]
    pCordy <- parki[parki$Genus == "Cordylanthus",]
    pOrtho <- parki[parki$Genus == "Orthocarpus",]
    pDasi <- parki[parki$Genus == "Dasistoma",]
    pGeo <- parki[parki$Genus == "Geocaulon",]
    pTri <- parki[parki$Genus == "Triphysaria",]
    pEu <- parki[parki$Genus == "Euphrasia",]
    pSey <- parki[parki$Genus == "Seymeria",]
    Cvec <- tapply(pCastilleja$Pct_Cov, as.character(pCastilleja$Species), sum)
    Pvec <- tapply(pPedicularis$Pct_Cov, as.character(pPedicularis$Species), sum)
    Kvec <- tapply(pKrameria$Pct_Cov, as.character(pKrameria$Species), sum)
    Covec <- tapply(pComandra$Pct_Cov, as.character(pComandra$Species), sum)
    Auvec <-  tapply(pAur$Pct_Cov, as.character(pAur$Species), sum)
    Cordyvec <- tapply(pCordy$Pct_Cov, as.character(pCordy$Species), sum)
    Orthovec <- tapply(pOrtho$Pct_Cov, as.character(pOrtho$Species), sum)
    Dasivec <-  tapply(pDasi$Pct_Cov, as.character(pDasi$Species), sum)
    Trivec <- tapply(pTri$Pct_Cov, as.character(pTri$Species), sum)
    Geovec <- tapply(pGeo$Pct_Cov, as.character(pGeo$Species), sum)
    Euvec <- tapply(pEu$Pct_Cov, as.character(pEu$Species), sum)
    Seyvec <- tapply(pSey$Pct_Cov, as.character(pSey$Species), sum)
    parkhemi$Castilleja[i] <- sum(Cvec)
    parkhemi$Pedicularis[i] <- sum(Pvec)
    parkhemi$Krameria[i] <- sum(Kvec)
    parkhemi$Comandra[i] <- sum(Covec)
    parkhemi$Aureolaria[i] <-sum(Auvec)
    parkhemi$Dasistoma[i] <-sum(Dasivec)
    parkhemi$Cordylanthus[i] <-sum(Cordyvec)
    parkhemi$Geocaulon[i] <-sum(Geovec)
    parkhemi$Triphysaria[i]<-sum(Trivec)
    parkhemi$Seymeria[i] <-sum(Seyvec)
    parkhemi$Orthocarpus <- sum(Orthovec)
    parkhemi$Euphrasia <-sum(Euvec)
    parkhemi$Castilleja.rich[i] <- length(Cvec)
    parkhemi$Pedicularis.rich[i] <- length(Pvec)
    parkhemi$Krameria.rich[i] <- length(Kvec)
    parkhemi$Comandra.rich[i] <- length(Covec)
    parkhemi$Aureolaria.rich[i] <- length(Auvec)
    parkhemi$Cordy.rich[i] <- length(Cordyvec)
    parkhemi$Orthocarpus.rich[i] <-length(Orthovec)
    parkhemi$Dasi.rich[i] <- length(Dasivec)
    parkhemi$Geo.rich[i] <- length(Geovec)
    parkhemi$Tri.rich[i] <- length(Trivec)
    parkhemi$Eu.rich[i] <-length(Euvec)
    parkhemi$Sey.rich[i] <- length(Seyvec)
    
  }
  pnohemi<-parki[!parki$Species %in% hnam & parki$Growth_habit %in% c("Forb/herb", "Subshrub", "Graminoid"),] ### or just use parkiherb
  if (nrow(pnohemi)>0){
  pnohemivec<-tapply(pnohemi$Pct_Cov, as.character(pnohemi$Species),sum)
  parkhemi$NOHEMIRICH[i]<-specnumber(pnohemivec)
  parkhemi$NOHEMIAB[i]<-sum(pnohemivec)
  parkhemi$NOHEMISH[i]<-diversity(pnohemivec, index = 'shannon')
  parkhemi$NOHEMIEVEN[i]<- parkhemi$NOHEMISH[i]/log(parkhemi$NOHEMIRICH[i])
  }
  }
}


summary(parkhemi)




write.csv(parkhemi, "PARKHEMIFINAL/parkhemi.genus.csv")
parkhemi <- read.csv("PARKHEMIFINAL/parkhemi.genus.csv", header = TRUE)


## merging with ecodata from prior dataframe

data.eco <- read.csv("PARKHEMIFINAL/parkhemi.ecoall.020421.csv", header = TRUE)

eco <- data.eco %>% select(SITEPLOT, L3_KEY, L4_KEY)
parkhemi.new <- merge(eco, parkhemi, by = "SITEPLOT"); head(parkhemi.new)

write.csv(parkhemi.new, "PARKHEMIFINAL/parkhemi.eco.csv")


