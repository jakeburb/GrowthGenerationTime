rm(list = ls())
library(gulf)
library(FSA)
library(ggplot2)
library(car)
library(dplyr)
cat("\f")
clg()

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

fp <- getwd()

start.year <- 1988
current.year <- 2021

#* means scientific equipment
#OTM*: Trawls: Midwater trawls;  - otter trawls
#OTB*: Trawls: Bottom trawl;  - otter trawls
#GNS : Gillnets and entagling nets: Set gillnets (anchored)
#GNS*: Experimental Gillnets and entagling nets: Set gillnets (anchored)
#FPN : Traps: Stationary uncovered pount-nets
#PS  : Surrounding Nets: with purse lines (purse seines)
#GND : Gillnets and entagling nets: Drift nets

# Load herring Oracle data:

chan = oracle.open()

sampling.detail=
  paste0("select i.SAMPLING_INFORMATION_ID,i.NAFO_GEAR_CODE,i.SAMPLE_DATE,i.SAMPLE_NUMBER,i.LAB_CODE,i.SAMPLING_PROTOCOL_CODE,i.NAFO_UNIT_AREA_CODE,i.HERRING_FISHING_AREA_CODE,i.LATITUDE,i.LONGITUDE,i.VESSEL_CODE,i.CARRIER_CODE,i.COLLECTED_BY,i.AGED_BY,i.REMARK,i.PROCESSING_PLANT_CODE,i.SET_NUMBER,i.RESEARCH_CODE,i.FOREIGN_VESSEL_CODE,i.LENGTH_FREQUENCY_DEVICE_CODE,i.MESH_SIZE,i.PORT_CODE,i.PORT_NAME,i.COORDINATE_METHOD_CODE,i.MISSION,i.LANDED_WEIGHT,i.DETAIL_DEVICE_CODE,i.WHARF_CODE,i.GEAR_CODE,i.NUMBER_OF_FISH_MEASURED,i.COLLECTED_BY_ID,i.AGED_BY_ID,i.GEAR_FISH_DURATION_HOUR,i.FAT_PCT_METHOD_CODE,i.FAT_METER_CALIB_LOW_READING,i.FAT_METER_CALIB_HIGH_READING,i.PROCESSED_DATE,i.NUMBER_OF_FISH_SAVED,i.COORDINATE_EPSG_ID, d.*
from glf_herring.v_herring_sampling_information i, glf_herring.v_herring_detail d
where i.sampling_information_id=d.sampling_information_id
and i.nafo_gear_code in ('OTM*', 'OTB*', 'OTB', 'GNS', 'GNS*', 'FPN', 'PS','GND')
and extract(year from sample_date) BETWEEN ", start.year, " AND ", current.year, "
order by i.sample_date")


q2 <- oracle.query(channel = chan, query = sampling.detail)

q2$year <- as.numeric(substr(q2$SAMPLE_DATE, 1, 4))
q2$month <- as.numeric(substr(q2$SAMPLE_DATE, 6, 7))
q2$day <- as.numeric(substr(q2$SAMPLE_DATE, 9, 10))
q2 = q2[,c(1,2,4,5,7,9,10,11,21,22,23,26,29,30,38,40,41,42,43,44,45,46,47,48,50,51,52,53,54,61,62,63)]

head(q2)


##Clean data
det<-subset(q2, FISH_AGE != "99")  ##remove all fish that couldnt be aged
str(det)
det$year<-as.character(det$year)

#################SPRING###########
##subset for spring
detS<-subset(det, GSI_SPAWNING_SEASON == "P")
detS<-subset(detS, FROZEN_FISH_LENGTH <= 500)
head(detS)
str(detS)
detS$length<-detS$FRESH_FISH_LENGTH
detS$age<-detS$FISH_AGE

##Fit VB
vb <- vbFuns()
predict2 <- function(x) predict(x,data.frame(age=ages))

##Get age range for each year
agesum<-group_by(detS, year) %>%
  summarize(minage=min(age), maxage=max(age))


#Specify year factor levels and length
#years<-unique(as.factor(detS$year))
years<-c("1988", "1989", "1990", "1991", "1992", "1993", "1994"
         , "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011"
         , "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021")
nyears<-length(years)


cfs <- cis <- preds1 <- preds2 <- NULL

for (i in 1:nyears) {
  ## Loop notification (for peace of mind)
  cat(years[i],"Loop\n")
  ## Isolate years's data
  tmp1 <- filter(detS,year==years[i])
  ## Fit von B to that year
  sv1 <- vbStarts(length~age,data=tmp1)    #Use 'detS' to have a common start value, and 'tmp1' for year specific start value for each VB model
  fit1 <- nls(length~vb(age,Linf,K,t0),data=tmp1,start=sv1)
  ## Extract and store parameter estimates and CIs
  cfs <- rbind(cfs,coef(fit1))
  boot1 <- Boot(fit1)
  tmp2 <-  confint(boot1)
  cis <- rbind(cis,c(tmp2["Linf",],tmp2["K",],tmp2["t0",]))
  ## Predict mean lengths-at-age with CIs
  ##   preds1 -> across all ages
  ##   preds2 -> across observed ages only
  ages <- seq(-1,15,0.2)
  boot2 <- Boot(fit1,f=predict2)
  tmp2 <- data.frame(year=years[i],age=ages,
                     predict(fit1,data.frame(age=ages)),
                     confint(boot2))
  preds1 <- rbind(preds1,tmp2)
  tmp2 <- filter(tmp2,age>=agesum$minage[i],age<=agesum$maxage[i])
  preds2 <- rbind(preds2,tmp2)
}

##Name rows and check VB parameters and predictions
rownames(cfs) <- rownames(cis) <- years
colnames(cis) <- paste(rep(c("Linf","K","t0"),each=2),
                       rep(c("LCI","UCI"),times=2),sep=".")
colnames(preds1) <- colnames(preds2) <- c("year","age","fit","LCI","UCI")
head(preds1) #preditions across all ages
head(preds2) #predictions across observed ages
head(cis) #VB parameters

###Merge year, VB paramters and VB parameter CIs into one dataframe
cisS<-data.frame(years,cfs, cis)
excel(cisS)
#Max and Min Linf over time series
cisS[which.max(cisS$Linf),]
cisS[which.min(cisS$Linf),]


##QUickly plot Linf over time series
plot(cisS$years, cisS$Linf)
lines(cisS$years, cisS$Linf)

#plot all VB over timeseries
ggplot() +
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI),alpha=0.25) +
  geom_point(data=detS,aes(x=age,y=length),
             size=2,alpha=0.3) +
  geom_line(data=preds1,aes(x=age,y=fit),
            size=1,linetype="dashed") +
  geom_line(data=preds2,aes(x=age,y=fit),
            size=1) +
  scale_y_continuous(name="Total Length (mm)",limits=c(0,450)) +
  scale_x_continuous(name="Age (years)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14, colour= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.position="none") +
  facet_wrap(vars(year))
ggsave("Spring VB CurvesR.png", dpi=600, width=30, height= 30, units= "cm")

###Plot all VB on same graph to show changes over years
ggplot() +
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI,fill=year),alpha=0.25) +
  geom_point(data=detS,aes(x=age,y=length,color=year),
             size=2,alpha=0.3) +
  geom_line(data=preds1,aes(x=age,y=fit,color=year),
            size=1,linetype="dashed") +
  geom_line(data=preds2,aes(x=age,y=fit,color=year),
            size=1) +
  scale_y_continuous(name="Total Length (mm)",limits=c(0,450)) +
  scale_x_continuous(name="Age (years)",breaks=0:15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Spring VB Curves 1PlotR.png", dpi=600, width=30, height= 30, units= "cm")

###Plot Linf over time
cisS$years<-as.numeric(cisS$years)
ggplot(data=cisS,aes(x=years,y=Linf, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=Linf.LCI,ymax=Linf.UCI),alpha=0.25) +
  scale_y_continuous(name= "Asymptotic Length (mm)", breaks = seq(290, 360, by = 10)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

ggsave("Spring LinfR.png", dpi=600, width=30, height= 30, units= "cm")


#######Pull M estimates from Stock assessment model#################
####Use Abundance weighted mortality computed by Francois
source("mseRtools.r")
h<-lisread("smq.rep")
#get biomass at age maximum likelihood estimates
baa = h$Bta
#baa

#get natural mortality
m1mle = h$Mt[1,]
#m1
m2mle = h$Mt[2,]
#m2

mctmp <- read.table("mcoutM1.dat", sep = "")
m1 <- apply( mctmp, 2, median)
m1lci <- apply( mctmp, 2, quantile, 0.025)
m1uci <- apply( mctmp, 2, quantile, 0.975)

mctmp <- read.table("mcoutM2.dat", sep = "")
m2 <- apply( mctmp, 2, median)
m2lci <- apply( mctmp, 2, quantile, 0.025)
m2uci <- apply( mctmp, 2, quantile, 0.975)


#build df
m = baa #just to get proper sized object
m[,1:5] = m1
m[,6:10] = m2
#m

muci = baa #just to get proper sized object
muci[,1:5] = m1uci
muci[,6:10] = m2uci

mlci = baa #just to get proper sized object
mlci[,1:5] = m1lci
mlci[,6:10] = m2lci

#do abundance weighted average
mprop = m*(baa/rowSums(baa))
Mtot = rowSums(mprop)

mprop = muci*(baa/rowSums(baa))
MtotUCI = rowSums(mprop)

mprop = mlci*(baa/rowSums(baa))
MtotLCI = rowSums(mprop)

yr = 1978:2021
plot(yr, Mtot, ylim = c(0.2,0.6), type = "l")
lines(yr, MtotUCI, col = "red")
lines(yr, MtotLCI, col = "red")


Myears<-c("1978", "1979", "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989", "1990", "1991", "1992", "1993", "1994"
          , "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011"
          , "2012", "2013", "2014", "2015", "2016", "2017", "2018","2019", "2020", "2021")
SpringM<-data.frame(Myears,Mtot, MtotUCI, MtotLCI)
SpringM<-subset(SpringM, Myears >= "1988")
###Compute Lopt
LoptS<-cisS$Linf* (3/(3+(SpringM$Mtot/cisS$K)))
LoptUCIS<-cisS$Linf.UCI* (3/(3+(SpringM$MtotUCI/cisS$K.UCI)))
LoptLCIS<-cisS$Linf.LCI* (3/(3+(SpringM$MtotLCI/cisS$K.LCI)))
##Compute Generation Time with estimated T0 from VB models
#GTS<- -cisS$t0 - ((log(1-LoptS/cisS$Linf))/cisS$K)
#GTUCIS<- -cisS$t0.UCI - ((log(1-LoptUCIS/cisS$Linf.UCI))/cisS$K.UCI)
#GTLCIS<- -cisS$t0.LCI - ((log(1-LoptLCIS/cisS$Linf.LCI))/cisS$K.LCI)

##Pauly (1979) showed T0 is often unrealistic and biased so proposed a default computation for t0 based on Linf and K when estimating Generation time
##computes log(-t0) so take e to power to flip to true t0 estimate
PT0log<- -0.3922 - 0.2752*(log(cisS$Linf)) - 1.038*(log(cisS$K))
PT0<-exp(-PT0log)
PT0logUCI<- -0.3922 - 0.2752*(log(cisS$Linf.UCI)) - 1.038*(log(cisS$K.UCI))
PT0UCI<-exp(-PT0logUCI)
PT0logLCI<- -0.3922 - 0.2752*(log(cisS$Linf.LCI)) - 1.038*(log(cisS$K.LCI))
PT0LCI<-exp(-PT0logLCI)
###Now update Generation TIme calculation using Pauly (1979) T0, provides more realistic estimate and should be used
GTS<- PT0 - ((log(1-LoptS/cisS$Linf))/cisS$K)
GTUCIS<- PT0UCI - ((log(1-LoptUCIS/cisS$Linf.UCI))/cisS$K.UCI)
GTLCIS<- PT0LCI- ((log(1-LoptLCIS/cisS$Linf.LCI))/cisS$K.LCI)

GTSsum<-data.frame(years,LoptS,LoptLCIS,LoptUCIS,GTS,GTLCIS,GTUCIS)
excel(GTSsum)
meanGTS<-mean(GTS)
meanGTS
meanGTUCIS<-mean(GTUCIS)
meanGTUCIS
meanGTLCIS<-mean(GTLCIS)
meanGTLCIS
GTSsum[which.max(GTSsum$GTS),]
GTSsum[which.min(GTSsum$GTS),]
##Plot Lopt through time
GTSsum$years<-as.numeric(GTSsum$years)
ggplot(data=GTSsum,aes(x=years,y=LoptS, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=LoptLCIS,ymax=LoptUCIS),alpha=0.25) +
  scale_y_continuous(name= "Optimal Length (mm)", breaks = seq(140, 290, by = 10)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Spring LoptR.png", dpi=600, width=30, height= 30, units= "cm")



###Plot Generation Time through time
ggplot(data=GTSsum,aes(x=years,y=GTS, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=GTLCIS,ymax=GTUCIS),alpha=0.25) +
  scale_y_continuous(name= "Generation Time (years)", breaks = seq(4.2, 8.0, by = 0.4)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Spring Generation TimeR.png", dpi=600, width=30, height= 30, units= "cm")




#####Can fit the same VB models by fitting 1 model with year as a factor as suggested by Dan, this is the exact same as what was done above###
##This way I can extract nls results such as DF, RSE and p-vals
##Not necessary to use, only run if I want to report nls results##
#detS$year<-as.factor(detS$year)
#vbLKt <- length~Linf[year]*(1-exp(-K[year]*(age-t0[year])))
#grps<-unique(as.factor(detS$year))
#ngrps<-length(grps)
#sv0<-vbStarts(length~age, data=detS)
#svLKt<-Map(rep,sv0,c(34,34,34))
#fitLKt <- nls(vbLKt,data=detS,start=svLKt)
#summary(fitLKt)


##Fit a global VB with all data
vbLKtS <- length~Linf*(1-exp(-K*(age-t0)))
sv0S<-vbStarts(length~age, data=detS)
fitLKtS <- nls(vbLKtS,data=detS,start=sv0S)
summary(fitLKtS)
confint(fitLKtS)

################################################################################
####################################FALL########################################
################################################################################
###Subset for Fall and clean
detF<-subset(det, GSI_SPAWNING_SEASON == "A")
detF<-subset(detF, FROZEN_FISH_LENGTH <= 500)
detF<-subset(detF, FROZEN_FISH_LENGTH >= 2)
head(detF)
str(detF)
detF$length<-detF$FRESH_FISH_LENGTH
detF$age<-detF$FISH_AGE

##Fit VB
vb <- vbFuns()
predict2 <- function(x) predict(x,data.frame(age=ages))

##Get age range for each year
agesum<-group_by(detF, year) %>%
  summarize(minage=min(age), maxage=max(age))


#Specify year factor levels and length
#years<-unique(as.factor(detF$year))
years<-c("1988", "1989", "1990", "1991", "1992", "1993", "1994"
         , "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011"
         , "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021")
nyears<-length(years)


cfs <- cis <- preds1 <- preds2 <- NULL

for (i in 1:nyears) {
  ## Loop notification (for peace of mind)
  cat(years[i],"Loop\n")
  ## Isolate years's data
  tmp1 <- filter(detF,year==years[i])
  ## Fit von B to that year
  sv1 <- vbStarts(length~age,data=tmp1)    #Use 'detS' to have a common start value, and 'tmp1' for year specific start value for each VB model
  fit1 <- nls(length~vb(age,Linf,K,t0),data=tmp1,start=sv1)
  ## Extract and store parameter estimates and CIs
  cfs <- rbind(cfs,coef(fit1))
  boot1 <- Boot(fit1)
  tmp2 <-  confint(boot1)
  cis <- rbind(cis,c(tmp2["Linf",],tmp2["K",],tmp2["t0",]))
  ## Predict mean lengths-at-age with CIs
  ##   preds1 -> across all ages
  ##   preds2 -> across observed ages only
  ages <- seq(-1,17,0.2)
  boot2 <- Boot(fit1,f=predict2)
  tmp2 <- data.frame(year=years[i],age=ages,
                     predict(fit1,data.frame(age=ages)),
                     confint(boot2))
  preds1 <- rbind(preds1,tmp2)
  tmp2 <- filter(tmp2,age>=agesum$minage[i],age<=agesum$maxage[i])
  preds2 <- rbind(preds2,tmp2)
}

##Name rows and check VB parameters and predictions
rownames(cfs) <- rownames(cis) <- years
colnames(cis) <- paste(rep(c("Linf","K","t0"),each=2),
                       rep(c("LCI","UCI"),times=2),sep=".")
colnames(preds1) <- colnames(preds2) <- c("year","age","fit","LCI","UCI")
head(preds1) #preditions across all ages
head(preds2) #predictions across observed ages
head(cis) #VB parameters

###Merge year, VB paramters and VB parameter CIs into one dataframe
cisF<-data.frame(years,cfs, cis)
excel(cisF)
#Max and Min Linf over time series
cisF[which.max(cisF$Linf),]
cisF[which.min(cisF$Linf),]

##QUickly plot Linf over time series
plot(cisF$years, cisF$Linf)
lines(cisF$years, cisF$Linf)

#plot all VB over timeseries
ggplot() +
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI),alpha=0.25) +
  geom_point(data=detF,aes(x=age,y=length),
             size=2,alpha=0.3) +
  geom_line(data=preds1,aes(x=age,y=fit),
            size=1,linetype="dashed") +
  geom_line(data=preds2,aes(x=age,y=fit),
            size=1) +
  scale_y_continuous(name="Total Length (mm)",limits=c(0,450)) +
  scale_x_continuous(name="Age (years)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20),
        axis.text = element_text(size=14, colour= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.position="none") +
  facet_wrap(vars(year))
ggsave("Fall VB CurvesR.png", dpi=600, width=30, height= 30, units= "cm")

###Plot all VB on same graph to show changes over years
ggplot() +
  geom_ribbon(data=preds2,aes(x=age,ymin=LCI,ymax=UCI,fill=year),alpha=0.25) +
  geom_point(data=detF,aes(x=age,y=length,color=year),
             size=2,alpha=0.3) +
  geom_line(data=preds1,aes(x=age,y=fit,color=year),
            size=1,linetype="dashed") +
  geom_line(data=preds2,aes(x=age,y=fit,color=year),
            size=1) +
  scale_y_continuous(name="Total Length (mm)",limits=c(0,450)) +
  scale_x_continuous(name="Age (years)",breaks=0:17) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Fall VB Curves 1PlotR.png", dpi=600, width=30, height= 30, units= "cm")

###Plot Linf over time
cisF$years<-as.numeric(cisF$years)
ggplot(data=cisF,aes(x=years,y=Linf, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=Linf.LCI,ymax=Linf.UCI),alpha=0.25) +
  scale_y_continuous(name= "Asymptotic Length (mm)", breaks = seq(290, 460, by = 10)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

ggsave("Fall LinfR.png", dpi=600, width=30, height= 30, units= "cm")


##Load in Mortality estimates, which are seperated y N, M and S for Fall Herring
#Must load them then take average of all for use in Generation time estimates
f=lisread("scaFqm.rep")

#North
baa = f$BtaN

mctmp <- read.table("mcoutMn1.dat", sep = "")
m1 <- apply( mctmp, 2, median)
m1lci <- apply( mctmp, 2, quantile, 0.025)
m1uci <- apply( mctmp, 2, quantile, 0.975)

mctmp <- read.table("mcoutMn2.dat", sep = "")
m2 <- apply( mctmp, 2, median)
m2lci <- apply( mctmp, 2, quantile, 0.025)
m2uci <- apply( mctmp, 2, quantile, 0.975)

#build df
m = baa #just to get proper sized object
m[,1:5] = m1
m[,6:10] = m2
#m

muci = baa #just to get proper sized object
muci[,1:5] = m1uci
muci[,6:10] = m2uci

mlci = baa #just to get proper sized object
mlci[,1:5] = m1lci
mlci[,6:10] = m2lci

#do abundance weighted average
mprop = m*(baa/rowSums(baa))
m.aw = rowSums(mprop)

mprop = muci*(baa/rowSums(baa))
muci.aw = rowSums(mprop)

mprop = mlci*(baa/rowSums(baa))
mlci.aw = rowSums(mprop)

yr = 1978:2021
plot(yr, m.aw, ylim = c(0,0.6), type = "l")
lines(yr, muci.aw, col = "red")
lines(yr, mlci.aw, col = "red")

north = data.frame(median = m.aw, ucl = muci.aw, lcl = mlci.aw)
#



#middle
baa = f$BtaM

mctmp <- read.table("mcoutMm1.dat", sep = "")
m1 <- apply( mctmp, 2, median)
m1lci <- apply( mctmp, 2, quantile, 0.025)
m1uci <- apply( mctmp, 2, quantile, 0.975)

mctmp <- read.table("mcoutMm2.dat", sep = "")
m2 <- apply( mctmp, 2, median)
m2lci <- apply( mctmp, 2, quantile, 0.025)
m2uci <- apply( mctmp, 2, quantile, 0.975)

#build df
m = baa #just to get proper sized object
m[,1:5] = m1
m[,6:10] = m2
#m

muci = baa #just to get proper sized object
muci[,1:5] = m1uci
muci[,6:10] = m2uci

mlci = baa #just to get proper sized object
mlci[,1:5] = m1lci
mlci[,6:10] = m2lci

#do abundance weighted average
mprop = m*(baa/rowSums(baa))
m.aw = rowSums(mprop)

mprop = muci*(baa/rowSums(baa))
muci.aw = rowSums(mprop)

mprop = mlci*(baa/rowSums(baa))
mlci.aw = rowSums(mprop)

yr = 1978:2021
plot(yr, m.aw, ylim = c(0,0.6), type = "l")
lines(yr, muci.aw, col = "red")
lines(yr, mlci.aw, col = "red")

middle = data.frame(median = m.aw, ucl = muci.aw, lcl = mlci.aw)


#####south
baa = f$BtaS

mctmp <- read.table("mcoutMs1.dat", sep = "")
m1 <- apply( mctmp, 2, median)
m1lci <- apply( mctmp, 2, quantile, 0.025)
m1uci <- apply( mctmp, 2, quantile, 0.975)

mctmp <- read.table("mcoutMs2.dat", sep = "")
m2 <- apply( mctmp, 2, median)
m2lci <- apply( mctmp, 2, quantile, 0.025)
m2uci <- apply( mctmp, 2, quantile, 0.975)

#build df
m = baa #just to get proper sized object
m[,1:5] = m1
m[,6:10] = m2
#m

muci = baa #just to get proper sized object
muci[,1:5] = m1uci
muci[,6:10] = m2uci

mlci = baa #just to get proper sized object
mlci[,1:5] = m1lci
mlci[,6:10] = m2lci

#do abundance weighted average
mprop = m*(baa/rowSums(baa))
m.aw = rowSums(mprop)

mprop = muci*(baa/rowSums(baa))
muci.aw = rowSums(mprop)

mprop = mlci*(baa/rowSums(baa))
mlci.aw = rowSums(mprop)

yr = 1978:2021
plot(yr, m.aw, ylim = c(0,0.6), type = "l")
lines(yr, muci.aw, col = "red")
lines(yr, mlci.aw, col = "red")

south = data.frame(median = m.aw, ucl = muci.aw, lcl = mlci.aw)


#do biomass weighted avergae M over the 3 regions
Mtot = (north$median*f$SSBt[1,]/colSums(f$SSBt))+(middle$median*f$SSBt[2,]/colSums(f$SSBt))+(south$median*f$SSBt[3,]/colSums(f$SSBt))   

MtotUCI = (north$ucl*f$SSBt[1,]/colSums(f$SSBt))+(middle$ucl*f$SSBt[2,]/colSums(f$SSBt))+(south$ucl*f$SSBt[3,]/colSums(f$SSBt))   
MtotLCI = (north$lcl*f$SSBt[1,]/colSums(f$SSBt))+(middle$lcl*f$SSBt[2,]/colSums(f$SSBt))+(south$lcl*f$SSBt[3,]/colSums(f$SSBt))   

plot(yr, fall$median, ylim = c(0, 0.5), type = "l")
lines(yr, fall$ucl, col = "red")
lines(yr, fall$lcl, col = "red")

fall
Myears<-c("1978", "1979", "1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989", "1990", "1991", "1992", "1993", "1994"
          , "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011"
          , "2012", "2013", "2014", "2015", "2016", "2017", "2018","2019", "2020", "2021")

FallM<-data.frame(Myears,Mtot, MtotUCI, MtotLCI)
FallM<-subset(FallM, Myears >= "1988")
###Compute Lopt
LoptF<-cisF$Linf* (3/(3+(FallM$Mtot/cisF$K)))
LoptUCIF<-cisF$Linf.UCI* (3/(3+(FallM$MtotUCI/cisF$K.UCI)))
LoptLCIF<-cisF$Linf.LCI* (3/(3+(FallM$MtotLCI/cisF$K.LCI)))
##Compute Generation Time with estimated T0 from VB models
#GTF<- -cisF$t0 - ((log(1-LoptF/cisF$Linf))/cisF$K)
#GTUCIF<- -cisF$t0.UCI - ((log(1-LoptUCIF/cisF$Linf.UCI))/cisF$K.UCI)
#GTLCIF<- -cisF$t0.LCI - ((log(1-LoptLCIF/cisF$Linf.LCI))/cisF$K.LCI)

##Pauly (1979) showed T0 is often unrealistic and biased so proposed a default computation for t0 based on Linf and K when estimating Generation time
##computes log(-t0) so take e to power to flip to true t0 estimate
PT0log<- -0.3922 - 0.2752*(log(cisF$Linf)) - 1.038*(log(cisF$K))
PT0<-exp(-PT0log)
PT0logUCI<- -0.3922 - 0.2752*(log(cisF$Linf.UCI)) - 1.038*(log(cisF$K.UCI))
PT0UCI<-exp(-PT0logUCI)
PT0logLCI<- -0.3922 - 0.2752*(log(cisF$Linf.LCI)) - 1.038*(log(cisF$K.LCI))
PT0LCI<-exp(-PT0logLCI)
###Now update Generation TIme calculation using Pauly (1979) T0, provides more realistic estimate and should be used
GTF<- PT0 - ((log(1-LoptF/cisF$Linf))/cisF$K)
GTUCIF<- PT0UCI - ((log(1-LoptUCIF/cisF$Linf.UCI))/cisF$K.UCI)
GTLCIF<- PT0LCI - ((log(1-LoptLCIF/cisF$Linf.LCI))/cisF$K.LCI)
#Compute mean across time series
meanGTF<-mean(GTF)
meanGTF
meanGTUCIF<-mean(GTUCIF)
meanGTUCIF
meanGTLCIF<-mean(GTLCIF)
meanGTLCIF
#Collate into data frame
GTFsum<-data.frame(years,LoptF,LoptLCIF,LoptUCIF,GTF,GTLCIF,GTUCIF)
excel(GTFsum)
#Max and Min across time series
GTFsum[which.max(GTFsum$GTF),]
GTFsum[which.min(GTFsum$GTF),]

##Plot Lopt through time
GTFsum$years<-as.numeric(GTFsum$years)
ggplot(data=GTFsum,aes(x=years,y=LoptF, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=LoptLCIF,ymax=LoptUCIF),alpha=0.25) +
  scale_y_continuous(name= "Optimal Length (mm)", breaks = seq(140, 340, by = 10)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Fall LoptR.png", dpi=600, width=30, height= 30, units= "cm")



###Plot Generation Time through time
ggplot(data=GTFsum,aes(x=years,y=GTF, group=1)) +
  geom_point(size=2) +
  geom_line(size=1.2) +
  geom_ribbon(aes(x=years,ymin=GTLCIF,ymax=GTUCIF),alpha=0.25) +
  scale_y_continuous(name= "Generation Time (years)", breaks = seq(5.0, 11.0, by = 0.4)) +
  scale_x_continuous(name = "Year", breaks = seq(1988, 2021, by = 3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        axis.text=element_text(size=14, colour="black"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave("Fall Generation TimeR.png", dpi=600, width=30, height= 30, units= "cm")



#####Can fit the same VB models by fitting 1 model with year as a factor as suggested by Dan, this is the exact same as what was done above###
##This way I can extract nls results such as DF, RSE and p-vals
##Not necessary to use, only run if I want to report nls results##
#detF$year<-as.factor(detF$year)
#vbLKt <- length~Linf[year]*(1-exp(-K[year]*(age-t0[year])))
#grps<-unique(as.factor(detF$year))
#ngrps<-length(grps)
#sv0<-vbStarts(length~age, data=detF)
#svLKt<-Map(rep,sv0,c(34,34,34))
#fitLKt <- nls(vbLKt,data=detF,start=svLKt)
#summary(fitLKt)

##Fit a global VB with all data
vbLKtF <- length~Linf*(1-exp(-K*(age-t0)))
sv0F<-vbStarts(length~age, data=detF)
fitLKtF <- nls(vbLKtF,data=detF,start=sv0F)
summary(fitLKtF)
confint(fitLKtF)


