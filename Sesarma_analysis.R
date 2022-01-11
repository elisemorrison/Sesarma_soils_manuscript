# Script for analyzing elemental and stable isotope data for Sesarma manuscript

library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggsci)
library(dunn.test)
library(lmerTest)
library(plotrix)

## Import data
raw_EA_isotope_data<-read.csv("./Data/ElementalAnalysis_Isotopes.csv")

#create a new object to work with, so you don't manipulate the raw data
bulk<-raw_EA_isotope_data

#add creek ID for manuscript, based on creek name
bulk<-mutate(bulk, CreekID=case_when(Site == "North" ~ "A", 
                                                   Site == "East" ~ "B", 
                                                   Site == "Airport" ~ "C", 
                                                   Site == "BR" ~ "D")) %>% 
  group_by(SampleID) %>% 
  mutate(CN=case_when(Location=='Benthic diatoms'~wt..TC/wt..TN, 
                      Location!='Benthic diatoms'~acidified.wt..C/wt..TN), 
         d13C=case_when(Location=='Benthic diatoms'~d13C_not_acidified, 
                        Location!='Benthic diatoms'~acidified.d13C..permil..vs.VPDB.))



################# Analyze USS

USS<-filter(bulk, !is.na(USS.depth.cm)) %>% 
  select(SampleID, Position, Location, CreekID, Status, Replicate, SoilType, USS.depth.cm) %>% 
  unique()

#check variance
var(USS$USS.depth.cm)
summary(USS$USS.depth.cm)
hist(USS$USS.depth.cm)

#given the structure of the dataset, best to use Kruskal Wallis test
USS_KW<-kruskal.test(USS.depth.cm~Status, data=USS)

dunn.test(USS$USS.depth.cm, USS$Status)

summary_USS<-USS %>% 
  group_by(Position, Status) %>%
  summarize(depthmean=mean(USS.depth.cm, na.rm=TRUE), 
            depthsd=sd(USS.depth.cm, na.rm=TRUE), 
            depthmin=min(USS.depth.cm, na.rm=TRUE), 
            depthmax=max(USS.depth.cm, na.rm=TRUE))

plotA<-ggplot(USS, aes(Position, USS.depth.cm, color=Status))+
  theme_minimal()+
  ylab("USS depth (cm)")+
  xlab("Position along Creek")+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7"))+
  theme(strip.text.x = element_text(size = 10))+
  panel_border(color="black")+
  geom_point(aes(Position, USS.depth.cm), position=position_dodge(width=0.8), alpha=0.2)+
  stat_summary(aes(Position, USS.depth.cm), position=position_dodge(width=0.8), size=1)+
  scale_color_npg()


plotB<-ggplot(USS, aes(Position, USS.depth.cm, color=Status))+
  theme_minimal()+
  ylab("USS depth (cm)")+
  xlab("Position along Creek")+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7"))+
  theme(strip.text.x = element_text(size = 10))+
  panel_border(color="black")+
  geom_point(aes(Position, USS.depth.cm), position=position_dodge(width=0.8), alpha=0.2)+
  stat_summary(aes(Position, USS.depth.cm), position=position_dodge(width=0.8), size=0.5)+
  facet_grid(.~CreekID)+
  scale_color_npg()+
  theme(legend.position="none")

plot_grid(plotA, plotB, labels=c("A", "B"), ncol = 1)

############### Create summary table for bulk measures

#rather than grouping by soil type, get mean values from all soil types
table1_allcreeks_aggregated<-select(bulk, SampleID, SoilType, Status, Position, Location, 
                                    Site, Status, Replicate, USS.depth.cm, d15N, 
                                    acidified.d13C..permil..vs.VPDB., d13C,
                                    d13C_not_acidified,
                                    wt..TN, wt..TC, acidified.wt..C, CN) %>% 
  group_by(Status, SoilType, Position) %>% 
  summarize(meanCN=round(mean(CN, na.rm=TRUE), 2), 
            sdCN=round(sd(CN, na.rm=TRUE), 2),
            mean15N=round(mean(d15N, na.rm=TRUE), 2), 
            sd15N=round(sd(d15N, na.rm=TRUE), 2), 
            meanTC=round(mean(wt..TC, na.rm=TRUE), 2), 
            sdTC=round(sd(wt..TC, na.rm=TRUE), 2), 
            meanTN=round(mean(wt..TN, na.rm=TRUE), 2), 
            sdTN=round(sd(wt..TN, na.rm=TRUE), 2), 
            mean13C=round(mean(d13C, na.rm=TRUE), 2), 
            sd13C=round(sd(d13C, na.rm=TRUE), 2), 
            meanTOC=round(mean(acidified.wt..C, na.rm=TRUE), 2), 
            sdTOC=round(sd(acidified.wt..C, na.rm=TRUE), 2), 
            n=n()) 

write.csv(table1_allcreeks_aggregated, "./Tables/TableS1.csv")

############### Test for significant differences between fore- and hind-guts

guts<-filter(bulk, Position=="Fore-gut"|Position=="Hind-gut") %>% 
  group_by(Position) %>% 
  mutate(CN=acidified.wt..C/wt..TC)
#First run with TC
#F test to compare variances
var.test(wt..TC~Position, data=guts) #ratio of variances=0.4259082 p=0.425
t.test(wt..TC~Position, data=guts) #p<2.2 e-16
wilcox.test(wt..TC~Position, data=guts) #W = 208, p-value = 2.947e-08, Wilcoxon exact

#Now TOC
var.test(acidified.wt..C~Position, data=guts) #ratio of variances=0.266312 p=0.0299
wilcox.test(acidified.wt..C~Position, data=guts) #W = 169, p-value = 1.644e-05, continuity correction

#Now TN
var.test(wt..TN~Position, data=guts) #ratio of variances=28.60338 p=8.838e-07
wilcox.test(wt..TN~Position, data=guts) #W = 208, p-value = 5.613e-06, continuity correction

#Now d13C 
var.test(acidified.d13C..permil..vs.VPDB.~Position, data=guts) #ratio of variances=1.28663 p-value = 0.6694
wilcox.test(acidified.d13C..permil..vs.VPDB.~Position, data=guts) #W = 26, p-value = 0.002907

#Now d15N
var.test(d15N~Position, data=guts) #ratio of variances=1.581623 p-value = 0.429
wilcox.test(d15N~Position, data=guts) 

#Now C:N
var.test(CN~Position, data=guts) #ratio of variances=0.05375183 p-value = 1.503e-06
wilcox.test(CN~Position, data=guts) #W = 0, p-value = 2.947e-08

############### Plot Figure 3


############### Evaluate differences between grazed and ungrazed sediments

############### Use a simple two-endmember mixing model (Middelburg, 2014)
############### to evaluate the proportion of C from spartina (Resource A)

mean.d13C<-bulk %>% 
  filter(Site=="Soft tissue"|Site=="Benthic diatoms"|Site=="Spartina") %>% 
  select(Site, d13C_not_acidified, acidified.d13C..permil..vs.VPDB.) %>% 
  group_by(Site) %>% 
  mutate(meand13C_na=mean(d13C_not_acidified, na.rm=TRUE),
         meand13C_a=mean(acidified.d13C..permil..vs.VPDB., na.rm=TRUE)) %>% 
  select(-d13C_not_acidified, -'acidified.d13C..permil..vs.VPDB.') %>% 
  unique()


Ciso<-mean.d13C %>% 
  mutate(numerator= -14.31578--19.55496, 
         denom=-13.78658--19.55496, 
         pc=numerator/denom)

## % Carbon contributed by Spartina 90.8%
