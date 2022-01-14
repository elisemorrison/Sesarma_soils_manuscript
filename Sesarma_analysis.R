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
                        Location!='Benthic diatoms'~acidified.d13C..permil..vs.VPDB.)) %>% 
  mutate(Category=case_when(SoilType=="Upper Horizon" ~ "USS",
                            SoilType=="Lower Horizon" ~ "Sediment", 
                            SoilType=="No Horizon" ~ "Sediment",
                            SoilType=="Benthic diatoms"~"Benthic diatoms", 
                            SoilType=="Whole crab"~"Sesarma", 
                            SoilType=="Soft tissue"~"Sesarma",
                            SoilType=="Gut" ~ "Sesarma", 
                            SoilType=="Spartina"~"Spartina",
                            SoilType=="Visceral mass"~"Mussel", 
                            SoilType=="Gills"~"Mussel"))



################# Analyze USS

USS<-filter(bulk, !is.na(USS.depth.cm)) %>% 
  select(SampleID, Position, Location, CreekID, Status, Replicate, Category, USS.depth.cm, 
         wt..TN, wt..TC, d15N, d13C, acidified.wt..C, CN) %>% 
  unique() %>% 
  rename(Sample=SampleID)

#check variance
var(USS$USS.depth.cm)
summary(USS$USS.depth.cm)
hist(USS$USS.depth.cm)

#given the structure of the dataset, best to use Kruskal Wallis test
USS_KW<-kruskal.test(USS.depth.cm~Status, data=USS)

dunn.test(USS$USS.depth.cm, USS$Status)

#summary table (Table S3) for USS
summary_USS<-USS %>% 
  group_by(Position, Status) %>%
  summarize(depthmean=mean(USS.depth.cm, na.rm=TRUE), 
            depthsd=sd(USS.depth.cm, na.rm=TRUE), 
            depthmin=min(USS.depth.cm, na.rm=TRUE), 
            depthmax=max(USS.depth.cm, na.rm=TRUE))

write.csv(summary_USS, "./Tables/TableS3.csv")

#supplementary plot (S1) for USS
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
table_allcreeks_aggregated<-select(bulk, SampleID, Category, Status, Position, Location, 
                                    Site, Replicate, USS.depth.cm, d15N, 
                                    d13C,wt..TN, wt..TC, acidified.wt..C, CN) %>% 
  group_by(Status, Category, Position) %>% 
  summarize(meanUSSDepth=round(mean(USS.depth.cm, na.rm=TRUE), 1), 
            sdUSSDepth=round(sd(USS.depth.cm, na.rm=TRUE), 1),
            meanCN=round(mean(CN, na.rm=TRUE), 2), 
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

write.csv(table_allcreeks_aggregated, "./Tables/TableS_bulk.csv")

############### Random effects model for sediments

#subset for sediments first
sediment<-filter(bulk, Category=="Sediment")

TOC_sediment_lmer<-lmerTest::lmer(acidified.wt..C~Status+(1|Position), data=sediment)
anova(TOC_sediment_lmer)

#test for TN
TN_sediment_lmer<-lmerTest::lmer(wt..TN~Status+(1|Position), data=sediment)
anova(TN_sediment_lmer)

#test for TC
TC_sediment_lmer<-lmerTest::lmer(wt..TC~Status+(1|Position), data=sediment)
anova(TC_sediment_lmer)

#test for C:N
CN_sediment_lmer<-lmerTest::lmer(CN~Status+(1|Position), data=sediment)
anova(CN_sediment_lmer)

#test for d13C
d13C_sediment_lmer<-lmerTest::lmer(acidified.d13C..permil..vs.VPDB.~Status+(1|Position), data=sediment)
anova(d13C_sediment_lmer)

#test for d15N
d15N_sediment_lmer<-lmerTest::lmer(d15N~Status+(1|Position), data=sediment)
anova(d15N_sediment_lmer)

#sediment aggregated across zones
table_nozones<-select(sediment, SampleID, Category, Status, Position, Location, 
                                   Site, Replicate, USS.depth.cm, d15N, 
                                   d13C,wt..TN, wt..TC, acidified.wt..C, CN) %>% 
  group_by(Status) %>% 
  summarize(meanUSSDepth=round(mean(USS.depth.cm, na.rm=TRUE), 1), 
            sdUSSDepth=round(sd(USS.depth.cm, na.rm=TRUE), 1),
            meanCN=round(mean(CN, na.rm=TRUE), 2), 
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

############### Plot figures for sediments and endmembers

## take out mussel samples
renamed_for_plots<-filter(bulk, Category!="Mussel")

#update order that Biodeposit and Type will plot
renamed_for_plots$Position<-gsub("Benthic diatoms", "", renamed_for_plots$Position) #workaround to not have duplicate labels
renamed_for_plots$Position<-gsub("Spartina", "", renamed_for_plots$Position) 
renamed_for_plots$Status<-ordered(renamed_for_plots$Status, c("Grazed", "Un-grazed", "Benthic diatoms", "Whole crab", "Soft tissue", "Fore-gut", 
                                                          "Hind-gut", "Spartina"))
renamed_for_plots$Category<-ordered(renamed_for_plots$Category, c("USS", "Sediment", "Sesarma", "Spartina", "Benthic diatoms"))

# ## first want to create legend for shared "bulk data" plot...just need a temp plot to extract legend from
temp_plot<-ggplot(renamed_for_plots, aes(Position, d15N, fill=Category))+
  geom_boxplot(position="dodge", color="black")+ 
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  guides(fill=guide_legend(title="Sample Type"))+
  ggtitle(expression(paste(delta,""^15,"N ")))+
  theme(strip.text.x = element_text())+
  panel_border(color="black")+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  theme(text = element_text(size = 20))

bulk_legend<-get_legend(temp_plot)


CN_all<-ggplot(renamed_for_plots, aes(Position, CN, fill=Category))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  ylab("C:N")+
  xlab("")+ #only have x label for bottom plot; d13C will be bottom plot
  theme(text = element_text(size = 20))+
  theme(axis.text.x = element_blank())+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free_x", nrow=1)+ 
  theme(legend.position = "none") #have no legend so bulk legend can be added to full bulk plot
 

CN_all


d15N_all<-ggplot(renamed_for_plots, aes(Position, d15N, fill=Category))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  ylab(expression(paste(delta,""^15,"N (‰)")))+
  xlab("")+
  theme(text = element_text(size = 20))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free_x", nrow=1)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())

d15N_all

TC_all<-ggplot(renamed_for_plots, aes(Position, wt..TC, fill=Category))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  ylab("Total Carbon (%)")+
  xlab("")+
  theme(text = element_text(size = 20))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free", nrow=1)+ 
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))

TC_all

TN_all<-ggplot(renamed_for_plots, aes(Position, wt..TN, fill=Category))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  ylab("Total Nitrogen (%)")+
  xlab("")+
  theme(text = element_text(size = 20))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free", nrow=1)+ 
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))

TN_all

d13C_all<-ggplot(renamed_for_plots, aes(Position, d13C, fill=Category))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  ylab(expression(paste(delta,""^13,"C (‰)")))+
  xlab("")+
  theme(text = element_text(size = 20))+
  panel_border(color="black")+
  facet_grid(.~Category, scales="free_x")+ 
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

d13C_all

### plot all bulk measures together
### plot one with stable isotopes, another with C, N, C:N

plot_grid(TC_all, TN_all, CN_all, bulk_legend, ncol=1, 
          rel_heights = c(1,1,1.1,0.1))


# plot one with stable isotopes, another with C, N, C:N
plot_grid(CN_all, d15N_all, d13C_all, bulk_legend,  ncol=1, 
          rel_heights = c(1,1,1.3,0.1), labels=c("A", "B", "C", ""))

############### Plot C and N isotopes

isoplot_data<-renamed_for_plots %>% 
  mutate(Status=case_when(Position=="Fore-gut" ~ "Fore-gut", 
                   Position=="Hind-gut" ~ "Hind-gut",
                   Position=="Whole crab" ~ "Whole crab",
                   Position=="Soft tissue" ~ "Soft tissue",
                   Status=="Grazed" ~ "Grazed", 
                   Status=="Un-grazed" ~ "Un-grazed",
                   Status=="Spartina" ~ "Spartina",
                   Status=="Benthic diatoms" ~ "Benthic diatoms"
                   )) %>% 
  group_by(Category, Status) %>% 
  mutate(meand15N=mean(d15N, na.rm=TRUE), 
         meand13C=mean(d13C, na.rm=TRUE), 
         sd15N=sd(d15N, na.rm=TRUE), 
         sd13C=sd(d13C, na.rm=TRUE)) %>% 
  select(Status, Category, d13C, d15N, sd15N, sd13C, meand13C, meand15N) %>% 
  unique() %>% 
  rename(`Sample Type`=Category, Description=Status)

#reorder so it plots more intuitively
isoplot_data$Description <-ordered(isoplot_data$Description, c("Grazed", "Un-grazed", "Benthic diatoms", "Whole crab", "Soft tissue", "Fore-gut", 
                                                              "Hind-gut", "Spartina"))

## get the legend
isoplot_legend<-get_legend(ggplot(isoplot_data)+
                             geom_point(aes(d13C, d15N, fill=Description, color=Description), alpha=0.5, size=1)+
                             geom_errorbar(aes(x=meand13C, ymin=meand15N-sd15N, ymax=meand15N+sd15N, color=Description), size=0.5)+
                             geom_errorbarh(aes(y=meand15N, xmin=meand13C-sd13C, xmax=meand13C+sd13C, color=Description), size=0.5)+
                             geom_point(aes(meand13C, meand15N, fill=Description, color=Description, shape=`Sample Type`), size=3)+ #keep this line for now, for plotting the legend
                             theme_minimal()+
                             scale_fill_npg()+
                             scale_color_npg()+
                             scale_shape_manual(values=c(21, 22, 23, 24, 25, 12))+
                             theme(text = element_text(size = 15)))

## now save actual plot, with black outlines around mean values (if we kept the legend, it doesn't plot the colors clearly)
isoplot<-ggplot(isoplot_data)+
  geom_point(aes(d13C, d15N, fill=Description, color=Description, shape=`Sample Type`), alpha=0.5, size=3)+
  geom_errorbar(aes(x=meand13C, ymin=meand15N-sd15N, ymax=meand15N+sd15N, color=Description), size=0.5)+
  geom_errorbarh(aes(y=meand15N, xmin=meand13C-sd13C, xmax=meand13C+sd13C, color=Description), size=0.5)+
  geom_point(aes(meand13C, meand15N, fill=Description, shape=`Sample Type`), color="black",size=5)+ #keep this line for now, for plotting the legend
  theme_minimal()+
  scale_fill_npg()+
  scale_color_npg()+
  scale_shape_manual(values=c(21, 22, 23, 24, 25, 12))+
  theme(legend.position="none")+
  ylab(expression(paste(delta,""^15,"N (‰)")))+
  xlab(expression(paste(delta,""^13,"C (‰)")))+
  panel_border(color = "grey85", size = 1, linetype = 1)+
  theme(text = element_text(size = 20))

### final isotope plot
plot_grid(isoplot, isoplot_legend, rel_widths = c(1, 0.25))


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


############### Amino acid analyses
###############



