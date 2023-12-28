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
                            SoilType=="Spartina"~"Spartina"))



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

plot_grid(plotA, plotB, labels=c("(a)", "(b)"), ncol = 1)

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

write.csv(table_allcreeks_aggregated, "./Tables/TableS_bulk_creeksAggregated.csv")

#get endmember supplementary table, need to make sex/sample codes more reader-friendly
endmembers<-filter(bulk, Status!="Grazed"&Status!="Un-grazed") %>% 
  select(SoilType, Location, Position, Status, CN, wt..TC, acidified.wt..C, wt..TN, d13C, d15N) %>% 
  mutate(Other=case_when(Location=="LM"~"Large male", 
                         Location=="LG"~"Large gravid", 
                         Location=="MM"~"Medium male", 
                         Location=="MG"~"Medium gravid", 
                         Location=="MF"~"Medium female",
                         Location=="SM"~"Small male", 
                         Location=="SG"~"Small gravid", 
                         Location=="SF"~"Small female", 
                         Status=="acid rinsed" ~"acid rinsed", 
                         Status=="non-acid rinsed" ~"non-acid rinsed", 
  )) %>% 
  select(SoilType, Position, Other, CN, wt..TC, acidified.wt..C, wt..TN, d13C, d15N)

table_endmember<-select(endmembers, SampleID, SoilType, Position, Other, d15N, 
                                   d13C,wt..TN, wt..TC, acidified.wt..C, CN) %>% 
  group_by(SoilType, Position, Other) %>% 
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

endmember_formatted<-mutate(table_endmember,
                             CN=paste0(meanCN, " (", sdCN, ")"),
                             TC=paste0(meanTC, " (", sdTC, ")"), 
                             TOC=paste0(meanTOC, " (", sdTOC, ")"), 
                             TN=paste0(meanTN, " (", sdTN, ")"), 
                             d13C=paste0(mean13C, " (", sd13C, ")"), 
                             d15N=paste0(mean15N, " (", sd15N, ")")
) %>% 
  select(SoilType, Position, Other, CN, TC, TOC, TN, d13C, d15N, n)

write.csv(endmember_formatted, "./Tables/TableS_bulk_endmember.csv")

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

## Give a color code to sample types
renamed_for_plots<-bulk %>% 
  mutate(ColorCode=case_when(Position=="Fore-gut" ~ "Fore-gut", 
                                             Position=="Hind-gut" ~ "Hind-gut",
                                             Position=="Whole crab" ~ "Whole crab",
                                             Position=="Soft tissue" ~ "Soft tissue",
                                             Category=="USS" ~ "USS", 
                                             Category=="Sediment" ~ "Sediment",
                                             Status=="Spartina" ~ "Spartina",
                                             Status=="Benthic diatoms" ~ "Benthic diatoms"))

#update order that Biodeposit and Type will plot
renamed_for_plots$Position<-gsub("Benthic diatoms", "", renamed_for_plots$Position) #workaround to not have duplicate labels
renamed_for_plots$Position<-gsub("Spartina", "", renamed_for_plots$Position) 
renamed_for_plots$Status<-ordered(renamed_for_plots$Status, c("Grazed", "Un-grazed", "Benthic diatoms", "Whole crab", "Soft tissue", "Fore-gut", 
                                                          "Hind-gut", "Spartina"))
renamed_for_plots$ColorCode<-ordered(renamed_for_plots$ColorCode, c("USS", "Sediment", "Benthic diatoms", "Whole crab", "Soft tissue", "Fore-gut", 
                                                              "Hind-gut", "Spartina"))
renamed_for_plots$Category<-ordered(renamed_for_plots$Category, c("USS", "Sediment", "Benthic diatoms", "Sesarma", "Spartina"))


#CN has one fewer category, causing colors to not align
#hardcode colors for this plot
Palette <- c("#3C5488FF", "#7E6148FF", "#42B540FF", "#B09C85FF",  
             "#F39B7FFF", "#8491B4FF", "#00A087FF", "#91D1C2FF")

CNPalette<- c("#3C5488FF", "#7E6148FF", "#42B540FF", "#B09C85FF",  
              "#8491B4FF", "#00A087FF", "#91D1C2FF")

# ## first want to create legend for shared "bulk data" plot...just need a temp plot to extract legend from
temp_plot<-ggplot(renamed_for_plots, aes(Position, d15N, fill=ColorCode))+
  geom_boxplot(position="dodge", color="black")+ 
  theme_minimal()+
  scale_fill_manual(values=Palette)+
  ggtitle(expression(paste(delta,""^15,"N ")))+
  theme(strip.text.x = element_text())+
  panel_border(color="black")+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  theme(text = element_text(size = 20))+
  guides(fill=guide_legend(title="Sample Type"))

bulk_legend<-get_legend(temp_plot)

CN_all<-ggplot(renamed_for_plots, aes(Position, CN, fill=ColorCode))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_manual(values=CNPalette)+
  ylab("C:N")+
  xlab("")+ #only have x label for bottom plot; d13C will be bottom plot
  theme(text = element_text(size = 20))+
  theme(axis.text.x = element_blank(), 
        axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free_x", nrow=1)+ 
  theme(legend.position = "none") #have no legend so bulk legend can be added to full bulk plot

CN_all


d15N_all<-ggplot(renamed_for_plots, aes(Position, d15N, fill=ColorCode))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_manual(values=Palette)+
  ylab(expression(paste(delta,""^15,"N (‰)")))+
  xlab("")+
  theme(text = element_text(size = 20))+
  theme(axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free_x", nrow=1)+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())

d15N_all

TC_all<-ggplot(renamed_for_plots, aes(Position, wt..TC, fill=ColorCode))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_manual(values=Palette)+
  ylab("Total Carbon (%)")+
  xlab("")+
  theme(text = element_text(size = 20))+
  theme(axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free", nrow=1)+ 
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))

TC_all

TN_all<-ggplot(renamed_for_plots, aes(Position, wt..TN, fill=ColorCode))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+
  theme_minimal()+
  scale_fill_manual(values=Palette)+
  ylab("Total Nitrogen (%)")+
  xlab("")+
  theme(text = element_text(size = 20))+
  panel_border(color="black")+
  facet_wrap(.~Category, scales="free", nrow=1)+ 
  theme(legend.position = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))

TN_all

d13C_all<-ggplot(renamed_for_plots, aes(Position, d13C, fill=ColorCode))+
  geom_boxplot(position=position_dodge(preserve="single"), color="black")+ 
  theme_minimal()+
  scale_fill_manual(values=Palette)+
  ylab(expression(paste(delta,""^13,"C (‰)")))+
  xlab("")+
  theme(text = element_text(size = 20))+
  theme(axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))+
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
plot_grid(CN_all, d15N_all, d13C_all, ncol=1, 
          rel_heights = c(1,1,1.3,0.2), labels=c("(a)", "(b)", "(c)", ""))



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

IsoPalette <- c("#E64B35FF", "#4DBBD5FF", "#42B540FF", "#B09C85FF", 
                "#F39B7FFF", "#8491B4FF", "#00A087FF", "#91D1C2FF")
                

## get the legend
isoplot_legend<-get_legend(ggplot(isoplot_data)+
                             geom_point(aes(d13C, d15N, fill=Description, color=Description), alpha=0.5, size=1)+
                             geom_errorbar(aes(x=meand13C, ymin=meand15N-sd15N, ymax=meand15N+sd15N, color=Description), size=0.5)+
                             geom_errorbarh(aes(y=meand15N, xmin=meand13C-sd13C, xmax=meand13C+sd13C, color=Description), size=0.5)+
                             geom_point(aes(meand13C, meand15N, fill=Description, color=Description, shape=`Sample Type`), size=3)+ #keep this line for now, for plotting the legend
                             theme_minimal()+
                             scale_color_manual(values=IsoPalette)+
                             scale_shape_manual(values=c(21, 22, 23, 24, 25, 12))+
                             theme(text = element_text(size = 20)))

## now save actual plot, with black outlines around mean values (if we kept the legend, it doesn't plot the colors clearly)
isoplot<-ggplot(isoplot_data)+
  geom_point(aes(d13C, d15N, fill=Description, color=Description, shape=`Sample Type`), alpha=0.5, size=3)+
  geom_errorbar(aes(x=meand13C, ymin=meand15N-sd15N, ymax=meand15N+sd15N, color=Description), size=0.5)+
  geom_errorbarh(aes(y=meand15N, xmin=meand13C-sd13C, xmax=meand13C+sd13C, color=Description), size=0.5)+
  geom_point(aes(meand13C, meand15N, fill=Description, color=Description, shape=`Sample Type`), color="black",size=5)+
  theme_minimal()+
  scale_color_manual(values=IsoPalette)+
  scale_fill_manual(values=IsoPalette)+
  scale_shape_manual(values=c(21, 22, 23, 24, 25, 12))+
  theme(legend.position="none")+
  ylab(expression(paste(delta,""^15,"N (‰)")))+
  xlab(expression(paste(delta,""^13,"C (‰)")))+
  panel_border(color = "grey85", size = 1, linetype = 1)+
  theme(text = element_text(size = 20))+
  theme(axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))+
  panel_border(color="black")

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

library(tidyverse)
library(ggplot2)
library(data.table)
library(readxl)
library(writexl)
library(ggsci) #has color palettes for various journals
library(dunn.test)

rawAA <- read_excel("../Sesarma_soils_manuscript/Data/Sapelo_AA_data.xlsx")

#split the unique_id column into "Creek type", "Status", "Zone", "Rep", "Soil Type"
data<-mutate(rawAA, sample=unique_id) %>% #still retain a column with sampleid
  separate(unique_id, c("Location", "Site", "Status", 
                        "Zone", "Rep", "SoilType"), sep="-") %>% 
  mutate(Zone = if_else(is.na(Zone), "Gut", Zone)) #added this to remove NAs to prevent issues downstream

data<-data[-1,]

#the crab samples don't have a zone, rep, soiltype, update these fields
data<-mutate(data, SampleType=case_when(
  Location=="FG" | Location=="HG" ~ "Sesarma", 
  Location=="LS" ~ "Sediment"
))

#reorder the Zones so they make sense. In order of creek to marsh
data$Zone<-ordered(data$Zone, c("MC", "MF", "HF", "MUD", "CH", "LP", "HP", "Gut"))

######### calcuate the amino acid degradation index according to Dauwe et al., 1999 and Philben et al., 2016
### sum (vari-avgi)/sdi * PC1i
## var = relative molar percentage of the amino acid, i
## avg = average of amino acid, i, for whole dataset
## sd = standard deviation of amino acid, i, for the whole dataset

#don't include Hyp, DAPA, MurA, GABA because they are source-specific
#their mol percents are already NA because they are not used to calculate the mole percentage

DI_data<-select(data, sample, starts_with("mole")) %>% 
  select(-mole.percent_Hyp, -mole.percent_MurA, #omit Hyp and DAPA (per Philben et al. 2016)
         -mole.percent_DAPA, -mole.percent_GABA, -mole.percent_NA) %>% #have var, the relative molar percentage of each amino acid
  mutate(across(starts_with("mole"), list(mean=~mean(.x, na.rm=TRUE), 
                                          sd=~sd(.x, na.rm=TRUE)))) #get the average and standard deviation of each amino acid for the whole dataset

# to do this, need to restructure the data

mean<-select(DI_data, sample, ends_with("mean"))
mean_reshaped<-pivot_longer(mean, -sample, names_to="comp", values_to="mean")
mean_reshaped$comp<-gsub("mole.percent_", "", mean_reshaped$comp)
mean_reshaped$comp<-gsub("_mean", "", mean_reshaped$comp)


sd<-select(DI_data, sample, ends_with("sd"))
sd_reshaped<-pivot_longer(sd, -sample, names_to="comp", values_to="sd")
sd_reshaped$comp<-gsub("mole.percent_", "", sd_reshaped$comp)
sd_reshaped$comp<-gsub("_sd", "", sd_reshaped$comp)

mole<-select(DI_data, 1:16)
mole_reshaped<-pivot_longer(mole, -sample, names_to="comp", values_to="mole")
mole_reshaped$comp<-gsub("mole.percent_", "", mole_reshaped$comp)
mole_reshaped$comp<-gsub("_mole", "", mole_reshaped$comp)

## use the mole data frame for the PCA
## convert NAs to 0 
mole[is.na(mole)] <- 0


## use the brackets -1 to remove the first column (sampleIDs)
pca<-princomp(mole[,-1])
biplot(pca)
pca_summ<-summary(pca, loadings=TRUE)
pca_summ_noGut<-summary(pca, loadings=TRUE)


#want to extract loading scores for each amino acid from the pca output
#likely not best way to do this, but works for now
loadings<-loadings(pca_summ) 
comp<-row.names(loadings)
col1<-loadings[1:15]
PC1_loadings<-as.data.frame(cbind(comp, col1))
PC1_loadings$col1<-as.numeric(PC1_loadings$col1)
PC1_loadings$comp<-gsub("mole.percent_", "", PC1_loadings$comp)

## join by sample and comp
DI<-full_join(mean_reshaped, sd_reshaped, by=c("sample", "comp")) %>% 
  full_join(mole_reshaped, by=c("sample", "comp")) %>% 
  full_join(PC1_loadings, by="comp") %>%  #col1 is the PC1 loading score for each comp
  dplyr::group_by(sample, comp) %>% 
  dplyr::mutate(a=((mole-mean)/sd)*col1) %>% #then get ( var - avg )/sd *PC1 loading for each sample, comp
  ungroup() %>% 
  select(sample, comp, a, col1) %>% 
  dplyr::group_by(sample) %>% #plyr was interfering with dplyr grouping, dplyr:: added to resolve
  dplyr::mutate(DI=sum(a, na.rm = TRUE)) #%>%  #see above re: dplyr::
# select(sample, DI) %>% 
# unique()

### join the DI with the rest of the dataset

data_DI<-full_join(data, DI,by="sample")

#first restructure the dataset 
reformat_data_DI<-data_DI %>% 
  select(sample, Location, Site, Status, Zone, Rep, SoilType, SampleType, DI, THAA_umgc, starts_with("umgc")) %>% 
  pivot_longer(cols=starts_with("umgc"), names_to="comp", 
               values_to="umgc")
#remove the extra "umgc_" from comp IDs
reformat_data_DI$comp<-gsub("umgc_", "", reformat_data_DI$comp)

#just get Hyp, MurA, and DAPA to plot bacterial and vegetation markers
sourcemarkers<-filter(reformat_data_DI, comp=="Hyp"|comp=="DAPA"|comp=="MurA") %>% 
  mutate(Type=case_when(Location=="FG" ~ "Fore-gut", 
                        Location=="HG" ~ "Hind-gut", 
                        Status=="G" ~ "Grazed", 
                        Status=="H" ~ "Un-grazed"), 
         Position=case_when(Zone=="MC" ~ "1", 
                            Zone=="MF" ~ "2", 
                            Zone=="HF" ~ "3", 
                            Zone=="MUD" ~ "4", 
                            Zone=="CH" ~ "5", 
                            Zone=="LP" ~ "6", 
                            Zone=="HP" ~ "7", 
                            Location=="FG" ~ "Fore-gut", 
                            Location=="HG" ~ "Hind-gut"), 
         SoilType=case_when(SoilType=="M" ~ "Upper Horizon", 
                            SoilType=="P" ~ "Lower Horizon", 
                            SoilType=="O" ~ "No Horizon", 
                            Location=="FG" ~ "Fore-gut", 
                            Location=="HG" ~ "Hind-gut"), 
         Biodeposit=case_when(SoilType=="Upper Horizon" ~ "Biodeposit",
                              SoilType=="Lower Horizon" ~ "Sediment", 
                              SoilType=="No Horizon" ~ "Sediment", 
                              SampleType=="Sesarma" ~ "Gut", 
                              SampleType=="Sesarma" ~ "Gut")) %>% 
  select(SoilType, Biodeposit, DI, THAA_umgc, comp, umgc, Type, Position) %>% 
  unique()

#rename DAPA, MurA and Hyp to more intutive names for plotting
sourcemarkers$comp<-gsub("DAPA", "Bacterial Marker (DAPA)", sourcemarkers$comp)
sourcemarkers$comp<-gsub("MurA", "Bacterial Marker (MurA)", sourcemarkers$comp)
sourcemarkers$comp<-gsub("Hyp", "Vegetation Marker (Hyp)", sourcemarkers$comp)

#rearrange the order of type, so they show up on the legend in a more intuitive way
sourcemarkers$Type<-ordered(sourcemarkers$Type, c("Grazed", "Un-grazed", "Fore-gut", "Hind-gut"))
sourcemarkers$SoilType<-ordered(sourcemarkers$SoilType, c("No Horizon", "Upper Horizon", "Lower Horizon", "Fore-gut", "Hind-gut"))
sourcemarkers$Biodeposit<-ordered(sourcemarkers$Biodeposit, c("Biodeposit", "Sediment", "Gut"))

##get mean value when possible, also get n to add to the plot
sourcemarkers_ordered<-dplyr::group_by(sourcemarkers, 
                                       Type, comp, Biodeposit) %>% 
  dplyr::mutate(mean_umgc=mean(umgc, na.rm=TRUE), ## have to add dplyr:: again - plyr causes issues
                n_samples=n())

sourcemarkers_ordered$Biodeposit<-gsub("Biodeposit", "USS", sourcemarkers_ordered$Biodeposit)

ggplot(sourcemarkers_ordered, aes(Type, mean_umgc, fill=Type))+
  geom_bar(stat="identity", position="dodge", color="black")+
  geom_point(aes(Type, umgc), position=position_dodge(width=0.78), alpha=0.2, size=1.5)+
  theme_minimal()+
  ylab(expression(paste(mu,"mol gC"^-1)))+
  scale_fill_npg()+
  scale_color_npg()+
  xlab("Sample Type")+
  theme(strip.text.x = element_text(size = 10))+
  panel_border(color="black")+
  facet_grid(comp~Biodeposit, scales="free_x")+ 
  guides(fill=guide_legend(title="Sample Type"))+
  ggtitle("Vegetation and Bacterial Biomarkers")+
  theme(strip.text.x = element_text(size = 12))


############ Run simple plot to see if there are any other amino acid differences
renamed_fortotalAA_plot<-reformat_data_DI %>% #also create Status and Position to match names on plots
  mutate(Type=case_when(Location=="FG" ~ "Fore-gut", 
                        Location=="HG" ~ "Hind-gut", 
                        Status=="G" ~ "Grazed", 
                        Status=="H" ~ "Un-grazed"), 
         Position=case_when(Zone=="MC" ~ "1", 
                            Zone=="MF" ~ "2", 
                            Zone=="HF" ~ "3", 
                            Zone=="MUD" ~ "4", 
                            Zone=="CH" ~ "5", 
                            Zone=="LP" ~ "6", 
                            Zone=="HP" ~ "7", 
                            Location=="FG" ~ "Fore-gut", 
                            Location=="HG" ~ "Hind-gut"), 
         SoilType=case_when(SoilType=="M" ~ "Upper Horizon",
                            SoilType=="P" ~ "Lower Horizon", 
                            SoilType=="O" ~ "No Horizon", 
                            Location=="FG" ~ "Fore-gut", 
                            Location=="HG" ~ "Hind-gut"), 
         Biodeposit=case_when(SoilType=="Upper Horizon" ~ "Biodeposit",
                              SoilType=="Lower Horizon" ~ "Sediment", 
                              SoilType=="No Horizon" ~ "Sediment", 
                              Location=="FG" ~ "Fore-gut", 
                              Location=="HG" ~ "Hind-gut"))

renamed_fortotalAA_plot$Type<-ordered(renamed_fortotalAA_plot$Type, c("Grazed", "Un-grazed", "Fore-gut", "Hind-gut"))
renamed_fortotalAA_plot$SoilType<-ordered(renamed_fortotalAA_plot$SoilType, c("No Horizon", "Upper Horizon", "Lower Horizon", "Fore-gut", "Hind-gut"))
renamed_fortotalAA_plot$Biodeposit<-ordered(renamed_fortotalAA_plot$Biodeposit, c("Biodeposit", "Sediment", "Fore-gut", "Hind-gut"))

renamed_fortotalAA_plot<-filter(renamed_fortotalAA_plot, comp!="NA")

renamed_fortotalAA_plot$Biodeposit<-gsub("Biodeposit", "USS", renamed_fortotalAA_plot$Biodeposit)
ggplot(renamed_fortotalAA_plot, aes(Biodeposit, umgc, fill=Type))+
  geom_boxplot(position="dodge")+
  #geom_point(aes(Biodeposit, umgc), position=position_dodge(width=0.78), alpha=0.2, size=1.5)+
  facet_wrap(.~comp, scales="free", ncol=4)+
  theme_minimal()+
  scale_fill_npg()+
  ylab(expression(paste(mu,"mol gC"^-1)))+
  xlab("")+
  guides(fill=guide_legend(title="Sample Type"))+
  ggtitle("")+
  theme(strip.text.x = element_text(size = 10))+
  panel_border(color="black")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.ticks=element_line(linewidth=1), 
        axis.ticks.length=unit(0.25, "cm"))
  
