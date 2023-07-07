##Hackerott et al 2023 Coral Phenotypic Plasticity and DNA Methylation Script
##Updated June 2023 by Serena Hackerott

#-----*Set Up-----------------------------------------------------------

setwd("C:/Users/snhac/OneDrive/Documents/FIU/Bonaire Project/Bonaire Year 1 Manuscript/Molecular Ecology/Revised Version/Final Resubmit Data and Scripts/Final Final")

####Load Packages####
library(adegenet)
library(car)
library(ComplexHeatmap)
library(cowplot)
library(cluster)
library(effectsize)
library(emmeans) 
library(geodata)
library(ggplot2)
library(ggpubr)
library(lattice)
library(lubridate)
library(msap)
library(pairwiseAdonis)
library(plotrix)
library(reshape2)
library(Rmisc)
library(scales)
library(vegan)
library(tidyr)


####Group Colors####
##Using Viridis palettes https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html 

##Genotypes
show_col(viridis_pal(option = "plasma")(30))
Genotype.colors.o<-c("#F5E926FF", "#FDB130FF", "#EB7556FF", "#D7566CFF", "#B52F8CFF", "#6E00A8FF", "#240691FF")
#All Genotypes: AC8, 10, 12, 15, AP1, 8, 10

AC.genotype.colors.o <- c("#F5E926FF", "#FDB130FF", "#EB7556FF", "#D7566CFF")
#AC Genotypes: AC8, 10, 12, 15

AP.genotype.colors.o <- c("#B52F8CFF", "#6E00A8FF", "#240691FF")
#AP Genotypes: AP1, 8, 10


##Sites
show_col(viridis_pal(option = "viridis")(20))
Site.colors.o <- c("#482677FF", "#1F968BFF", "#32648EFF", "#56C667FF")
#Sites: BD, SS, KL, OL

AC.site.colors.o <- c("#482677FF", "#1F968BFF", "#32648EFF")
#AC Sites: BD, SS, KL

AP.site.colors.o <- c("#482677FF", "#1F968BFF", "#56C667FF")
#AP Sites: BD, SS, OL


##Seasons
show_col(viridis_pal(option = "turbo")(20))
Season.colors.o <- c("#FD8A26FF", "#455ED2FF", "#EFCD3AFF", "#A51301FF")
#Season: Fall, Winter, Spring, Summer

Period.colors.o <- c("#455ED2FF", "#A51301FF")
#Seasonal Periods: Cooling: (Fall and Winter), Warming: (Spring and Summer)


##Additional
show_col(viridis_pal(option = "magma")(20))
Species.colors.o<-c("#E85362FF", "#802582FF")
#Species: Acer, Apal

show_col(viridis_pal(option = "turbo")(20))
GE.colors.o<-c("#FD8A26FF", "#18DEC1FF", "#4681F7FF")
#Factors: Genotype, Site, Season

show_col(viridis_pal(option = "mako")(20))
Met.colors.o<-c("#DEF5E5FF", "#5ACCADFF", "#3482A4FF", "#3C5397FF")
#Methylation Levels: "NMT", "HMM", "ICM", "HPM/NA"


####Plot Feature Sizes####
axis.title.sz=18
axis.txt.sz=14
leg.title.sz=15
leg.txt.sz=12
levels.sz=7
sig.sz=5
panel.lab.sz=25


#-----*Study Site Map-----------------------------------------------------------

####**Map Part 1: Bonaire Map and Study Sites####

##Load Bonaire Map data
BES_remote <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_BES_1_sf.rds"

BES_rds <- file.path(tempdir(), "gadm36_BES_1_sf.rds")

if (toupper(Sys.info()["sysname"]) == "WINDOWS") {
  download.file(
    url = BES_remote,
    destfile = BES_rds,
    method = "wininet",
    mode = "wb"
  )
} else {
  download.file(
    url = BES_remote,
    destfile = BES_rds,
    method = "auto"
  )
}


##Read
BES.sf<-readRDS(BES_rds)

##Subset Bonaire
Bonaire.sf<-BES.sf[1,]

##Load study site data
Bonaire.Sites<-read.csv("BonaireSites.csv", header=TRUE)
str(Bonaire.Sites)

##Factors
Bonaire.Sites$Site<-factor(Bonaire.Sites$Site, levels=c("BD", "SS", "KL", "OL", "KR"), ordered=TRUE)
Bonaire.Sites$Type<-factor(Bonaire.Sites$Type, levels=c("Study", "Capital"), ordered=TRUE)


####Plot map of Bonaire with study sites####
Study_Sites.plot<-ggplot(Bonaire.sf)+
  geom_sf(fill="white")+
  labs(x="Longitude", y="Latitude")+
  geom_point(data=Bonaire.Sites, 
             aes(x=Long, y=Lat, colour=Site, shape=Type), size=6.5)+
  scale_colour_manual(values=c(Site.colors.o, "darkred"))+
  scale_shape_manual(values=c(16, 17))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz))

Study_Sites.plot


#-----*Coral Genotype Selection Background Information-----------------------------------------------------------

####**Coral Genotype Selection Part 1: Outplant Survival####
##Note- Data collected by Reef Renewal Foundation Bonaire##

##Load Data
Outplant<-read.csv("RRFB_Outplant.csv", header=TRUE)
str(Outplant)

##Factors
Outplant$Species<-factor(Outplant$Species, levels=c("Acer", "Apal"), ordered=TRUE)
Outplant$Timepoint<-factor(Outplant$Timepoint, levels=c("Initial", "Mid", "Year"), ordered=TRUE)
Outplant$Genotype<-factor(Outplant$Genotype)
Outplant$Site<-factor(Outplant$Site)

##Live Tissue bin midpoints
Tissue.Bins<-c(NA, ((.25+.01)/2), ((.50+.26)/2), ((.75+.51)/2), ((.99+.76)/2), 1)

##Calculate weighted average % Live Tissue per Cluster
Outplant$LiveTissue<-(((Outplant$Tissue_1_25*Tissue.Bins[2]+Outplant$Tissue_26_50*Tissue.Bins[3]+
                          Outplant$Tissue_51_75*Tissue.Bins[4]+Outplant$Tissue_76_99*Tissue.Bins[5]+
                          Outplant$Tissue_100*Tissue.Bins[6])/Outplant$nPresent)*100)

##Average Live Tissue across Clusters, by Species, Site, Month, Timepoint, Genotype
names(Outplant)
Outplant.avg<-summarySE(Outplant, measurevar="LiveTissue", groupvars=c("Species", "Site", "Month", "Timepoint", "Genotype"), na.rm=TRUE)

##Add Site.Genotype grouping variable 
Outplant.avg$Site.Geno<-paste(Outplant.avg$Site, Outplant.avg$Genotype, sep=".")


####Outplant Performance over Time and Site- Acropora cervicornis####

##Subset data for A. cervicornis 
Outplant.Acer<-subset(Outplant, Species=="Acer")
Outplant.avg.Acer<-subset(Outplant.avg, Species=="Acer")

##Filter data for genotypes of interest (AC8, AC10, AC12, AC15)
Outplant.Acer.Study<-Outplant.Acer[which(Outplant.Acer$Genotype=="AC8" | Outplant.Acer$Genotype=="AC10" |
                                           Outplant.Acer$Genotype=="AC12" | Outplant.Acer$Genotype=="AC15"),]
Outplant.Acer.Study$Genotype<-factor(Outplant.Acer.Study$Genotype, 
                                     levels=c("AC8", "AC10", "AC12", "AC15"), ordered=TRUE)

Outplant.avg.Acer.Study<-Outplant.avg.Acer[which(Outplant.avg.Acer$Genotype=="AC8" | Outplant.avg.Acer$Genotype=="AC10" | 
                                                   Outplant.avg.Acer$Genotype=="AC12" | Outplant.avg.Acer$Genotype=="AC15"),]
Outplant.avg.Acer.Study$Genotype<-factor(Outplant.avg.Acer.Study$Genotype, 
                                         levels=c("AC8", "AC10", "AC12", "AC15"), ordered=TRUE)

##Colonies per cluster
mean(Outplant.Acer.Study$nCorals)
std.error(Outplant.Acer.Study$nCorals)
#10.23 +/- 0.21

##Number of clusters per site and genotype
Outplant.Acer.Study.Clusters<-data.frame(count(Outplant.Acer.Study, vars=c("Site", "Genotype")))

mean(Outplant.Acer.Study.Clusters$freq)
std.error(Outplant.Acer.Study.Clusters$freq)
#11.5 +/-1.85

##Compare Slopes of Tissue over Time between Genotypes and Sites
Acer.Tissue<-lm(LiveTissue~Month + Month:Genotype+ Month:Site + Month:Genotype:Site, data=Outplant.Acer.Study)
summary(Acer.Tissue)
anova(Acer.Tissue)
#                       Df  Sum Sq Mean Sq F value   Pr(>F)    
# Month                 1   165.5  165.48  1.5599 0.214490    
# Month:Genotype        3  4652.5 1550.82 14.6180 5.17e-08 ***
# Month:Site            2  2861.5 1430.75 13.4862 6.20e-06 ***
# Month:Genotype:Site   4  1581.7  395.41  3.7272 0.007093 ** 
# Residuals           104 11033.3  106.09  


####Plot Live Tissue over Time by Genotypes and Sites- A. cervicornis####
Outplant_Geno_AC.plot<-ggplot(Outplant.avg.Acer.Study, aes(x=Month, y=LiveTissue, colour=Genotype, group=Site.Geno)) + 
  geom_errorbar(aes(ymin=LiveTissue-se, ymax=LiveTissue+se), width=1, position=position_dodge(0.6), size=0.7) +
  geom_line(position=position_dodge(0.6), size=1.5, aes(linetype=Site)) +
  geom_point(position=position_dodge(0.6), size=4)+
  scale_colour_manual(values=AC.genotype.colors.o)+
  scale_linetype_manual(values=c("twodash", "solid", "dotted"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.position = c(0.9, 0.28), legend.box.background = element_rect(color = "black"))+
  labs(x="Months after Outplanting", y="Average Percent Live Tissue")+
  ylim(50,100)+
  xlim(0,19)+
  annotate("text", x = 0.1, y = 56, label = "Month x Genotype p < 0.001", size=sig.sz, hjust = 0, fontface="bold.italic")+
  annotate("text", x = 0.1, y = 53, label = "Month x Site p < 0.001", size=sig.sz, hjust = 0, fontface="bold.italic")+
  annotate("text", x = 0.1, y = 50, label = "Month x Genotype x Site p = 0.007", size=sig.sz, hjust = 0, fontface="bold.italic")

Outplant_Geno_AC.plot


####Outplant Performance over Time and Site- Acropora palmata####

##Subset data for A. palmata 
Outplant.Apal<-subset(Outplant, Species=="Apal")
Outplant.avg.Apal<-subset(Outplant.avg, Species=="Apal")

##Filter data for genotypes of interest (AP1, AP8, AP10)
Outplant.Apal.Study<-Outplant.Apal[which(Outplant.Apal$Genotype=="AP1" | Outplant.Apal$Genotype=="AP8" | 
                                           Outplant.Apal$Genotype=="AP10"),]
Outplant.Apal.Study$Genotype<-factor(Outplant.Apal.Study$Genotype, 
                                     levels=c("AP1", "AP8", "AP10"), ordered=TRUE)

Outplant.avg.Apal.Study<-Outplant.avg.Apal[which(Outplant.avg.Apal$Genotype=="AP1" | Outplant.avg.Apal$Genotype=="AP8" | 
                                                   Outplant.avg.Apal$Genotype=="AP10"),]
Outplant.avg.Apal.Study$Genotype<-factor(Outplant.avg.Apal.Study$Genotype, 
                                         levels=c("AP1", "AP8", "AP10"), ordered=TRUE)

##Colonies per cluster
mean(Outplant.Apal.Study$nCorals)
std.error(Outplant.Apal.Study$nCorals)
#5.02 +/- 0.02

##Number of clusters per site and genotype
Outplant.Apal.Study.Clusters<-data.frame(count(Outplant.Apal.Study, vars=c("Site", "Genotype")))

mean(Outplant.Apal.Study.Clusters$freq)
std.error(Outplant.Apal.Study.Clusters$freq)
# 8.43 +/-1.32

##Compare Slopes of Tissue over Time between Genotypes and Sites
Apal.Tissue<-lm(LiveTissue~Month + Month:Genotype+ Month:Site + Month:Genotype:Site, data=Outplant.Apal.Study)
summary(Apal.Tissue)
anova(Apal.Tissue)
#                     Df  Sum Sq Mean Sq F value   Pr(>F)    
# Month                1   970.5   970.5  2.3129 0.134477    
# Month:Genotype       2  7856.8  3928.4  9.3627 0.000344 ***
# Month:Site           2  4871.0  2435.5  5.8046 0.005355 ** 
# Month:Genotype:Site  2  1292.2   646.1  1.5399 0.224207    
# Residuals           51 21398.4   419.6 


####Plot Live Tissue over Time by Genotypes and Sites- A. palmata####
Outplant_Geno_AP.plot<-ggplot(Outplant.avg.Apal.Study, aes(x=Month, y=LiveTissue, colour=Genotype, group=Site.Geno)) + 
  geom_errorbar(aes(ymin=LiveTissue-se, ymax=LiveTissue+se), width=1, position=position_dodge(0.6), size=0.7) +
  geom_line(position=position_dodge(0.6), size=1.5, aes(linetype=Site)) +
  geom_point(position=position_dodge(0.6), size=4)+
  scale_colour_manual(values=AP.genotype.colors.o)+
  scale_linetype_manual(values=c("dotted", "solid", "twodash"))+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.position = c(0.9, 0.26), legend.box.background = element_rect(color = "black"))+
  labs(x="Months after Outplanting", y="Average Percent Live Tissue")+
  ylim(-10,100)+
  xlim(0,19)+
  annotate("text", x = 0.1, y = 2, label = "Month x Genotype p < 0.001", size=sig.sz, hjust = 0, fontface="bold.italic")+
  annotate("text", x = 0.1, y = -4, label = "Month x Site p = 0.005", size=sig.sz, hjust = 0, fontface="bold.italic")+
  annotate("text", x = 0.1, y = -10, label = "Month x Genotype x Site p = 0.224", size=sig.sz, hjust = 0, fontface="italic")

Outplant_Geno_AP.plot

####***Figure S1: Plot Map of Study Sites and Outplant Performance for Genotype Selection####

##Figure S1 Study Sites Map and Live Tissue over Time by Genotypes and Sites

Outplant_fig<-plot_grid(Outplant_Geno_AC.plot, Outplant_Geno_AP.plot, 
                        rel_widths=c(1, 1), rel_heights=c(1, 1), 
                        nrow=2, ncol=1, byrow=T, labels = c('B', 'C'), 
                        label_size=panel.lab.sz, align="v")

Study_fig<-plot_grid(Study_Sites.plot, Outplant_fig, 
                     rel_widths=c(1, 0.8), rel_heights=c(1, 1), 
                     nrow=1, ncol=2, byrow=T, labels = c('A', ''), 
                     label_size=panel.lab.sz, align="h")

ggsave(filename="FigS1_Study_Sites_Genotype_Selection.png", plot=Study_fig, dpi=300, width=18, height=12, units="in")
ggsave(filename="FigS1_Study_Sites_Genotype_Selection.pdf", plot=Study_fig, dpi=300, width=18, height=12, units="in")


#-----*Environmental Conditions---------------------------------------------------


####**Environment Part 1: Seasonal Sampling####


####Temperature and Conductivity####

##Load Data
T.S<-read.csv("Temp_Cond.csv", header=TRUE)
str(T.S)

##Remove M13
T.S<-T.S[-c(which(T.S$Month.n=="M13")),]

##Factors
T.S$Site<-factor(T.S$Site, levels=c("BD", "SS", "KL", "OL"), ordered=TRUE)
T.S$Month.a<-factor(T.S$Month.a, levels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"), ordered=TRUE)
T.S$Month.n<-factor(T.S$Month.n, levels=c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13"), ordered=TRUE)

##Clean Date Data
T.S$Month <- sprintf("%02d",T.S$Month) # fix to 2 characters
T.S$Day <- sprintf("%02d",T.S$Day) # fix to 2 characters
Date2<-paste(T.S$Year, T.S$Month, T.S$Day, sep="/")
T.S$Date<-c(as.Date(Date2, format='%Y/%m/%d'))
T.S$Date.n<-paste(T.S$Month.n, T.S$Day, sep=".")

##Clean Date Time Data
T.S$Dtime<-ymd_hms(c(paste(as.character(T.S$Date), as.character(T.S$Time))), tz="America/Curacao")

##Add Column for Season
T.S$Season<-c(rep("Fall", nrow(T.S)))
T.S$Season[which(T.S$Month.a=="Dec" | T.S$Month.a=="Jan" | T.S$Month.a=="Feb")]<-"Winter"
T.S$Season[which(T.S$Month.a=="Mar" | T.S$Month.a=="Apr" | T.S$Month.a=="May")]<-"Spring"
T.S$Season[which(T.S$Month.a=="Jun" | T.S$Month.a=="Jul" | T.S$Month.a=="Aug")]<-"Summer"
T.S$Season<-factor(T.S$Season, levels=c("Fall", "Winter", "Spring", "Summer"), ordered=TRUE)

##Add column for Site and Season
T.S$Site.Sea<-paste(T.S$Site, T.S$Season, sep=".")
T.S$Site.Sea<-factor(T.S$Site.Sea, levels=c("BD.Fall", "SS.Fall", "KL.Fall", "OL.Fall",
                                            "BD.Winter", "SS.Winter", "KL.Winter", "OL.Winter",
                                            "BD.Spring", "SS.Spring", "KL.Spring", "OL.Spring",
                                            "BD.Summer", "SS.Summer", "KL.Summer", "OL.Summer"), ordered=TRUE)

####Daily Temperature and Conductivity Values####

##Calculate Daily Mean, Min, Max, and Standard Deviation
TS.mean<-aggregate(T.S[8:9], list(T.S$Site, T.S$Date, T.S$Month.n, T.S$Month.a, T.S$Date.n, T.S$Season, T.S$Site.Sea), mean)
names(TS.mean)<-c("Site", "Date", "Month.n","Month.a", "Date.n", "Season", "Site.Sea", "Temp.a", "Cond.a")

TS.min<-aggregate(T.S[8:9], list(T.S$Site, T.S$Date, T.S$Month.n, T.S$Month.a, T.S$Date.n, T.S$Season, T.S$Site.Sea), min)
names(TS.min)<-c("Site", "Date", "Month.n", "Month.a", "Date.n", "Season", "Site.Sea", "Temp.min", "Cond.min")

TS.max<-aggregate(T.S[8:9], list(T.S$Site, T.S$Date, T.S$Month.n, T.S$Month.a, T.S$Date.n, T.S$Season, T.S$Site.Sea), max)
names(TS.max)<-c("Site", "Date", "Month.n", "Month.a", "Date.n", "Season", "Site.Sea",  "Temp.max", "Cond.max")

TS.sd<-aggregate(T.S[8:9], list(T.S$Site, T.S$Date, T.S$Month.n, T.S$Month.a, T.S$Date.n, T.S$Season, T.S$Site.Sea), sd)
names(TS.sd)<-c("Site", "Date", "Month.n", "Month.a", "Date.n", "Season", "Site.Sea", "Temp.sd", "Cond.sd")

##Dataframe of Summary Statistics
TS.sum<- merge(TS.mean, TS.min)
TS.sum<- merge(TS.sum, TS.max)
TS.sum<- merge(TS.sum, TS.sd)
rownames(TS.sum)<-paste(TS.sum$Site, TS.sum$Date.n, sep=".")


####Dates Complete for All and AC Only Sites####

##Find Dates Complete in Each Sites
TS.sum.rm<-na.omit(TS.sum)
names(TS.sum.rm)
Date.BD<-subset(TS.sum.rm[,c(2,5)], TS.sum.rm$Site=="BD")
Date.SS<-subset(TS.sum.rm[,c(2,5)], TS.sum.rm$Site=="SS")
Date.KL<-subset(TS.sum.rm[,c(2,5)], TS.sum.rm$Site=="KL")
Date.OL<-subset(TS.sum.rm[,c(2,5)], TS.sum.rm$Site=="OL")

##Dates Complete for All Sites
Date.C<-merge(Date.OL, Date.BD, all.x=FALSE, all.y=FALSE)
Date.C<-merge(Date.C, Date.SS, all.x=FALSE, all.y=FALSE)
Date.C<-merge(Date.C, Date.KL, all.x=FALSE, all.y=FALSE)
length(unique(Date.C$Date))
#80 

##Dates Complete for AC Sites
Date.AC<-merge(Date.BD, Date.SS, all.x=FALSE, all.y=FALSE)
Date.AC<-merge(Date.AC, Date.KL, all.x=FALSE, all.y=FALSE)
length(unique(Date.AC$Date))
#196 

##Dates Complete for AP Sites
Date.AP<-merge(Date.OL, Date.BD, all.x=FALSE, all.y=FALSE)
Date.AP<-merge(Date.AP, Date.OL, all.x=FALSE, all.y=FALSE)
length(unique(Date.AP$Date))
#80 

##Daily values with Complete Dates
TS.sum.C<-TS.sum.rm[TS.sum.rm$Date %in% Date.C$Date,]

##AC Sites
TS.sum.AC<-TS.sum.rm[-c(TS.sum.rm$Site=="OL"),]
TS.sum.AC$Site<-factor(TS.sum.AC$Site, levels=c("BD", "SS", "KL"), ordered=TRUE)
TS.sum.C.AC<-TS.sum.AC[TS.sum.AC$Date %in% Date.AC$Date,]

##AP Sites
TS.sum.AP<-TS.sum.rm[-c(TS.sum.rm$Site=="KL"),]
TS.sum.AP$Site<-factor(TS.sum.AP$Site, levels=c("BD", "SS", "OL"), ordered=TRUE)
TS.sum.C.AP<-TS.sum.AP[TS.sum.AP$Date %in% Date.AP$Date,]


####Seasonal Sampling####

####Sampling Period Dates
SamP<-read.csv("Samples.csv", header=TRUE)

##Factors
SamP$Period<-factor(SamP$Period, levels=c("Fall", "Winter", "Spring", "Summer"), ordered=TRUE)

##Clean Date Data
SamP$Month <- sprintf("%02d",SamP$Month) # fix to 2 characters
SamP$Day <- sprintf("%02d",SamP$Day) # fix to 2 characters
SamP$Date<-c(as.Date(paste(SamP$Year, SamP$Month, SamP$Day, sep="/"), format='%Y/%m/%d'))

##Add Date as character for plotting
SamP$Date.c<-ymd(as.character(SamP$Date), tz="America/Curacao")


####NOAA SST Data for Missing Temperature
NOAA_T<-read.csv("NOAA.ABC_T_fill.csv", header=TRUE)

##Factors
NOAA_T$Month.a<-factor(NOAA_T$Month.a, levels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"), ordered=TRUE)
NOAA_T$Month.n<-factor(NOAA_T$Month.n, levels=c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"), ordered=TRUE)
NOAA_T$Site<-factor(NOAA_T$Site, levels=c("NOAA", "BD", "SS", "KL", "OL"), ordered=TRUE)

##Clean Date Data
NOAA_T$Month <- sprintf("%02d",NOAA_T$Month) # fix to 2 characters
NOAA_T$Day <- sprintf("%02d",NOAA_T$Day) # fix to 2 characters
NOAA_T$Date<-c(as.Date(paste(NOAA_T$Year, NOAA_T$Month, NOAA_T$Day, sep="/"), format='%Y/%m/%d'))
NOAA_T$Date.n<-paste(NOAA_T$Month.n, NOAA_T$Day, sep=".")

##Merge with Temp Summary Data (Site, Date, Temp.a, Temp.min, Temp.max)
names(TS.sum)
T.sum<-TS.sum[,c(1, 2, 8, 10, 12)]
str(T.sum)
T.sum$Site<-factor(T.sum$Site, levels=c("NOAA", "BD", "SS", "KL", "OL"), ordered=TRUE)

names(NOAA_T)
str(NOAA_T[,c(1, 10, 9, 7, 8)])
T.sum<-rbind(NOAA_T[,c(1, 10, 9, 7, 8)], T.sum)

##Add Date as character for plotting
T.sum$Date.c<-ymd(as.character(T.sum$Date), tz="America/Curacao")


####Plot Temperature by Site with Sampling Dates and NOAA Fill####
Temp.ribbon.plot<-ggplot() +
  geom_ribbon(data=T.sum, aes(x=Date.c, ymin=Temp.min, ymax=Temp.max, fill=Site), alpha=0.5) +
  geom_line(data=T.sum, aes(color=Site, x=Date.c, y=Temp.a), size=1.25) +
  scale_color_manual(values = c("gray80", Site.colors.o))+
  scale_fill_manual(values = c("gray80", Site.colors.o))+
  theme_classic()+
  ylim(25, 30.5)+
  labs(title=NULL,x=NULL, y = "Temperature (\u00B0C)")+
  theme(axis.text.x = element_text(size=axis.txt.sz, colour="black"), axis.text.y = element_text(size=axis.txt.sz, colour="black"), 
        axis.title.y = element_text(size = axis.title.sz), 
        legend.position=c(.5, 0.05), legend.direction = "horizontal", 
        legend.text=element_text(size=axis.title.sz), legend.title=element_text(size=axis.title.sz))+
  geom_point(data=SamP, aes(x=Date.c,y=c(29.5, 28, 27.75, 29)),color=c(Season.colors.o), size=5.5)+
  annotate("text", x = SamP$Date.c[1], y = 29.9, label = "  Fall", size=levels.sz, hjust=0)+
  annotate("text", x = SamP$Date.c[2], y = 28.4, label = "Winter", size=levels.sz, hjust=0)+
  annotate("text", x = SamP$Date.c[3], y = 28.4, label = "Spring", size=levels.sz)+
  annotate("text", x = SamP$Date.c[4], y = 29.4, label = "Summer", size=levels.sz, hjust=1)

Temp.ribbon.plot


####**Environment Part 2: Temperature####


####Daily Average Temperature by Site and Season####

####Dates complete at all Sites

##Check normality
hist(TS.sum.C$Temp.a)
qqnorm(TS.sum.C$Temp.a)
shapiro.test(TS.sum.C$Temp.a)
#Not normal

##Compare models
Temp.a.glm.gaus<-glm(Temp.a~Site*Season, data=TS.sum.C, family=gaussian)
Temp.a.glm.gam<-glm(Temp.a~Site*Season, data=TS.sum.C, family=Gamma)
Temp.a.glm.inv<-glm(Temp.a~Site*Season, data=TS.sum.C, family=inverse.gaussian)

AIC(Temp.a.glm.gaus, Temp.a.glm.gam, Temp.a.glm.inv)
##Inverse Gaussian distribution fits best

##Check residuals
plot(fitted(Temp.a.glm.inv), resid(Temp.a.glm.inv))
abline(0,0)

qqnorm(resid(Temp.a.glm.inv))
qqline(resid(Temp.a.glm.inv))

plot(density(resid(Temp.a.glm.inv)))

##Model results
summary(Temp.a.glm.inv)
Anova(Temp.a.glm.inv)
#             LR Chisq Df Pr(>Chisq)    
# Site            7.34  3    0.06175 .  
# Season        323.84  1    < 2e-16 ***
# Site:Season     0.34  3    0.95261   

emmeans(Temp.a.glm.inv, pairwise~Season)
# contrast       estimate       SE  df z.ratio p.value
# Fall - Winter -5.07e-05 2.84e-06 Inf -17.874  <.0001


####Dates complete at AC Sites

##Check normality
hist(TS.sum.C.AC$Temp.a)
qqnorm(TS.sum.C.AC$Temp.a)
shapiro.test(TS.sum.C.AC$Temp.a)
#Not normal

##Compare models
Temp.a.AC.glm.gaus<-glm(Temp.a~Site*Season, data=TS.sum.C.AC, family=gaussian)
Temp.a.AC.glm.gam<-glm(Temp.a~Site*Season, data=TS.sum.C.AC, family=Gamma)
Temp.a.AC.glm.inv<-glm(Temp.a~Site*Season, data=TS.sum.C.AC, family=inverse.gaussian)

AIC(Temp.a.AC.glm.gaus, Temp.a.AC.glm.gam, Temp.a.AC.glm.inv)
##Gaussian distribution fits best

##Check residuals
plot(fitted(Temp.a.AC.glm.gaus), resid(Temp.a.AC.glm.gaus))
abline(0,0)

qqnorm(resid(Temp.a.AC.glm.gaus))
qqline(resid(Temp.a.AC.glm.gaus))

plot(density(resid(Temp.a.AC.glm.gaus)))

##Model results
summary(Temp.a.AC.glm.gaus)
Anova(Temp.a.AC.glm.gaus)
#             LR Chisq Df Pr(>Chisq)    
# Site            2.68  2     0.2621    
# Season       1378.58  2     <2e-16 ***
# Site:Season     0.25  4     0.9928     

emmeans(Temp.a.AC.glm.gaus, pairwise~Season)
# contrast        estimate     SE  df z.ratio p.value
# Fall - Winter      1.717 0.0465 Inf  36.921  <.0001
# Fall - Summer      0.709 0.0476 Inf  14.895  <.0001
# Winter - Summer   -1.008 0.0473 Inf -21.329  <.0001


####Plot Average Temperature by Site and Season####
Temp.a.plot<-ggplot(TS.sum, aes(x=Site.Sea, y=Temp.a, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), axis.text.x = element_blank(), 
        axis.text.y = element_text(size=axis.txt.sz, colour="black"), legend.position="none" )+
  ylab("Daily Average Temperature (\u00B0C)")+
  xlab(NULL)+
  ylim(26, 30.2)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 30.2, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 30.2, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 30.2, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 30.2, label = "Summer", size=sig.sz)+
  annotate("text", x = 2.5, y = 30, label = "A", size=sig.sz, fontface="bold.italic")+
  annotate("text", x = 6.5, y = 30, label = "B", size=sig.sz, fontface="bold.italic")+
  annotate("text", x = 14.5, y = 30, label = "C", size=sig.sz, fontface="bold.italic")

Temp.a.plot

####Daily Temperature Standard Deviation by Site and Season####

####Dates complete at all Sites

##Check normality
hist(TS.sum.C$Temp.sd)
qqnorm(TS.sum.C$Temp.sd)
shapiro.test(TS.sum.C$Temp.sd)
#Not normal

##Compare models
Temp.sd.glm.gaus<-glm(Temp.sd~Site*Season, data=TS.sum.C, family=gaussian)
Temp.sd.glm.gam<-glm(Temp.sd~Site*Season, data=TS.sum.C, family=Gamma)
Temp.sd.glm.inv<-glm(Temp.sd~Site*Season, data=TS.sum.C, family=inverse.gaussian)

AIC(Temp.sd.glm.gaus, Temp.sd.glm.gam, Temp.sd.glm.inv)
##Inverse Gaussian distribution fits best

##Check residuals
plot(fitted(Temp.sd.glm.inv), resid(Temp.sd.glm.inv))
abline(0,0)

qqnorm(resid(Temp.sd.glm.inv))
qqline(resid(Temp.sd.glm.inv))

plot(density(resid(Temp.sd.glm.inv)))

##Model results
summary(Temp.sd.glm.inv)
Anova(Temp.sd.glm.inv)
#             LR Chisq Df Pr(>Chisq)    
# Site         28.3484  3  3.069e-06 ***
# Season        1.7913  1     0.1808    
# Site:Season   0.8870  3     0.8286    

emmeans(Temp.sd.glm.inv, pairwise~Site)
# contrast estimate   SE  df z.ratio p.value
# BD - SS     -2.41 18.8 Inf  -0.128  0.9992
# BD - KL     13.89 18.3 Inf   0.759  0.8728
# BD - OL    -73.41 23.5 Inf  -3.123  0.0097
# SS - KL     16.31 18.3 Inf   0.892  0.8089
# SS - OL    -71.00 23.5 Inf  -3.023  0.0133
# KL - OL    -87.30 23.1 Inf  -3.787  0.0009


####Dates complete at AC Sites

##Check normality
hist(TS.sum.C.AC$Temp.sd)
qqnorm(TS.sum.C.AC$Temp.sd)
shapiro.test(TS.sum.C.AC$Temp.sd)
#Not normal

##Compare models
Temp.sd.AC.glm.gaus<-glm(Temp.sd~Site*Season, data=TS.sum.C.AC, family=gaussian)
Temp.sd.AC.glm.gam<-glm(Temp.sd~Site*Season, data=TS.sum.C.AC, family=Gamma)
Temp.sd.AC.glm.inv<-glm(Temp.sd~Site*Season, data=TS.sum.C.AC, family=inverse.gaussian)

AIC(Temp.sd.AC.glm.gaus, Temp.sd.AC.glm.gam, Temp.sd.AC.glm.inv)
##Inverse gaussian distribution fits best

##Check residuals
plot(fitted(Temp.sd.AC.glm.inv), resid(Temp.sd.AC.glm.inv))
abline(0,0)

qqnorm(resid(Temp.sd.AC.glm.inv))
qqline(resid(Temp.sd.AC.glm.inv))

plot(density(resid(Temp.sd.AC.glm.inv)))

##Model results
summary(Temp.sd.AC.glm.inv)
Anova(Temp.sd.AC.glm.inv)
#             LR Chisq Df Pr(>Chisq)    
# Site          6.2479  2    0.04398 *  
# Season       20.2967  2  3.914e-05 ***
# Site:Season   1.9144  4    0.75150   

emmeans(Temp.sd.AC.glm.inv, pairwise~Season)
# contrast        estimate   SE  df z.ratio p.value
# Fall - Winter      -11.2 9.44 Inf  -1.181  0.4642
# Fall - Summer       25.7 8.57 Inf   2.999  0.0076
# Winter - Summer     36.9 8.83 Inf   4.172  0.0001

emmeans(Temp.sd.AC.glm.inv, pairwise~Site)
# contrast estimate   SE  df z.ratio p.value
# BD - SS      1.53 9.28 Inf   0.165  0.9851
# BD - KL     17.75 8.81 Inf   2.014  0.1087
# SS - KL     16.22 8.77 Inf   1.850  0.1536


####Plot Temperature Standard Deviation by Site and Season####
Temp.sd.plot<-ggplot(TS.sum, aes(x=Site.Sea, y=Temp.sd, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), legend.position="none", 
        axis.text.x = element_text(size=axis.txt.sz, colour="black"), axis.text.y = element_text(size=axis.txt.sz, colour="black"))+
  ylab("Daily Temperature St. Dev. (\u00B0C)")+
  xlab(NULL)+
  ylim(0, 0.5)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 0.5, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 0.5, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 0.5, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 0.5, label = "Summer", size=sig.sz)+
  annotate("text", x = 2.5, y = 0.475, label = "A", size=sig.sz, fontface="bold.italic")+
  annotate("text", x = 6.5, y = 0.475, label = "A", size=sig.sz, fontface="bold.italic")+
  annotate("text", x = 14.4, y = 0.475, label = "B", size=sig.sz, fontface="bold.italic")

Temp.sd.plot


####**Environment Part 2: Conductivity####

####Daily Average Conductivity by Site and Season####

####Dates complete at all Sites

##Check normality
hist(TS.sum.C$Cond.a)
qqnorm(TS.sum.C$Cond.a)
shapiro.test(TS.sum.C$Cond.a)
#Not normal

##Compare models
Cond.a.glm.gaus<-glm(Cond.a~Site*Season, data=TS.sum.C, family=gaussian)
Cond.a.glm.gam<-glm(Cond.a~Site*Season, data=TS.sum.C, family=Gamma)
Cond.a.glm.inv<-glm(Cond.a~Site*Season, data=TS.sum.C, family=inverse.gaussian)

AIC(Cond.a.glm.gaus, Cond.a.glm.gam, Cond.a.glm.inv)
##Gaussian distribution fits best

##Check residuals
plot(fitted(Cond.a.glm.gaus), resid(Cond.a.glm.gaus))
abline(0,0)

qqnorm(resid(Cond.a.glm.gaus))
qqline(resid(Cond.a.glm.gaus))

plot(density(resid(Cond.a.glm.gaus)))

##Model results
summary(Cond.a.glm.gaus)
Anova(Cond.a.glm.gaus)
#             LR Chisq Df Pr(>Chisq)    
# Site         211.624  3  < 2.2e-16 ***
# Season        57.190  1  3.956e-14 ***
# Site:Season   15.902  3   0.001188 ** 


####Dates complete at AC Sites

##Check normality
hist(TS.sum.C.AC$Cond.a)
qqnorm(TS.sum.C.AC$Cond.a)
shapiro.test(TS.sum.C.AC$Cond.a)
#Not normal

##Compare models
Cond.a.AC.glm.gaus<-glm(Cond.a~Site*Season, data=TS.sum.C.AC, family=gaussian)
Cond.a.AC.glm.gam<-glm(Cond.a~Site*Season, data=TS.sum.C.AC, family=Gamma)
Cond.a.AC.glm.inv<-glm(Cond.a~Site*Season, data=TS.sum.C.AC, family=inverse.gaussian)

AIC(Cond.a.AC.glm.gaus, Cond.a.AC.glm.gam, Cond.a.AC.glm.inv)
##Gaussian distribution fits best

##Check residuals
plot(fitted(Cond.a.AC.glm.gaus), resid(Cond.a.AC.glm.gaus))
abline(0,0)

qqnorm(resid(Cond.a.AC.glm.gaus))
qqline(resid(Cond.a.AC.glm.gaus))

plot(density(resid(Cond.a.AC.glm.gaus)))

##Model results
summary(Cond.a.AC.glm.gaus)
Anova(Cond.a.AC.glm.gaus)
#             LR Chisq Df Pr(>Chisq)    
# Site         290.870  2  < 2.2e-16 ***
# Season        34.981  2  2.535e-08 ***
# Site:Season    9.496  4    0.04983 *     


####Plot Average Conductivity by Site and Season####
Cond.a.plot<-ggplot(TS.sum, aes(x=Site.Sea, y=Cond.a, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), axis.text.x = element_blank(), 
        axis.text.y = element_text(size=axis.txt.sz, colour="black"), legend.position="none" )+
  ylab("Daily Average Conductivity (\u03BCS/cm)")+
  xlab(NULL)+
  ylim(44, 56.5)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 56.5, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 56.5, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 56.5, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 56.5, label = "Summer", size=sig.sz)

Cond.a.plot


####Daily Conductivity Standard Deviation by Site and Season####

####Dates complete at all Sites

##Check normality
hist(TS.sum.C$Cond.sd)
qqnorm(TS.sum.C$Cond.sd)
shapiro.test(TS.sum.C$Cond.sd)
#Not normal

##Compare models
Cond.sd.glm.gaus<-glm(Cond.sd~Site*Season, data=TS.sum.C, family=gaussian)
Cond.sd.glm.gam<-glm(Cond.sd~Site*Season, data=TS.sum.C, family=Gamma)
Cond.sd.glm.inv<-glm(Cond.sd~Site*Season, data=TS.sum.C, family=inverse.gaussian)

AIC(Cond.sd.glm.gaus, Cond.sd.glm.gam, Cond.sd.glm.inv)
##Gamma distribution fits best

##Check residuals
plot(fitted(Cond.sd.glm.gam), resid(Cond.sd.glm.gam))
abline(0,0)

qqnorm(resid(Cond.sd.glm.gam))
qqline(resid(Cond.sd.glm.gam))

plot(density(resid(Cond.sd.glm.gam)))

##Model results
summary(Cond.sd.glm.gam)
Anova(Cond.sd.glm.gam)
#             LR Chisq Df Pr(>Chisq)    
# Site         211.615  3  < 2.2e-16 ***
# Season         0.227  1   0.633770    
# Site:Season   14.392  3   0.002417 ** 


####Dates complete at AC Sites

##Check normality
hist(TS.sum.C.AC$Cond.sd)
qqnorm(TS.sum.C.AC$Cond.sd)
shapiro.test(TS.sum.C.AC$Cond.sd)
#Not normal

##Compare models
Cond.sd.AC.glm.gaus<-glm(Cond.sd~Site*Season, data=TS.sum.C.AC, family=gaussian)
Cond.sd.AC.glm.gam<-glm(Cond.sd~Site*Season, data=TS.sum.C.AC, family=Gamma)
Cond.sd.AC.glm.inv<-glm(Cond.sd~Site*Season, data=TS.sum.C.AC, family=inverse.gaussian)

AIC(Cond.sd.AC.glm.gaus, Cond.sd.AC.glm.gam, Cond.sd.AC.glm.inv)
##Gamma distribution fits best

##Check residuals
plot(fitted(Cond.sd.AC.glm.gam), resid(Cond.sd.AC.glm.gam))
abline(0,0)

qqnorm(resid(Cond.sd.AC.glm.gam))
qqline(resid(Cond.sd.AC.glm.gam))

plot(density(resid(Cond.sd.AC.glm.gam)))

##Model results
summary(Cond.sd.AC.glm.gam)
Anova(Cond.sd.AC.glm.gam)
#             LR Chisq Df Pr(>Chisq)    
# Site          399.82  2  < 2.2e-16 ***
# Season         29.97  2  3.099e-07 ***
# Site:Season    37.81  4  1.224e-07 ***   


####Plot Conductivity Standard Deviation by Site and Season####
Cond.sd.plot<-ggplot(TS.sum, aes(x=Site.Sea, y=Cond.sd, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), legend.position="none", 
        axis.text.x = element_text(size=axis.txt.sz, colour="black"), axis.text.y = element_text(size=axis.txt.sz, colour="black"))+
  ylab("Daily Conductivity St. Dev. (\u03BCS/cm)")+
  xlab(NULL)+
  ylim(0,4.2)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 4.2, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 4.2, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 4.2, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 4.2, label = "Summer", size=sig.sz)

Cond.sd.plot


####**Environment Part 3: Nutrients####

##Load and Clean Data
Nutr<-read.csv("Total_Nutrients.csv", header=TRUE)
Nutr$Site<-factor(Nutr$Site, levels=c("BD", "SS", "KL", "OL"), ordered=TRUE)
Nutr$Month.a<-factor(Nutr$Month.a, levels=c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"), ordered=TRUE)
Nutr$Month.n<-factor(Nutr$Month.n, levels=c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12"), ordered=TRUE)

##Add column for Season
Nutr$Season<-c(rep("Fall", nrow(Nutr)))
Nutr$Season[which(Nutr$Month.a=="Dec" | Nutr$Month.a=="Jan" | Nutr$Month.a=="Feb")]<-"Winter"
Nutr$Season[which(Nutr$Month.a=="Mar" | Nutr$Month.a=="Apr" | Nutr$Month.a=="May")]<-"Spring"
Nutr$Season[which(Nutr$Month.a=="Jun" | Nutr$Month.a=="Jul" | Nutr$Month.a=="Aug")]<-"Summer"
Nutr$Season<-factor(Nutr$Season, levels=c("Fall", "Winter", "Spring", "Summer"), ordered=TRUE)

##Add column for Site and Season
Nutr$Site.Sea<-paste(Nutr$Site, Nutr$Season, sep=".")
Nutr$Site.Sea<-factor(Nutr$Site.Sea, levels=c("BD.Fall", "SS.Fall", "KL.Fall", "OL.Fall",
                                              "BD.Winter", "SS.Winter", "KL.Winter", "OL.Winter",
                                              "BD.Spring", "SS.Spring", "KL.Spring", "OL.Spring",
                                              "BD.Summer", "SS.Summer", "KL.Summer", "OL.Summer"), ordered=TRUE)
##Add column for N:P Ratio
Nutr$N.P<-c(Nutr$TN_uM / Nutr$TP_uM)

##Check for outliers
boxplot(Nutr$TN_uM)
boxplot(Nutr$TP_uM)
boxplot(Nutr$N.P)

##Remove samples with TN > 20uM, TP > 0.6, and N.P > 100
Nutr.o<-Nutr[-c(which(Nutr$TN_uM>20 | Nutr$TP_uM>0.5 |  Nutr$N.P>100)),]

##Remove samples with TN = 0
Nutr.o<-Nutr.o[-c(which(Nutr.o$TN_uM==0)),]


####Total Nitrogen by Site and Season####

##Check normality
hist(Nutr.o$TN_uM)
qqnorm(Nutr.o$TN_uM)
shapiro.test(Nutr.o$TN_uM)
#Not normal

##Compare models
N.glm.gaus<-glm(TN_uM~Site*Season, data=Nutr.o, family=gaussian)
N.glm.gam<-glm(TN_uM~Site*Season, data=Nutr.o, family=Gamma)
N.glm.inv<-glm(TN_uM~Site*Season, data=Nutr.o, family=inverse.gaussian)

AIC(N.glm.gaus, N.glm.gam, N.glm.inv)
##Gamma distribution fits best

##Check residuals
plot(fitted(N.glm.gam), resid(N.glm.gam))
abline(0,0)

qqnorm(resid(N.glm.gam))
qqline(resid(N.glm.gam))

plot(density(resid(N.glm.gam)))

##Model results
summary(N.glm.gam)
Anova(N.glm.gam)
#             LR Chisq Df Pr(>Chisq)
# Site          4.2701  3     0.2337
# Season        5.2725  3     0.1529
# Site:Season   7.3021  9     0.6057


####Plot Total Nitrogen by Site and Season####
Nitrogen.plot<-ggplot(Nutr.o, aes(x=Site.Sea, y=TN_uM, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), legend.position="none", 
        axis.text.x = element_text(size=axis.txt.sz, colour="black"), axis.text.y = element_text(size=axis.txt.sz, colour="black"))+
  ylab("Total Nitrogen (\u03BCM)")+
  xlab(NULL)+
  ylim(0,12.5)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 12.5, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 12.5, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 12.5, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 12.5, label = "Summer", size=sig.sz)

Nitrogen.plot


####Total Phosphorus by Site and Season####

##Check normality
hist(Nutr.o$TP_uM)
qqnorm(Nutr.o$TP_uM)
shapiro.test(Nutr.o$TP_uM)
#Not normal

##Compare models
P.glm.gaus<-glm(TP_uM~Site*Season, data=Nutr.o, family=gaussian)
P.glm.gam<-glm(TP_uM~Site*Season, data=Nutr.o, family=Gamma)
P.glm.inv<-glm(TP_uM~Site*Season, data=Nutr.o, family=inverse.gaussian)

AIC(P.glm.gaus, P.glm.gam, P.glm.inv)
##Inverse Gaussian distribution fits best

##Check residuals
plot(fitted(P.glm.inv), resid(P.glm.inv))
abline(0,0)

qqnorm(resid(P.glm.inv))
qqline(resid(P.glm.inv))

plot(density(resid(P.glm.inv)))

##Model results
summary(P.glm.inv)
Anova(P.glm.inv)
#             LR Chisq Df Pr(>Chisq)
# Site          2.5441  3     0.4674
# Season        5.1972  3     0.1579
# Site:Season   5.8643  9     0.7534


####Plot Total Phosphorus by Site and Season####
Phosphorus.plot<-ggplot(Nutr.o, aes(x=Site.Sea, y=TP_uM, fill=Site)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), legend.position="none", 
        axis.text.x = element_text(size=axis.txt.sz, colour="black"), axis.text.y = element_text(size=axis.txt.sz, colour="black"))+
  ylab("Total Phosphorus (\u03BCM)")+
  xlab(NULL)+
  ylim(0.1,0.45)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 0.45, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 0.45, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 0.45, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 0.45, label = "Summer", size=sig.sz)

Phosphorus.plot


####Total N to P by Site and Season####

##Check normality
hist(Nutr.o$N.P)
qqnorm(Nutr.o$N.P)
shapiro.test(Nutr.o$N.P)
#Not normal

##Compare models
N.P.glm.gaus<-glm(N.P~Site*Season, data=Nutr.o, family=gaussian)
N.P.glm.gam<-glm(N.P~Site*Season, data=Nutr.o, family=Gamma)
N.P.glm.inv<-glm(N.P~Site*Season, data=Nutr.o, family=inverse.gaussian)

AIC(N.P.glm.gaus, N.P.glm.gam, N.P.glm.inv)
##Gamma distribution fits best

##Check residuals
plot(fitted(N.P.glm.gam), resid(N.P.glm.gam))
abline(0,0)

qqnorm(resid(N.P.glm.gam))
qqline(resid(N.P.glm.gam))

plot(density(resid(N.P.glm.gam)))

##Model results
summary(N.P.glm.gam)
Anova(N.P.glm.gam)
#             LR Chisq Df Pr(>Chisq)
# Site          5.4457  3     0.1419
# Season        4.3791  3     0.2233
# Site:Season   6.0526  9     0.7346


####Nitrogen Enrichment####

####N Enrichment N:P>16:1

##Instances of N Enrichment by Site 
N_Enriched.Site<-aggregate(Nutr.o[,9], list(Nutr.o$Site), function(x) length(which(x>16)))
names(N_Enriched.Site)<-c("Site", "N_Enriched")
N_Enriched.Site$Total<-aggregate(Nutr.o[,9], list(Nutr.o$Site), function(x) length(x))[,2]

##Test difference in proportions between Sites  
prop.test(N_Enriched.Site$N_Enriched, N_Enriched.Site$Total)
# X-squared = 5.8011, df = 3, p-value = 0.1217
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.6521739 0.5652174 0.3478261 0.3809524 


##Instances of N Enrichment by Season 
N_Enriched.Seas<-aggregate(Nutr.o[,9], list(Nutr.o$Season), function(x) length(which(x>16)))
names(N_Enriched.Seas)<-c("Season", "N_Enriched")
N_Enriched.Seas$Total<-aggregate(Nutr.o[,9], list(Nutr.o$Season), function(x) length(x))[,2]

##Test difference in proportions between Seasons  
prop.test(N_Enriched.Seas$N_Enriched, N_Enriched.Seas$Total)
# X-squared = 7.6498, df = 3, p-value = 0.05383
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.7083333 0.5000000 0.3333333 0.3809524 


####P Starvation N:P>22:1

##Instances of P Starvation by Site 
P_Starved.Site<-aggregate(Nutr.o[,9], list(Nutr.o$Site), function(x) length(which(x>22)))
names(P_Starved.Site)<-c("Site", "P_Starved")
P_Starved.Site$Total<-aggregate(Nutr.o[,9], list(Nutr.o$Site), function(x) length(x))[,2]

##Test difference in proportions between Sites  
prop.test(P_Starved.Site$P_Starved, P_Starved.Site$Total)
# X-squared = 7.8839, df = 3, p-value = 0.04847
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.4782609 0.4347826 0.1739130 0.1904762 

pairwise.prop.test(P_Starved.Site$P_Starved, P_Starved.Site$Total, "bonferroni")
# 1    2    3   
# 2 1.00 -    -   
# 3 0.35 0.65 -   
# 4 0.54 0.94 1.00


##Instances of P Starvation by Season 
P_Starved.Seas<-aggregate(Nutr.o[,9], list(Nutr.o$Season), function(x) length(which(x>22)))
names(P_Starved.Seas)<-c("Season", "P_Starved")
P_Starved.Seas$Total<-aggregate(Nutr.o[,9], list(Nutr.o$Site), function(x) length(x))[,2]

##Test difference in proportions between Seasons  
prop.test(P_Starved.Seas$P_Starved, P_Starved.Seas$Total)
# X-squared = 4.0161, df = 3, p-value = 0.2597
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1    prop 2    prop 3    prop 4 
# 0.4347826 0.3043478 0.1739130 0.3809524 


####Plot Total Nitrogen:Phosphorous by Site and Season####
Nutrients.plot<-ggplot(Nutr.o, aes(x=Site.Sea, y=N.P, fill=Site)) + 
  annotate("rect", xmin=0.4, xmax=16.6, ymin=0, ymax=16, alpha=0.3, fill="grey")+
  annotate("rect", xmin=0.4, xmax=16.6, ymin=16, ymax=22, alpha=0.3, fill="#C43C4EFF")+
  geom_hline(yintercept=16,linetype=1, color="grey", size=1)+
  geom_hline(yintercept=22,linetype=1, color="#C43C4EFF", size=1)+
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7,  size=2)+
  geom_boxplot(alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values =Site.colors.o)+
  theme_classic()+
  scale_x_discrete(labels=c(rep(c("BD", "SS", "KL", "OL"), 4)))+
  theme(axis.title.y = element_text(size = axis.title.sz), axis.text.x = element_blank(), 
        axis.text.y = element_text(size=axis.txt.sz, colour="black"), legend.position="none" )+
  ylab("Total N : Total P (\u03BCM)")+
  xlab(NULL)+
  ylim(0, 66)+
  geom_vline(xintercept=4.5,linetype=1)+
  geom_vline(xintercept=8.5,linetype=1)+
  geom_vline(xintercept=12.5,linetype=1)+
  annotate("text", x = 2.5, y = 66, label = "Fall", size=sig.sz)+
  annotate("text", x = 6.5, y = 66, label = "Winter", size=sig.sz)+
  annotate("text", x = 10.5, y = 66, label = "Spring", size=sig.sz)+
  annotate("text", x = 14.5, y = 66, label = "Summer", size=sig.sz)


Nutrients.plot

####***Figure 1 and S2: Plot Environmental Conditions####

##Figure 1 Sampling, Average Temp, Average Cond, N:P
Env.alignment<-align_plots(Temp.ribbon.plot, Temp.a.plot, align="v", axis="1")

Env.Main_fig.bottom<-plot_grid(Env.alignment[[2]], Cond.a.plot, Nutrients.plot, 
                               rel_widths=c(1, 1, 1), rel_heights=c(1, 1, 1), 
                               ncol=3, byrow=T, labels = c('B', 'C', 'D'), 
                               label_size=panel.lab.sz, align="h")

Env.Main_fig<-plot_grid(Env.alignment[[1]], Env.Main_fig.bottom, 
                        rel_widths=c(1, 1), rel_heights=c(1, 0.9), 
                        nrow=2, byrow=T, labels = c('A', ''), 
                        label_size=panel.lab.sz)
ggsave(filename="Fig1_Environmental_Conditions.png", plot=Env.Main_fig, dpi=300, width=14, height=10, units="in")
ggsave(filename="Fig1_Environmental_Conditions.pdf", plot=Env.Main_fig, dpi=300, width=14, height=10, units="in")


##Figure S2 Temp SD, Cond SD, Total N, Total P
Env.Supp_fig<-plot_grid(Temp.sd.plot, Cond.sd.plot, Nitrogen.plot, Phosphorus.plot, 
                        rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                        nrow=2, ncol=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                        label_size=panel.lab.sz, align="hv")

ggsave(filename="FigS2_Environmental_Conditions.png", plot=Env.Supp_fig, dpi=300, width=12, height=10, units="in")
ggsave(filename="FigS2_Environmental_Conditions.pdf", plot=Env.Supp_fig, dpi=300, width=12, height=10, units="in")



#-----*Coral Physiology-----------------------------------------------------------


####**Physiology Part 1: Physiological Metric Calculations####


####Coral Meta Data####

##Load Data
CoralMeta<-read.csv("CoralMeta.csv",  header=TRUE)
str(CoralMeta)

##Factors
CoralMeta$Species<-factor(CoralMeta$Species, levels=c("Acer", "Apal"))
CoralMeta$Genotype<-factor(CoralMeta$Genotype, levels=c("AC8", "AC10", "AC12", "AC15", "AP1", "AP8", "AP10"), ordered=TRUE)
CoralMeta$Site<-factor(CoralMeta$Site, levels=c("BD", "SS", "KL", "OL"), ordered=TRUE)
CoralMeta$Season<-factor(CoralMeta$Season, levels=c("Fall", "Winter", "Spring", "Summer"), ordered=TRUE)
str(CoralMeta)

##Add Variable of Sample Set
CoralMeta$SampSet<-paste(CoralMeta$Site, CoralMeta$Genotype, CoralMeta$Season, sep=".")


####Surface Area####

##Load Standards Data
WaxStands<-read.csv("WaxStandards.csv", header=TRUE, stringsAsFactors = TRUE)
str(WaxStands)

##Calculate Surface Area of Cylinders
WaxStands$SA_cm2<-2*pi*(WaxStands$Dm_cm/2)*WaxStands$Ht_cm+2*pi*(WaxStands$Dm_cm/2)^2

##Calculate Difference in Wax Weight
WaxStands$Wax.D_g<- WaxStands$Wax.F_g-WaxStands$Wax.I_g

##Linear Model of SA as a Function of Wax Weight
wax.SA.lm  <- lm(SA_cm2~Wax.D_g, data=WaxStands)

summary(wax.SA.lm)
coef(wax.SA.lm)

SA.mod <- function(wax.weight) {
  coefs <- coef(wax.SA.lm)
  #y = mx + b
  SA <- ((coefs[2] * wax.weight) + coefs[1])
  return(SA)}

##Load Sample Data
SampDat<-read.csv("CoralSampleData.csv", header=TRUE)
str(SampDat)

##Calculate Difference in Wax Weight (g)
SampDat$Wax.D_g<- SampDat$Wax.F_g-SampDat$Wax.I_g

##Calculate Surface Area (cm^2)
SampDat$SA_cm2<-SA.mod(SampDat$Wax.D_g)

##Add Surface Area to Coral Data
names(SampDat)
CoralData<-merge(CoralMeta, SampDat[,c(1:3,7)])


####Protein####

##Load Data
Protein<-read.csv("Protein.csv", header=TRUE, stringsAsFactors = TRUE)
str(Protein)

##Standard Curve

##Subset Standards Data
Prot.St<- subset(Protein, Input=="Standard")

##Load Standard MetaData and Merge 
BCAStands.Meta<-read.csv("BCAStandards_Meta.csv", header=TRUE, stringsAsFactors = TRUE)
str(BCAStands.Meta)
BCAStands<-merge(Prot.St, BCAStands.Meta, all.x=TRUE)

##Plot Data
plot(BCAStands$A562, BCAStands$Protein_ug.ml, pch=19)

##Fit Models of Protein Concentration as a function of Absorbance 

#Fit first degree polynomial equation:
BCA.lm  <- lm(Protein_ug.ml~A562, data=BCAStands)
#Second degree
BCA.ploy2 <- lm(Protein_ug.ml~poly(A562,2,raw=TRUE), data=BCAStands)
#Third degree
BCA.ploy3 <- lm(Protein_ug.ml~poly(A562,3,raw=TRUE), data=BCAStands)

##Fit Models to Plot
xx <- seq(0,3, length=50)

lines(xx, predict(BCA.lm, newdata=data.frame(A562=xx)), col="red")
lines(xx, predict(BCA.ploy2, newdata=data.frame(A562=xx)), col="green")
lines(xx, predict(BCA.ploy3, newdata=data.frame(A562=xx)), col="blue")

##Compare Models for Best Fit
anova(BCA.lm, BCA.ploy2)
#Poly2 is a significantly better fit than linear

anova(BCA.ploy2, BCA.ploy3)
#Poly3 is a significantly better fit than Ploy2

##Create a function with the Poly3 model Equation
summary(BCA.ploy3)
coef(BCA.ploy3)

TP.mod <- function(absorbance) {
  coefs <- coef(BCA.ploy3)
  #y = d + cx + bx^2 + ax^3
  protein <- coefs[1] + (coefs[2] * absorbance) + (coefs[3] * absorbance^2) + (coefs[4] * absorbance^3)
  return(protein)}

##Subset Sample Data (Unknown Sample)
Protein.Un<- subset(Protein, Input=="Un")

##Calculate Total Protein Concentration (ug/ml)
Protein.Un$TP_ug.ml<-TP.mod(Protein.Un$A562)

##Merge with Sample Data to Calculate Total Protein per Surface Area
Protein.Un<-merge(Protein.Un, SampDat, all=TRUE)

##Calculate Total Protein (ug) 
Protein.Un$TP_ug<-Protein.Un$TP_ug.ml*Protein.Un$Vol_ml

##Calculate Total Protein per Surface Area (ug/cm^2)
Protein.Un$TP_ug.cm2<-Protein.Un$TP_ug/Protein.Un$SA_cm2

##Check for Outliers within a Sample Set
Protein.Un<-merge(Protein.Un, CoralMeta, all.x=TRUE)

##Separate Coral and Symbiont Fractions
Protein.C<- subset(Protein.Un, Fraction=="C")
names(Protein.C)
names(Protein.C)[19]<-paste(names(Protein.C)[19], "C",  sep="_")

Protein.S<- subset(Protein.Un, Fraction=="Z")
names(Protein.S)
names(Protein.S)[19]<-paste(names(Protein.S)[19], "S",  sep="_")

##Boxplots to Check for Outliers Coral Host
boxplot(Protein.C$TP_ug.cm2_C~Protein.C$Species)
abline(h=800, col="red", lty=2)
Protein.C$ID[which(Protein.C$TP_ug.cm2_C>800)]

##Remove Readings above 800 ug/cm2
Protein.Co<-Protein.C[-c(which(Protein.C$TP_ug.cm2_C>800)),]

##Boxplots to Check for Outliers Symbiont
boxplot(Protein.S$TP_ug.cm2_S~Protein.S$Species)
abline(h=900, col="red", lty=2)
Protein.S$ID[which(Protein.S$TP_ug.cm2_S>900)]

##Remove Readings above 900 ug/cm
Protein.So<-Protein.S[-c(which(Protein.S$TP_ug.cm2_S>900)),]

##Average Across Readings
names(Protein.Co)
Protein.Ca<-aggregate(Protein.Co$TP_ug.cm2_C, list(Protein.Co$RandN, Protein.Co$ID), mean)
names(Protein.Ca)<-c("RandN", "ID", "TP_ug.cm2_C")

names(Protein.So)
Protein.Sa<-aggregate(Protein.So$TP_ug.cm2_S, list(Protein.So$RandN, Protein.So$ID), mean)
names(Protein.Sa)<-c("RandN", "ID", "TP_ug.cm2_S")

##Add Total Protein of Coral Host (C) and Symbiont (S) to Coral Data
CoralData<-merge(CoralData, Protein.Ca, all.x=TRUE)
CoralData<-merge(CoralData, Protein.Sa, all.x=TRUE)


####Chlorophyll####
# Equations for Dinos from Jeffrey and Humphrey 1975 in 100% acetone
# chla = 11.43*A663 - 0.64*A630
# chlc2 = 27.09*A630 - 3.63*A663

##Load Data
Chl<-read.csv("Chlorophyll.csv",  header=TRUE)
str(Chl)

##Subtract Background A750 from A630 and A663
Chl$A630.c<-Chl$A630-Chl$A750
Chl$A663.c<-Chl$A663-Chl$A750

##Divide by Pathlength (0.5cm pathlength for 175ul sample in UVStar Plate) 
Chl$A630.c<-c(Chl$A630.c/0.5)
Chl$A663.c<-c(Chl$A663.c/0.5)

##Calculate Chl-a and Chl-c2 in µg/ml
Chl$Chl.a_ug.ml<-11.43*Chl$A663.c- 0.64*Chl$A630.c
Chl$Chl.c2_ug.ml <- 27.09*Chl$A630.c - 3.63*Chl$A663.c

##Merge with Sample Data to Calculate Chlorophyll per Surface Area
Chl<-merge(Chl, SampDat,  all=TRUE)

##Calculate Total Chlorophyll-a and c2 (ug) 
Chl$Chl.a_ug<-Chl$Chl.a_ug.ml*Chl$Vol_ml
Chl$Chl.c2_ug<-Chl$Chl.c2_ug.ml*Chl$Vol_ml

##Calculate Chlorophyll-a and c2 per Surface Area (ug/cm^2)
Chl$Chl.a_ug.cm2<-Chl$Chl.a_ug/Chl$SA_cm2
Chl$Chl.c2_ug.cm2<-Chl$Chl.c2_ug/Chl$SA_cm2

##Check for Outliers within a Sample Set
Chl<-merge(Chl, CoralMeta, all.x=TRUE)

##Boxplots to Check for Outliers

#Chl-a
boxplot(Chl$Chl.a_ug.cm2~Chl$Species)
abline(h=5, col="red", lty=2)
Chl$ID[which(Chl$Chl.a_ug.cm2>5)]

##Remove Readings with Chl-a above 5 ug/cm2
Chl.o<-Chl[-c(which(Chl$Chl.a_ug.cm2>5)),]

#Chl-c2
boxplot(Chl$Chl.c2_ug.cm2~Chl$Species)

abline(h=2, col="red", lty=2)
Chl$ID[which(Chl$Chl.c2_ug.cm2>2 & Chl$Species=="Acer")]

abline(h=7, col="red", lty=2)
Chl$ID[which(Chl$Chl.c2_ug.cm2>7 & Chl$Species=="Apal")]

##Remove Readings with Chl-c2 above 2 ug/cm2 (Acer) and 7ug/cm2 (Apal)
Chl.o<-Chl.o[-c(which(Chl.o$Chl.c2_ug.cm2>2 & Chl.o$Species=="Acer")),]
Chl.o<-Chl.o[-c(which(Chl.o$Chl.c2_ug.cm2>7 & Chl.o$Species=="Apal")),]

#Also remove readings with negative c2
Chl.o<-Chl.o[-c(which(Chl.o$Chl.c2_ug.cm2<0)),]

##Average Across Readings
names(Chl.o)
Chl.a<-aggregate(Chl.o[,c(22:23)], list(Chl.o$RandN, Chl.o$ID), mean, na.action=na.omit)
names(Chl.a)<-c("RandN", "ID", "Chl.a_ug.cm2", "Chl.c2_ug.cm2")

##Add Chlorophyll-a and c2 to Coral Data
CoralData<-merge(CoralData, Chl.a)


####Biomass####

##Load Data
Biomass<-read.csv("Biomass.csv",  header=TRUE)
str(Biomass)
Biomass$Fraction<-factor(Biomass$Fraction, levels=c("C", "Z"))

##Calculate Ash Free Dry Weight (g/ml)
Biomass$DeltaBurn_g<-Biomass$Dried_g-Biomass$Burned_g
Biomass$AFDW_g.ml<-Biomass$DeltaBurn_g/Biomass$InVol_ml

##Merge with Sample Data to Calculate Biomass per Surface Area
Biomass<-merge(Biomass, SampDat, all.x=TRUE)

##Calculate Total AFDW (g) 
Biomass$AFDW_g<-Biomass$AFDW_g.ml*Biomass$Vol_ml

##Calculate AFDW per Surface Area (mg/cm^2)
Biomass$AFDW_g.cm2<-Biomass$AFDW_g/Biomass$SA_cm2
Biomass$AFDW_mg.cm2 <- Biomass$AFDW_g.cm2 * 1000

##Separate Coral and Symbiont Fractions
Biomass.C<-subset(Biomass, Fraction=="C")
names(Biomass.C)
names(Biomass.C)[17]<-paste(names(Biomass.C)[17], "C",  sep="_")

Biomass.S<-subset(Biomass, Fraction=="Z")
names(Biomass.S)
names(Biomass.S)[17]<-paste(names(Biomass.S)[17], "S",  sep="_")

##Add AFDW of Coral and Symbiont to Coral Data
CoralData<-merge(CoralData, Biomass.C[,c(1, 17)], all.x=TRUE)
CoralData<-merge(CoralData, Biomass.S[,c(1, 17)], all.x=TRUE)



####**Physiology Part 2: Univariate Analysis####


####Physiology Univariate Analysis- Acropora cervicornis ####
CoralData.AC<-subset(CoralData, Species=="Acer")
rownames(CoralData.AC)<-paste("AC", CoralData.AC$RandN, sep="")


####Physiology- Protein- A. cervicornis ####

####Protein Host Fraction

##Check normality
hist(CoralData.AC$TP_ug.cm2_C)
qqnorm(CoralData.AC$TP_ug.cm2_C)
shapiro.test(CoralData.AC$TP_ug.cm2_C)
#Not normal

hist(log(CoralData.AC$TP_ug.cm2_C+1))
qqnorm(log(CoralData.AC$TP_ug.cm2_C+1))
shapiro.test(log(CoralData.AC$TP_ug.cm2_C+1))
#Normal

##Model with log+1 transformation
Prot.C.AC.log.int<-lm(log(TP_ug.cm2_C+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Prot.C.AC.log.int), resid(Prot.C.AC.log.int))
abline(0,0)

qqnorm(resid(Prot.C.AC.log.int))
qqline(resid(Prot.C.AC.log.int))

plot(density(resid(Prot.C.AC.log.int)))

##Model results
summary(Prot.C.AC.log.int)
# Multiple R-squared:  0.4053,	Adjusted R-squared:  0.3392 
# F-statistic: 6.133 on 23 and 207 DF,  p-value: 8.27e-14

anova(Prot.C.AC.log.int)
#                   Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype          3 0.3848 0.12827  3.8862  0.009888 ** 
# Site              2 1.3347 0.66733 20.2177 9.542e-09 ***
# Season            3 2.0155 0.67183 20.3539 1.342e-11 ***
# Genotype:Site     6 0.2656 0.04426  1.3410  0.240358    
# Genotype:Season   9 0.6554 0.07282  2.2062  0.022938 *  
# Residuals       207 6.8325 0.03301  

eta_squared(Prot.C.AC.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.03 | [0.00, 1.00]
# Site            | 0.12 | [0.05, 1.00]
# Season          | 0.18 | [0.10, 1.00]
# Genotype:Site   | 0.02 | [0.00, 1.00]
# Genotype:Season | 0.06 | [0.00, 1.00]

emmeans(Prot.C.AC.log.int, pairwise~Site)
# contrast estimate     SE  df t.ratio p.value
# BD - SS    0.0234 0.0300 207   0.780  0.7155
# BD - KL    0.1764 0.0300 207   5.880  <.0001 *
# SS - KL    0.1530 0.0287 207   5.325  <.0001 *

##Plot Reaction Norms
Prot.C.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="TP_ug.cm2_C", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Prot.C.AC_Geno.Site.plot<-ggplot(Prot.C.AC_Geno.Site, aes(x=Site, y=TP_ug.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_C-ci, ymax=TP_ug.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Protein (\u03BCg cm'^-2*") Host")))+
  ylim(225, 640)+
  annotate("text", x = 1, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 580, label = "B", size=levels.sz)+
  annotate("text", x=0.6, y=285, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.12"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=260, label=expression(italic(paste("Genotype x Site p = 0.240, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)

Prot.C.AC_Geno.Site.plot


Prot.C.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="TP_ug.cm2_C", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Prot.C.AC_Geno.Seas.plot<-ggplot(Prot.C.AC_Geno.Seas, aes(x=Season, y=TP_ug.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_C-ci, ymax=TP_ug.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Protein (\u03BCg cm'^-2*") Host")))+
  ylim(225, 640)+
  annotate("text", x=0.65, y=285, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.18"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.65, y=260, label=expression(bolditalic(paste("Genotype x Season p = 0.023, ", eta^2, " = 0.06"))), size=sig.sz, hjust = 0)

Prot.C.AC_Geno.Seas.plot


####Protein Symbiont Fraction

##Check normality
hist(CoralData.AC$TP_ug.cm2_S)
qqnorm(CoralData.AC$TP_ug.cm2_S)
shapiro.test(CoralData.AC$TP_ug.cm2_S)
#Not normal

hist(log(CoralData.AC$TP_ug.cm2_S+1))
qqnorm(log(CoralData.AC$TP_ug.cm2_S+1))
shapiro.test(log(CoralData.AC$TP_ug.cm2_S+1))
#Normal

##Model with log+1 transformation
Prot.S.AC.log.int<-lm(log(TP_ug.cm2_S+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Prot.S.AC.log.int), resid(Prot.S.AC.log.int))
abline(0,0)

qqnorm(resid(Prot.S.AC.log.int))
qqline(resid(Prot.S.AC.log.int))

plot(density(resid(Prot.S.AC.log.int)))

##Model results
summary(Prot.S.AC.log.int)
# Multiple R-squared:  0.5955,	Adjusted R-squared:  0.5506 
# F-statistic: 13.25 on 23 and 207 DF,  p-value: < 2.2e-16

anova(Prot.S.AC.log.int)
#                   Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype          3 0.6301 0.21003  4.8563  0.002756 ** 
# Site              2 2.6355 1.31773 30.4689  2.52e-12 ***
# Season            3 9.3549 3.11829 72.1020 < 2.2e-16 ***
# Genotype:Site     6 0.3157 0.05262  1.2167  0.298908    
# Genotype:Season   9 0.2449 0.02721  0.6292  0.771371    
# Residuals       207 8.9524 0.04325 

eta_squared(Prot.S.AC.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.03 | [0.00, 1.00]
# Site            | 0.12 | [0.05, 1.00]
# Season          | 0.42 | [0.34, 1.00]
# Genotype:Site   | 0.01 | [0.00, 1.00]
# Genotype:Season | 0.01 | [0.00, 1.00]

emmeans(Prot.S.AC.log.int, pairwise~Genotype)
# contrast     estimate     SE  df t.ratio p.value
# AC8 - AC10  -0.062269 0.0381 207  -1.632  0.3627
# AC8 - AC12  -0.138674 0.0392 207  -3.539  0.0028 *
# AC8 - AC15  -0.062771 0.0392 207  -1.602  0.3797
# AC10 - AC12 -0.076404 0.0390 207  -1.959  0.2072
# AC10 - AC15 -0.000501 0.0390 207  -0.013  1.0000
# AC12 - AC15  0.075903 0.0400 207   1.897  0.2327

emmeans(Prot.S.AC.log.int, pairwise~Site)
# contrast estimate     SE  df t.ratio p.value
# BD - SS   0.00784 0.0343 207   0.228  0.9717
# BD - KL   0.21394 0.0343 207   6.231  <.0001 *
# SS - KL   0.20610 0.0329 207   6.268  <.0001 *

emmeans(Prot.S.AC.log.int, pairwise~Season)
# contrast        estimate     SE  df t.ratio p.value
# Fall - Winter    -0.0254 0.0402 207  -0.631  0.9218
# Fall - Spring    -0.4841 0.0381 207 -12.690  <.0001 *
# Fall - Summer    -0.3107 0.0381 207  -8.146  <.0001 *
# Winter - Spring  -0.4587 0.0400 207 -11.460  <.0001 *
# Winter - Summer  -0.2854 0.0400 207  -7.130  <.0001 *
# Spring - Summer   0.1733 0.0380 207   4.565  0.0001 *

##Plot Reaction Norms
Prot.S.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="TP_ug.cm2_S", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Prot.S.AC_Geno.Site.plot<-ggplot(Prot.S.AC_Geno.Site, aes(x=Site, y=TP_ug.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_S-ci, ymax=TP_ug.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Protein (\u03BCg cm'^-2*") Symbiont")))+
  ylim(225, 640)+
  annotate("text", x = 1, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 580, label = "B", size=levels.sz)+
  annotate("text", x=0.6, y=285, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.12"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=260, label=expression(italic(paste("Genotype x Site p = 0.299, ", eta^2, " = 0.01"))), size=sig.sz, hjust = 0)

Prot.S.AC_Geno.Site.plot


Prot.S.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="TP_ug.cm2_S", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Prot.S.AC_Geno.Seas.plot<-ggplot(Prot.S.AC_Geno.Seas, aes(x=Season, y=TP_ug.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_S-ci, ymax=TP_ug.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Protein (\u03BCg cm'^-2*") Symbiont")))+
  ylim(225, 640)+
  annotate("text", x = 1, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 580, label = "A", size=levels.sz)+
  annotate("text", x = 2.9, y = 580, label = "B", size=levels.sz)+
  annotate("text", x = 4, y = 580, label = "C", size=levels.sz)+
  annotate("text", x=4.35, y=285, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.42"))), size=sig.sz, hjust = 1)+
  annotate("text", x=4.35, y=260, label=expression(italic(paste("Genotype x Season p = 0.771, ", eta^2, " = 0.01"))), size=sig.sz, hjust = 1)

Prot.S.AC_Geno.Seas.plot


####Physiology- Chlorophyll- A. cervicornis ####

####Chlorophyll-a

##Check normality
hist(CoralData.AC$Chl.a_ug.cm2)
qqnorm(CoralData.AC$Chl.a_ug.cm2)
shapiro.test(CoralData.AC$Chl.a_ug.cm2)
#Not normal

hist(log(CoralData.AC$Chl.a_ug.cm2+1))
qqnorm(log(CoralData.AC$Chl.a_ug.cm2+1))
shapiro.test(log(CoralData.AC$Chl.a_ug.cm2+1))
#Normal

##Model with log+1 transformation
Chl.a.AC.log.int<-lm(log(Chl.a_ug.cm2+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Chl.a.AC.log.int), resid(Chl.a.AC.log.int))
abline(0,0)

qqnorm(resid(Chl.a.AC.log.int))
qqline(resid(Chl.a.AC.log.int))

plot(density(resid(Chl.a.AC.log.int)))

##Model results
summary(Chl.a.AC.log.int)
# Multiple R-squared:  0.7044,	Adjusted R-squared:  0.6715 
# F-statistic: 21.44 on 23 and 207 DF,  p-value: < 2.2e-16

anova(Chl.a.AC.log.int)
#                   Df Sum Sq Mean Sq  F value    Pr(>F)    
# Genotype          3 1.4549 0.48497  24.2089 1.805e-13 ***
# Site              2 4.9921 2.49606 124.5984 < 2.2e-16 ***
# Season            3 2.6256 0.87519  43.6878 < 2.2e-16 ***
# Genotype:Site     6 0.3056 0.05093   2.5423  0.021384 *  
# Genotype:Season   9 0.5019 0.05577   2.7838  0.004231 ** 
# Residuals       207 4.1468 0.02003 

eta_squared(Chl.a.AC.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.10 | [0.04, 1.00]
# Site            | 0.36 | [0.27, 1.00]
# Season          | 0.19 | [0.11, 1.00]
# Genotype:Site   | 0.02 | [0.00, 1.00]
# Genotype:Season | 0.04 | [0.00, 1.00]

##Plot Reaction Norms
Chl.a.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="Chl.a_ug.cm2", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Chl.a.AC_Geno.Site.plot<-ggplot(Chl.a.AC_Geno.Site, aes(x=Site, y=Chl.a_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.a_ug.cm2-ci, ymax=Chl.a_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")))+
  ylim(0.55, 2.8)+
  annotate("text", x=0.6, y=0.9, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.36"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=0.75, label=expression(bolditalic(paste("Genotype x Site p = 0.021, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)

Chl.a.AC_Geno.Site.plot


Chl.a.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="Chl.a_ug.cm2", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Chl.a.AC_Geno.Seas.plot<-ggplot(Chl.a.AC_Geno.Seas, aes(x=Season, y=Chl.a_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.a_ug.cm2-ci, ymax=Chl.a_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")))+
  ylim(0.55, 2.8)+
  annotate("text", x=0.65, y=0.9, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.19"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.65, y=0.75, label=expression(bolditalic(paste("Genotype x Season p = 0.004, ", eta^2, " = 0.04"))), size=sig.sz, hjust = 0)

Chl.a.AC_Geno.Seas.plot


#####Chlorophyll-c2

##Check normality
hist(CoralData.AC$Chl.c2_ug.cm2)
qqnorm(CoralData.AC$Chl.c2_ug.cm2)
shapiro.test(CoralData.AC$Chl.c2_ug.cm2)
#Not normal

hist(log(CoralData.AC$Chl.c2_ug.cm2+1))
qqnorm(log(CoralData.AC$Chl.c2_ug.cm2+1))
shapiro.test(log(CoralData.AC$Chl.c2_ug.cm2+1))
#Still not normal, but less skewed

##Model with log+1 transformation
Chl.c2.AC.log.int<-lm(log(Chl.c2_ug.cm2+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Chl.c2.AC.log.int), resid(Chl.c2.AC.log.int))
abline(0,0)

qqnorm(resid(Chl.c2.AC.log.int))
qqline(resid(Chl.c2.AC.log.int))

plot(density(resid(Chl.c2.AC.log.int)))

##Model results
summary(Chl.c2.AC.log.int)
# Multiple R-squared:  0.6883,	Adjusted R-squared:  0.6537 
# F-statistic: 19.87 on 23 and 207 DF,  p-value: < 2.2e-16

anova(Chl.c2.AC.log.int)
#                   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype          3 0.62967 0.20989 22.2948 1.501e-12 ***
# Site              2 1.39534 0.69767 74.1074 < 2.2e-16 ***
# Season            3 1.92344 0.64115 68.1033 < 2.2e-16 ***
# Genotype:Site     6 0.12745 0.02124  2.2563  0.039390 *  
# Genotype:Season   9 0.22735 0.02526  2.6832  0.005716 ** 
# Residuals       207 1.94877 0.00941  

eta_squared(Chl.c2.AC.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.10 | [0.04, 1.00]
# Site            | 0.22 | [0.14, 1.00]
# Season          | 0.31 | [0.22, 1.00]
# Genotype:Site   | 0.02 | [0.00, 1.00]
# Genotype:Season | 0.04 | [0.00, 1.00]

##Plot Reaction Norms
Chl.c2.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="Chl.c2_ug.cm2", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Chl.c2.AC_Geno.Site.plot<-ggplot(Chl.c2.AC_Geno.Site, aes(x=Site, y=Chl.c2_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.c2_ug.cm2-ci, ymax=Chl.c2_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Chlorophyll-c2 (\u03BCg cm'^-2*")")))+
  ylim(0, 0.85)+
  annotate("text", x=0.6, y=0.14, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.22"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=0.08, label=expression(bolditalic(paste("Genotype x Site p = 0.039, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)

Chl.c2.AC_Geno.Site.plot


Chl.c2.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="Chl.c2_ug.cm2", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Chl.c2.AC_Geno.Seas.plot<-ggplot(Chl.c2.AC_Geno.Seas, aes(x=Season, y=Chl.c2_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.c2_ug.cm2-ci, ymax=Chl.c2_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Chlorophyll-c2 (\u03BCg cm'^-2*")")))+
  ylim(0, 0.85)+
  annotate("text", x=0.65, y=0.14, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.31"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.65, y=0.08, label=expression(bolditalic(paste("Genotype x Season p = 0.006, ", eta^2, " = 0.04"))), size=sig.sz, hjust = 0)

Chl.c2.AC_Geno.Seas.plot


####Physiology- Biomass- A. cervicornis ####

####Biomass Host Fraction

##Check normality
hist(CoralData.AC$AFDW_mg.cm2_C)
qqnorm(CoralData.AC$AFDW_mg.cm2_C)
shapiro.test(CoralData.AC$AFDW_mg.cm2_C)
#Normal

##Model with lm
Bio.C.AC.lm.int<-lm(AFDW_mg.cm2_C~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Bio.C.AC.lm.int), resid(Bio.C.AC.lm.int))
abline(0,0)

qqnorm(resid(Bio.C.AC.lm.int))
qqline(resid(Bio.C.AC.lm.int))

plot(density(resid(Bio.C.AC.lm.int)))

##Model results
summary(Bio.C.AC.lm.int)
# Multiple R-squared:  0.2933,	Adjusted R-squared:  0.2148 
# F-statistic: 3.736 on 23 and 207 DF,  p-value: 1.753e-07

anova(Bio.C.AC.lm.int)
#                   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype          3  1.6420 0.54733  6.3915 0.0003669 ***
# Site              2  1.9989 0.99944 11.6709 1.575e-05 ***
# Season            3  1.3579 0.45265  5.2858 0.0015657 ** 
# Genotype:Site     6  1.1857 0.19762  2.3077 0.0353323 *  
# Genotype:Season   9  1.1732 0.13035  1.5222 0.1417527    
# Residuals       207 17.7265 0.08564     

eta_squared(Bio.C.AC.lm.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.07 | [0.01, 1.00]
# Site            | 0.08 | [0.03, 1.00]
# Season          | 0.05 | [0.01, 1.00]
# Genotype:Site   | 0.05 | [0.00, 1.00]
# Genotype:Season | 0.05 | [0.00, 1.00]

emmeans(Bio.C.AC.lm.int, pairwise~Season)
# contrast        estimate     SE  df t.ratio p.value
# Fall - Winter     0.0233 0.0566 207   0.412  0.9764
# Fall - Spring    -0.1851 0.0537 207  -3.449  0.0038 *
# Fall - Summer    -0.0646 0.0537 207  -1.204  0.6248
# Winter - Spring  -0.2084 0.0563 207  -3.700  0.0016 *
# Winter - Summer  -0.0879 0.0563 207  -1.561  0.4030
# Spring - Summer   0.1205 0.0534 207   2.255  0.1123

##Plot Reaction Norms
Bio.C.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="AFDW_mg.cm2_C", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Bio.C.AC_Geno.Site.plot<-ggplot(Bio.C.AC_Geno.Site, aes(x=Site, y=AFDW_mg.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_C-ci, ymax=AFDW_mg.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Biomass (AFDW mg cm'^-2*") Host")))+
  ylim(0.75, 2.2)+
  annotate("text", x=0.6, y=1, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.08"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=0.9, label=expression(bolditalic(paste("Genotype x Site p = 0.035, ", eta^2, " = 0.05"))), size=sig.sz, hjust = 0)

Bio.C.AC_Geno.Site.plot


Bio.C.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="AFDW_mg.cm2_C", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Bio.C.AC_Geno.Seas.plot<-ggplot(Bio.C.AC_Geno.Seas, aes(x=Season, y=AFDW_mg.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_C-ci, ymax=AFDW_mg.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Biomass (AFDW mg cm'^-2*") Host")))+
  ylim(0.75,2.2)+
  annotate("text", x = 1, y = 2, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 2, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 2, label = "B", size=levels.sz)+
  annotate("text", x = 4, y = 2, label = "AB", size=levels.sz)+
  annotate("text", x=0.6, y=1, label=expression(bolditalic(paste("Season p = 0.002, ", eta^2, " = 0.05"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=0.9, label=expression(italic(paste("Genotype x Season p = 0.142, ", eta^2, " = 0.05"))), size=sig.sz, hjust = 0)

Bio.C.AC_Geno.Seas.plot


####Biomass Symbiont Fraction

##Check normality
hist(CoralData.AC$AFDW_mg.cm2_S)
qqnorm(CoralData.AC$AFDW_mg.cm2_S)
shapiro.test(CoralData.AC$AFDW_mg.cm2_S)
#Normal

##Model with lm 
Bio.S.AC.lm.int<-lm(AFDW_mg.cm2_S~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AC)

##Check residuals
plot(fitted(Bio.S.AC.lm.int), resid(Bio.S.AC.lm.int))
abline(0,0)

qqnorm(resid(Bio.S.AC.lm.int))
qqline(resid(Bio.S.AC.lm.int))

plot(density(resid(Bio.S.AC.lm.int)))

##Model results
summary(Bio.S.AC.lm.int)
# Multiple R-squared:  0.2321,	Adjusted R-squared:  0.1468 
# F-statistic:  2.72 on 23 and 207 DF,  p-value: 9.147e-05

anova(Bio.S.AC.lm.int)
#                   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype          3 0.1149 0.03829  0.9917   0.39768    
# Site              2 1.3787 0.68933 17.8515 7.041e-08 ***
# Season            3 0.3606 0.12019  3.1125   0.02730 *  
# Genotype:Site     6 0.4183 0.06972  1.8055   0.09946 .  
# Genotype:Season   9 0.1433 0.01593  0.4124   0.92762    
# Residuals       207 7.9932 0.03861      

eta_squared(Bio.S.AC.lm.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.01 | [0.00, 1.00]
# Site            | 0.13 | [0.07, 1.00]
# Season          | 0.03 | [0.00, 1.00]
# Genotype:Site   | 0.04 | [0.00, 1.00]
# Genotype:Season | 0.01 | [0.00, 1.00]

emmeans(Bio.S.AC.lm.int, pairwise~Site)
# contrast estimate     SE  df t.ratio p.value
# BD - SS    0.0506 0.0324 207   1.561  0.2651
# BD - KL    0.1828 0.0324 207   5.636  <.0001 *
# SS - KL    0.1322 0.0311 207   4.255  0.0001 *

emmeans(Bio.S.AC.lm.int, pairwise~Season)
# contrast        estimate     SE  df t.ratio p.value
# Fall - Winter   -0.00518 0.0380 207  -0.137  0.9991
# Fall - Spring   -0.04607 0.0360 207  -1.278  0.5778
# Fall - Summer    0.06212 0.0360 207   1.724  0.3141
# Winter - Spring -0.04089 0.0378 207  -1.081  0.7014
# Winter - Summer  0.06731 0.0378 207   1.780  0.2860
# Spring - Summer  0.10819 0.0359 207   3.016  0.0152 *

##Plot Reaction Norms
Bio.S.AC_Geno.Site<-summarySE(CoralData.AC, measurevar="AFDW_mg.cm2_S", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Bio.S.AC_Geno.Site.plot<-ggplot(Bio.S.AC_Geno.Site, aes(x=Site, y=AFDW_mg.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_S-ci, ymax=AFDW_mg.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Biomass (AFDW mg cm'^-2*") Symbiont")))+
  ylim(0.3, 1.1)+
  annotate("text", x = 1.1, y = 1, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 1, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 1, label = "B", size=levels.sz)+
  annotate("text", x=0.6, y=0.42, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.13"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.6, y=0.37, label=expression(italic(paste("Genotype x Site p = 0.099, ", eta^2, " = 0.04"))), size=sig.sz, hjust = 0)

Bio.S.AC_Geno.Site.plot


Bio.S.AC_Geno.Seas<-summarySE(CoralData.AC, measurevar="AFDW_mg.cm2_S", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Bio.S.AC_Geno.Seas.plot<-ggplot(Bio.S.AC_Geno.Seas, aes(x=Season, y=AFDW_mg.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AC.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_S-ci, ymax=AFDW_mg.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Biomass (AFDW mg cm'^-2*") Symbiont")))+
  ylim(0.3, 1.1)+
  annotate("text", x = 1, y = 1, label = "AB", size=levels.sz)+
  annotate("text", x = 2, y = 1, label = "AB", size=levels.sz)+
  annotate("text", x = 3, y = 1, label = "A", size=levels.sz)+
  annotate("text", x = 4, y = 1, label = "B", size=levels.sz)+
  annotate("text", x=0.65, y=0.42, label=expression(bolditalic(paste("Season p = 0.027, ", eta^2, " = 0.03"))), size=sig.sz, hjust = 0)+
  annotate("text", x=0.65, y=0.37, label=expression(italic(paste("Genotype x Season p = 0.928, ", eta^2, " = 0.01"))), size=sig.sz, hjust = 0)

Bio.S.AC_Geno.Seas.plot


####***Figure S3, S4, S5: Plot Physiology Univariate Panels- A. cervicornis ####

##Biomass
Phys.Bio.AC_fig<-plot_grid(Bio.C.AC_Geno.Seas.plot, Bio.C.AC_Geno.Site.plot,
                           Bio.S.AC_Geno.Seas.plot, Bio.S.AC_Geno.Site.plot, 
                           rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                           ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                           label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS3_Acer_Physiology_ReactionNorms_Bio.png", plot=Phys.Bio.AC_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS3_Acer_Physiology_ReactionNorms_Bio.pdf", plot=Phys.Bio.AC_fig, dpi=300, width=12, height=12, units="in")


##Protein
Phys.Prot.AC_fig<-plot_grid(Prot.C.AC_Geno.Seas.plot, Prot.C.AC_Geno.Site.plot,
                            Prot.S.AC_Geno.Seas.plot, Prot.S.AC_Geno.Site.plot, 
                            rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                            ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                            label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS4_Acer_Physiology_ReactionNorms_Prot.png", plot=Phys.Prot.AC_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS4_Acer_Physiology_ReactionNorms_Prot.pdf", plot=Phys.Prot.AC_fig, dpi=300, width=12, height=12, units="in")


##Chlorophyll
Phys.Chl.AC_fig<-plot_grid(Chl.a.AC_Geno.Seas.plot, Chl.a.AC_Geno.Site.plot,
                           Chl.c2.AC_Geno.Seas.plot, Chl.c2.AC_Geno.Site.plot, 
                           rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                           ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                           label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS5_Acer_Physiology_ReactionNorms_Chl.png", plot=Phys.Chl.AC_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS5_Acer_Physiology_ReactionNorms_Chl.pdf", plot=Phys.Chl.AC_fig, dpi=300, width=12, height=12, units="in")


####Physiology Univariate Analysis- Acropora palmata ####
CoralData.AP<-subset(CoralData, Species=="Apal")
rownames(CoralData.AP)<-paste("A", CoralData.AP$RandN, sep="")


####Physiology- Protein- A. palmata ####

####Protein Host Fraction

##Check normality
hist(CoralData.AP$TP_ug.cm2_C)
qqnorm(CoralData.AP$TP_ug.cm2_C)
shapiro.test(CoralData.AP$TP_ug.cm2_C)
#Normal

##Model with lm
Prot.C.AP.lm.int<-lm(TP_ug.cm2_C~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Prot.C.AP.lm.int), resid(Prot.C.AP.lm.int))
abline(0,0)

qqnorm(resid(Prot.C.AP.lm.int))
qqline(resid(Prot.C.AP.lm.int))

plot(density(resid(Prot.C.AP.lm.int)))

##Model results
summary(Prot.C.AP.lm.int)
# Multiple R-squared:  0.5457,	Adjusted R-squared:  0.459 
# F-statistic:  6.29 on 17 and 89 DF,  p-value: 2.236e-09

anova(Prot.C.AP.lm.int)
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  24015   12007  2.8436  0.063518 .  
# Site             2  61642   30821  7.2990  0.001161 ** 
# Season           3 324351  108117 25.6043 4.937e-12 ***
# Genotype:Site    4  15013    3753  0.8888  0.474083    
# Genotype:Season  6  26471    4412  1.0448  0.402003    
# Residuals       89 375812    4223 

eta_squared(Prot.C.AP.lm.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.03 | [0.00, 1.00]
# Site            | 0.07 | [0.00, 1.00]
# Season          | 0.39 | [0.25, 1.00]
# Genotype:Site   | 0.02 | [0.00, 1.00]
# Genotype:Season | 0.03 | [0.00, 1.00]

emmeans(Prot.C.AP.lm.int, pairwise~Genotype)
# contrast   estimate   SE df t.ratio p.value
# AP1 - AP8    -32.77 15.4 89  -2.122  0.0912
# AP1 - AP10   -34.62 15.4 89  -2.242  0.0698
# AP8 - AP10    -1.85 15.3 89  -0.121  0.9920

emmeans(Prot.C.AP.lm.int, pairwise~Site)
# contrast estimate   SE df t.ratio p.value
# BD - SS      56.3 15.4 89   3.644  0.0013 *
# BD - OL      28.3 15.4 89   1.831  0.1655
# SS - OL     -28.0 15.3 89  -1.828  0.1663

emmeans(Prot.C.AP.lm.int, pairwise~Season)
# contrast        estimate   SE df t.ratio p.value
# Fall - Winter      -45.4 17.7 89  -2.568  0.0567
# Fall - Spring       69.4 17.7 89   3.921  0.0010 *
# Fall - Summer       94.6 17.9 89   5.291  <.0001 *
# Winter - Spring    114.8 17.7 89   6.489  <.0001 *
# Winter - Summer    140.0 17.9 89   7.830  <.0001 *
# Spring - Summer     25.3 17.9 89   1.412  0.4951

##Plot Reaction Norms
Prot.C.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="TP_ug.cm2_C", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Prot.C.AP_Geno.Site.plot<-ggplot(Prot.C.AP_Geno.Site, aes(x=Site, y=TP_ug.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_C-ci, ymax=TP_ug.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Protein (\u03BCg cm'^-2*") Host")))+
  ylim(205, 620)+
  annotate("text", x = 1, y = 560, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 560, label = "B", size=levels.sz)+
  annotate("text", x = 3, y = 560, label = "AB", size=levels.sz)+
  annotate("text", x=0.85, y=265, label=expression(bolditalic(paste("Site p = 0.001, ", eta^2, " = 0.07"))), size=sig.sz, hjust = 0)+  
  annotate("text", x=0.85, y=240, label=expression(italic(paste("Genotype x Site p = 0.474, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)  

Prot.C.AP_Geno.Site.plot


Prot.C.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="TP_ug.cm2_C", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Prot.C.AP_Geno.Seas.plot<-ggplot(Prot.C.AP_Geno.Seas, aes(x=Season, y=TP_ug.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_C-ci, ymax=TP_ug.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Protein (\u03BCg cm'^-2*") Host")))+
  ylim(205, 620)+
  annotate("text", x = 1, y = 560, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 560, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 560, label = "B", size=levels.sz)+
  annotate("text", x = 4, y = 560, label = "B", size=levels.sz)+
  annotate("text", x=1, y=265, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.39"))), size=sig.sz, hjust = 0)+  
  annotate("text", x=1, y=240, label=expression(italic(paste("Genotype x Season p = 0.402, ", eta^2, " = 0.03"))), size=sig.sz, hjust = 0)

Prot.C.AP_Geno.Seas.plot


####Protein Symbiont Fraction

##Check normality
hist(CoralData.AP$TP_ug.cm2_S)
qqnorm(CoralData.AP$TP_ug.cm2_S)
shapiro.test(CoralData.AP$TP_ug.cm2_S)
#Normal

##Model with lm
Prot.S.AP.lm.int<-lm(TP_ug.cm2_S~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Prot.S.AP.lm.int), resid(Prot.S.AP.lm.int))
abline(0,0)

qqnorm(resid(Prot.S.AP.lm.int))
qqline(resid(Prot.S.AP.lm.int))

plot(density(resid(Prot.S.AP.lm.int)))

##Model results
summary(Prot.S.AP.lm.int)
# Multiple R-squared:  0.6064,	Adjusted R-squared:  0.5313 
# F-statistic: 8.067 on 17 and 89 DF,  p-value: 8.071e-12

anova(Prot.S.AP.lm.int)
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  24029   12014  6.4311  0.002462 ** 
# Site             2  15637    7819  4.1852  0.018316 *  
# Season           3 183476   61159 32.7373 2.363e-14 ***
# Genotype:Site    4  14152    3538  1.8938  0.118472    
# Genotype:Season  6  18900    3150  1.6861  0.133615    
# Residuals       89 166267    1868         

eta_squared(Prot.S.AP.lm.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.06 | [0.00, 1.00]
# Site            | 0.04 | [0.00, 1.00]
# Season          | 0.43 | [0.30, 1.00]
# Genotype:Site   | 0.03 | [0.00, 1.00]
# Genotype:Season | 0.04 | [0.00, 1.00]

emmeans(Prot.S.AP.lm.int, pairwise~Genotype)
# contrast   estimate   SE df t.ratio p.value
# AP1 - AP8     -35.3 10.3 89  -3.437  0.0026 *
# AP1 - AP10    -23.4 10.3 89  -2.274  0.0648
# AP8 - AP10     12.0 10.2 89   1.173  0.4721

emmeans(Prot.S.AP.lm.int, pairwise~Site)
# contrast estimate   SE df t.ratio p.value
# BD - SS     29.19 10.3 89   2.841  0.0152 *
# BD - OL     21.17 10.3 89   2.061  0.1039
# SS - OL     -8.01 10.2 89  -0.787  0.7122

emmeans(Prot.S.AP.lm.int, pairwise~Season)
# contrast        estimate   SE df t.ratio p.value
# Fall - Winter      -15.9 11.8 89  -1.355  0.5307
# Fall - Spring      -99.6 11.8 89  -8.470  <.0001 *
# Fall - Summer      -76.2 11.9 89  -6.406  <.0001 *
# Winter - Spring    -83.7 11.8 89  -7.115  <.0001 *
# Winter - Summer    -60.3 11.9 89  -5.066  <.0001 *
# Spring - Summer     23.4 11.9 89   1.971  0.2068

##Plot Reaction Norms
Prot.S.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="TP_ug.cm2_S", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Prot.S.AP_Geno.Site.plot<-ggplot(Prot.S.AP_Geno.Site, aes(x=Site, y=TP_ug.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_S-ci, ymax=TP_ug.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Protein (\u03BCg cm'^-2*") Symbiont")))+
  ylim(150, 400)+
  annotate("text", x = 1, y = 375, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 375, label = "B", size=levels.sz)+
  annotate("text", x = 3, y = 375, label = "AB", size=levels.sz)+
  annotate("text", x=0.85, y=185, label=expression(bolditalic(paste("Site p = 0.018, ", eta^2, " = 0.04"))), size=sig.sz, hjust = 0)+  
  annotate("text", x=0.85, y=170, label=expression(italic(paste("Genotype x Site p = 0.118, ", eta^2, " = 0.03"))), size=sig.sz, hjust = 0)

Prot.S.AP_Geno.Site.plot


Prot.S.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="TP_ug.cm2_S", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Prot.S.AP_Geno.Seas.plot<-ggplot(Prot.S.AP_Geno.Seas, aes(x=Season, y=TP_ug.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=TP_ug.cm2_S-ci, ymax=TP_ug.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Protein (\u03BCg cm'^-2*") Symbiont")))+
  ylim(150, 400)+
  annotate("text", x = 1, y = 375, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 375, label = "A", size=levels.sz)+
  annotate("text", x = 2.9, y = 375, label = "B", size=levels.sz)+
  annotate("text", x = 4, y = 375, label = "B", size=levels.sz)+
  annotate("text", x=4, y=185, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.43"))), size=sig.sz, hjust = 1)+  
  annotate("text", x=4, y=170, label=expression(italic(paste("Genotype x Season p = 0.133, ", eta^2, " = 0.04"))), size=sig.sz, hjust = 1)

Prot.S.AP_Geno.Seas.plot


####Physiology- Chlorophyll- A. palmata ####

####Chlorophyll-a

##Check normality
hist(CoralData.AP$Chl.a_ug.cm2)
qqnorm(CoralData.AP$Chl.a_ug.cm2)
shapiro.test(CoralData.AP$Chl.a_ug.cm2)
#Not normal

hist(log(CoralData.AP$Chl.a_ug.cm2+1))
qqnorm(log(CoralData.AP$Chl.a_ug.cm2+1))
shapiro.test(log(CoralData.AP$Chl.a_ug.cm2+1))
#Still Not Normal but less skewed

##Model with log+1 transformation
Chl.a.AP.log.int<-lm(log(Chl.a_ug.cm2+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Chl.a.AP.log.int), resid(Chl.a.AP.log.int))
abline(0,0)

qqnorm(resid(Chl.a.AP.log.int))
qqline(resid(Chl.a.AP.log.int))

plot(density(resid(Chl.a.AP.log.int)))

##Model results
summary(Chl.a.AP.log.int)
# Multiple R-squared:  0.7173,	Adjusted R-squared:  0.6633 
# F-statistic: 13.28 on 17 and 89 DF,  p-value: < 2.2e-16

anova(Chl.a.AP.log.int)
#                   Df Sum Sq Mean Sq  F value    Pr(>F)    
# Genotype         2 0.0479 0.02394  1.2198    0.3002    
# Site             2 0.5791 0.28954 14.7508 2.932e-06 ***
# Season           3 3.5076 1.16921 59.5658 < 2.2e-16 ***
# Genotype:Site    4 0.1019 0.02549  1.2984    0.2768    
# Genotype:Season  6 0.1957 0.03261  1.6613    0.1399    
# Residuals       89 1.7470 0.01963

eta_squared(Chl.a.AP.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 7.75e-03 | [0.00, 1.00]
# Site            |     0.09 | [0.01, 1.00]
# Season          |     0.57 | [0.45, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]
# Genotype:Season |     0.03 | [0.00, 1.00]

emmeans(Chl.a.AP.log.int, pairwise~Site)
# contrast estimate     SE df t.ratio p.value
# BD - SS    0.0776 0.0333 89   2.330  0.0568 *
# BD - OL   -0.1017 0.0333 89  -3.054  0.0083 *
# SS - OL   -0.1792 0.0330 89  -5.428  <.0001 *

emmeans(Chl.a.AP.log.int, pairwise~Season)
# contrast        estimate     SE df t.ratio p.value
# Fall - Winter    -0.1863 0.0381 89  -4.887  <.0001 *
# Fall - Spring     0.2507 0.0381 89   6.574  <.0001 *
# Fall - Summer     0.2389 0.0386 89   6.198  <.0001 *
# Winter - Spring   0.4370 0.0381 89  11.461  <.0001 *
# Winter - Summer   0.4253 0.0386 89  11.031  <.0001 *
# Spring - Summer  -0.0118 0.0386 89  -0.305  0.9901

##Plot Reaction Norms
Chl.a.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="Chl.a_ug.cm2", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Chl.a.AP_Geno.Site.plot<-ggplot(Chl.a.AP_Geno.Site, aes(x=Site, y=Chl.a_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.a_ug.cm2-ci, ymax=Chl.a_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")))+
  ylim(0.5, 3.3)+
  annotate("text", x = 1, y = 2.8, label = "A", size=levels.sz)+
  annotate("text", x = 2, y = 2.8, label = "A", size=levels.sz)+
  annotate("text", x = 3, y = 2.8, label = "B", size=levels.sz)+
  annotate("text", x = 0.85, y = 0.9, label=expression(bolditalic(paste("Site p < 0.001, ", eta^2, " = 0.09"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 0.85, y = 0.75, label=expression(italic(paste("Genotype x Site p = 0.277, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)

Chl.a.AP_Geno.Site.plot


Chl.a.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="Chl.a_ug.cm2", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Chl.a.AP_Geno.Seas.plot<-ggplot(Chl.a.AP_Geno.Seas, aes(x=Season, y=Chl.a_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.a_ug.cm2-ci, ymax=Chl.a_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")))+
  ylim(0.5, 3.3)+
  annotate("text", x = 1, y = 2.8, label = "A", size=levels.sz)+
  annotate("text", x = 1.85, y = 2.8, label = "B", size=levels.sz)+
  annotate("text", x = 3, y = 2.8, label = "C", size=levels.sz)+
  annotate("text", x = 4, y = 2.8, label = "C", size=levels.sz)+
  annotate("text", x = 1, y = 0.9, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.57"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 1, y = 0.75, label=expression(italic(paste("Genotype x Season p = 0.134, ", eta^2, " = 0.03"))), size=sig.sz, hjust = 0)

Chl.a.AP_Geno.Seas.plot


####Chlorophyll-c2

##Check normality
hist(CoralData.AP$Chl.c2_ug.cm2)
qqnorm(CoralData.AP$Chl.c2_ug.cm2)
shapiro.test(CoralData.AP$Chl.c2_ug.cm2)
#Not normal

hist(log(CoralData.AP$Chl.c2_ug.cm2+1))
qqnorm(log(CoralData.AP$Chl.c2_ug.cm2+1))
shapiro.test(log(CoralData.AP$Chl.c2_ug.cm2+1))
#Normal

##Model with log+1 transformation 
Chl.c2.AP.log.int<-lm(log(Chl.c2_ug.cm2+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Chl.c2.AP.log.int), resid(Chl.c2.AP.log.int))
abline(0,0)

qqnorm(resid(Chl.c2.AP.log.int))
qqline(resid(Chl.c2.AP.log.int))

plot(density(resid(Chl.c2.AP.log.int)))

##Model results
summary(Chl.c2.AP.log.int)
# Multiple R-squared:  0.639,	Adjusted R-squared:   0.57 
# F-statistic: 9.267 on 17 and 89 DF,  p-value: 2.521e-13

anova(Chl.c2.AP.log.int)
#                 Df Sum Sq Mean Sq F value  Pr(>F)    
# Genotype         2 0.0818 0.04089  0.6710 0.51375    
# Site             2 0.2261 0.11303  1.8550 0.16245    
# Season           3 8.3266 2.77553 45.5489 < 2e-16 ***
# Genotype:Site    4 0.2538 0.06345  1.0413 0.39052    
# Genotype:Season  6 0.7112 0.11854  1.9453 0.08211 .  
# Residuals       89 5.4232 0.06094  

eta_squared(Chl.c2.AP.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 5.44e-03 | [0.00, 1.00]
# Site            |     0.02 | [0.00, 1.00]
# Season          |     0.55 | [0.44, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]
# Genotype:Season |     0.05 | [0.00, 1.00]

emmeans(Chl.c2.AP.log.int, pairwise~Season)
# contrast        estimate     SE df t.ratio p.value
# Fall - Winter    -0.1787 0.0672 89  -2.660  0.0450 *
# Fall - Spring     0.4317 0.0672 89   6.426  <.0001 *
# Fall - Summer     0.4799 0.0679 89   7.065  <.0001 *
# Winter - Spring   0.6104 0.0672 89   9.086  <.0001 *
# Winter - Summer   0.6586 0.0679 89   9.696  <.0001 *
# Spring - Summer   0.0482 0.0679 89   0.709  0.8931

##Plot Reaction Norms
Chl.c2.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="Chl.c2_ug.cm2", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Chl.c2.AP_Geno.Site.plot<-ggplot(Chl.c2.AP_Geno.Site, aes(x=Site, y=Chl.c2_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.c2_ug.cm2-ci, ymax=Chl.c2_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Chlorophyll-c2 (\u03BCg cm'^-2*")")))+
  ylim(0, 3.8)+
  annotate("text", x = 0.85, y = 0.5, label=expression(italic(paste("Site p = 0.162, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 0.85, y = 0.3, label=expression(italic(paste("Genotype x Site p = 0.391, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)

Chl.c2.AP_Geno.Site.plot


Chl.c2.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="Chl.c2_ug.cm2", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Chl.c2.AP_Geno.Seas.plot<-ggplot(Chl.c2.AP_Geno.Seas, aes(x=Season, y=Chl.c2_ug.cm2, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=Chl.c2_ug.cm2-ci, ymax=Chl.c2_ug.cm2+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Chlorophyll-c2 (\u03BCg cm'^-2*")")))+
  ylim(0, 3.8)+
  annotate("text", x = 1, y = 3, label = "A", size=levels.sz)+
  annotate("text", x = 1.85, y = 3, label = "B", size=levels.sz)+
  annotate("text", x = 3, y = 3, label = "C", size=levels.sz)+
  annotate("text", x = 4, y = 3, label = "C", size=levels.sz)+
  annotate("text", x = 1, y = 0.5, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.55"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 1, y = 0.3, label=expression(italic(paste("Genotype x Season p = 0.082, ", eta^2, " = 0.05"))), size=sig.sz, hjust = 0)

Chl.c2.AP_Geno.Seas.plot


####Physiology- Biomass- A. palmata ####

####Biomass Host Fraction

##Check normality
hist(CoralData.AP$AFDW_mg.cm2_C)
qqnorm(CoralData.AP$AFDW_mg.cm2_C)
shapiro.test(CoralData.AP$AFDW_mg.cm2_C)
#Normal

##Model with lm
Bio.C.AP.lm.int<-lm(AFDW_mg.cm2_C~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Bio.C.AP.lm.int), resid(Bio.C.AP.lm.int))
abline(0,0)

qqnorm(resid(Bio.C.AP.lm.int))
qqline(resid(Bio.C.AP.lm.int))

plot(density(resid(Bio.C.AP.lm.int)))

##Model results
summary(Bio.C.AP.lm.int)
# Multiple R-squared:  0.3333,	Adjusted R-squared:  0.206 
# F-statistic: 2.618 on 17 and 89 DF,  p-value: 0.001801

anova(Bio.C.AP.lm.int)
#                   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.2892 0.14461  3.0246 0.0535974 .  
# Site             2 0.1135 0.05675  1.1870 0.3099245    
# Season           3 0.9920 0.33068  6.9162 0.0003075 ***
# Genotype:Site    4 0.3373 0.08433  1.7637 0.1432193    
# Genotype:Season  6 0.3957 0.06594  1.3792 0.2317234    
# Residuals       89 4.2552 0.04781 

eta_squared(Bio.C.AP.lm.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        | 0.05 | [0.00, 1.00]
# Site            | 0.02 | [0.00, 1.00]
# Season          | 0.16 | [0.04, 1.00]
# Genotype:Site   | 0.05 | [0.00, 1.00]
# Genotype:Season | 0.06 | [0.00, 1.00]

emmeans(Bio.C.AP.lm.int, pairwise~Season)
# contrast        estimate     SE df t.ratio p.value
# Fall - Winter   -0.15605 0.0595 89  -2.622  0.0495 *
# Fall - Spring   -0.00745 0.0595 89  -0.125  0.9993
# Fall - Summer    0.11973 0.0602 89   1.990  0.1996
# Winter - Spring  0.14860 0.0595 89   2.497  0.0673
# Winter - Summer  0.27578 0.0602 89   4.583  0.0001 *
# Spring - Summer  0.12718 0.0602 89   2.114  0.1567

##Plot Reaction Norms
Bio.C.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="AFDW_mg.cm2_C", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Bio.C.AP_Geno.Site.plot<-ggplot(Bio.C.AP_Geno.Site, aes(x=Site, y=AFDW_mg.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_C-ci, ymax=AFDW_mg.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Biomass (AFDW mg cm'^-2*") Host")))+
  ylim(0.65,1.65)+
  annotate("text", x = 0.85, y = 0.8, label=expression(italic(paste("Site p = 0.310, ", eta^2, " = 0.02"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 0.85, y = 0.73, label=expression(italic(paste("Genotype x Site p = 0.143, ", eta^2, " = 0.05"))), size=sig.sz, hjust = 0)  

Bio.C.AP_Geno.Site.plot


Bio.C.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="AFDW_mg.cm2_C", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Bio.C.AP_Geno.Seas.plot<-ggplot(Bio.C.AP_Geno.Seas, aes(x=Season, y=AFDW_mg.cm2_C, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_C-ci, ymax=AFDW_mg.cm2_C+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Biomass (AFDW mg cm'^-2*") Host")))+
  ylim(0.65,1.65)+
  annotate("text", x = 1, y = 1.5, label = "A", size=levels.sz)+
  annotate("text", x = 1.9, y = 1.5, label = "B", size=levels.sz)+
  annotate("text", x = 3, y = 1.5, label = "AB", size=levels.sz)+
  annotate("text", x = 4, y = 1.5, label = "A", size=levels.sz)+
  annotate("text", x = 1, y = 0.8, label=expression(bolditalic(paste("Season p < 0.001, ", eta^2, " = 0.16"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 1, y = 0.73, label=expression(italic(paste("Genotype x Season p = 0.232, ", eta^2, " = 0.06"))), size=sig.sz, hjust = 0)


Bio.C.AP_Geno.Seas.plot


####Biomass Symbiont Fraction

##Check normality
hist(CoralData.AP$AFDW_mg.cm2_S)
qqnorm(CoralData.AP$AFDW_mg.cm2_S)
shapiro.test(CoralData.AP$AFDW_mg.cm2_S)
#Not Normal

hist(log(CoralData.AP$AFDW_mg.cm2_S+1))
qqnorm(log(CoralData.AP$AFDW_mg.cm2_S+1))
shapiro.test(log(CoralData.AP$AFDW_mg.cm2_S+1))
#Still Not Normal but less skewed

##Model with log+1 transformation 
Bio.S.AP.log.int<-lm(log(AFDW_mg.cm2_S+1)~Genotype+Site+Season+Genotype:Site+Genotype:Season, data=CoralData.AP)

##Check residuals
plot(fitted(Bio.S.AP.log.int), resid(Bio.S.AP.log.int))
abline(0,0)

qqnorm(resid(Bio.S.AP.log.int))
qqline(resid(Bio.S.AP.log.int))

plot(density(resid(Bio.S.AP.log.int)))

##Model results
summary(Bio.S.AP.log.int)
# Multiple R-squared:  0.2455,	Adjusted R-squared:  0.1014 
# F-statistic: 1.703 on 17 and 89 DF,  p-value: 0.05665

anova(Bio.S.AP.log.int)
#                   Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.11545 0.057726  6.4088 0.00251 **
# Site             2 0.04474 0.022371  2.4836 0.08921 . 
# Season           3 0.00681 0.002269  0.2519 0.85976   
# Genotype:Site    4 0.06309 0.015772  1.7510 0.14589   
# Genotype:Season  6 0.03073 0.005122  0.5687 0.75423   
# Residuals       89 0.80165 0.009007     

eta_squared(Bio.S.AP.log.int, partial=FALSE)
# Parameter       | Eta2 |       95% CI
# Genotype        |     0.11 | [0.02, 1.00]
# Site            |     0.04 | [0.00, 1.00]
# Season          | 6.41e-03 | [0.00, 1.00]
# Genotype:Site   |     0.06 | [0.00, 1.00]
# Genotype:Season |     0.03 | [0.00, 1.00]]

emmeans(Bio.S.AP.log.int, pairwise~Genotype)
# contrast   estimate     SE df t.ratio p.value
# AP1 - AP8   -0.0787 0.0226 89  -3.487  0.0022 *
# AP1 - AP10  -0.0547 0.0226 89  -2.424  0.0453
# AP8 - AP10   0.0240 0.0224 89   1.072  0.5340

##Plot Reaction Norms
Bio.S.AP_Geno.Site<-summarySE(CoralData.AP, measurevar="AFDW_mg.cm2_S", groupvars=c("Species", "Genotype", "Site"), na.rm=TRUE)

Bio.S.AP_Geno.Site.plot<-ggplot(Bio.S.AP_Geno.Site, aes(x=Site, y=AFDW_mg.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_S-ci, ymax=AFDW_mg.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Site", y=expression(paste('Biomass (AFDW mg cm'^-2*") Symbiont")))+
  ylim(0.2, 0.9)+
  annotate("text", x = 0.85, y = 0.3, label=expression(italic(paste("Site p = 0.089, ", eta^2, " =  0.04"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 0.85, y = 0.26, label=expression(italic(paste("Genotype x Site p = 0.146, ", eta^2, " =  0.06"))), size=sig.sz, hjust = 0)

Bio.S.AP_Geno.Site.plot


Bio.S.AP_Geno.Seas<-summarySE(CoralData.AP, measurevar="AFDW_mg.cm2_S", groupvars=c("Species", "Genotype", "Season"), na.rm=TRUE)

Bio.S.AP_Geno.Seas.plot<-ggplot(Bio.S.AP_Geno.Seas, aes(x=Season, y=AFDW_mg.cm2_S, colour=Genotype, group=Genotype)) + 
  scale_colour_manual(values=AP.genotype.colors.o)+
  geom_errorbar(aes(ymin=AFDW_mg.cm2_S-ci, ymax=AFDW_mg.cm2_S+ci), width=0.2, position=position_dodge(0.2), size=0.6) +
  geom_line(position=position_dodge(0.2), size=1) +
  geom_point(position=position_dodge(0.2), size=3)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .05), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  labs(x="Season", y=expression(paste('Biomass (AFDW mg cm'^-2*") Symbiont")))+
  ylim(0.2, 0.9)+
  annotate("text", x = 1, y = 0.3, label=expression(italic(paste("Season p = 0.860, ", eta^2, " < 0.01"))), size=sig.sz, hjust = 0)+  
  annotate("text", x = 1, y = 0.26, label=expression(italic(paste("Genotype x Season p = 0.754, ", eta^2, " = 0.03"))), size=sig.sz, hjust = 0)

Bio.S.AP_Geno.Seas.plot


#####***Figure S6, S7, S8: Plot Physiology Univariate Panels- A. palmata ####

##Biomass
Phys.Bio.AP_fig<-plot_grid(Bio.C.AP_Geno.Seas.plot, Bio.C.AP_Geno.Site.plot,
                           Bio.S.AP_Geno.Seas.plot, Bio.S.AP_Geno.Site.plot, 
                           rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                           ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                           label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS6_Apal_Physiology_ReactionNorms_Bio.png", plot=Phys.Bio.AP_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS6_Apal_Physiology_ReactionNorms_Bio.pdf", plot=Phys.Bio.AP_fig, dpi=300, width=12, height=12, units="in")


##Protein
Phys.Prot.AP_fig<-plot_grid(Prot.C.AP_Geno.Seas.plot, Prot.C.AP_Geno.Site.plot,
                            Prot.S.AP_Geno.Seas.plot, Prot.S.AP_Geno.Site.plot, 
                            rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                            ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                            label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS7_Apal_Physiology_ReactionNorms_Prot.png", plot=Phys.Prot.AP_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS7_Apal_Physiology_ReactionNorms_Prot.pdf", plot=Phys.Prot.AP_fig, dpi=300, width=12, height=12, units="in")


##Chlorophyll
Phys.Chl.AP_fig<-plot_grid(Chl.a.AP_Geno.Seas.plot, Chl.a.AP_Geno.Site.plot,
                           Chl.c2.AP_Geno.Seas.plot, Chl.c2.AP_Geno.Site.plot, 
                           rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                           ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                           label_size=panel.lab.sz, align="hv")
ggsave(filename="FigS8_Apal_Physiology_ReactionNorms_Chl.png", plot=Phys.Chl.AP_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="FigS8_Apal_Physiology_ReactionNorms_Chl.pdf", plot=Phys.Chl.AP_fig, dpi=300, width=12, height=12, units="in")


####**Physiology Part 3: Multivariate Analysis####

##Remove NAs
CoralData.rm <- na.omit(CoralData)

##Log transform variables for analysis
names(CoralData.rm)
CoralData.log <- log(CoralData.rm[,10:15]+1)

CoralData.log.full <-CoralData.rm
CoralData.log.full[,10:15]<-log(CoralData.log.full[,10:15]+1)


####Physiology Permanova Analysis- A. cervicornis####
CoralData.log.full.AC<-subset(CoralData.log.full, Species=="Acer")
rownames(CoralData.log.full.AC)<-paste("AC", CoralData.log.full.AC$RandN, sep="")

names(CoralData.log.full.AC)
adonis(vegdist(CoralData.log.full.AC[,10:15], "euclidean")~CoralData.log.full.AC$Genotype+CoralData.log.full.AC$Site+CoralData.log.full.AC$Season+
         CoralData.log.full.AC$Genotype:CoralData.log.full.AC$Site+ CoralData.log.full.AC$Genotype:CoralData.log.full.AC$Season, data=CoralData.log.full.AC, method="euclidean")
#                                                               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# CoralData.log.full.AC$Genotype                                3     3.404  1.1348   8.570 0.05542  0.001 ***
# CoralData.log.full.AC$Site                                    2    11.151  5.5754  42.106 0.18153  0.001 ***
# CoralData.log.full.AC$Season                                  3    16.280  5.4268  40.985 0.26504  0.001 ***
# CoralData.log.full.AC$Genotype:CoralData.log.full.AC$Site     6     1.340  0.2233   1.687 0.02182  0.041 *  
# CoralData.log.full.AC$Genotype:CoralData.log.full.AC$Season   9     1.842  0.2046   1.545 0.02998  0.039 *  
# Residuals                                                   207    27.409  0.1324         0.44621           
# Total                                                       230    61.426                 1.00000 

####Genotype Comparisons

##Check dispersion by Genotype
Phys.Geno.disp.AC<-betadisper(vegdist(CoralData.log.full.AC[,10:15], "euclidean"), CoralData.log.full.AC$Genotype)
anova(Phys.Geno.disp.AC)
#             Df Sum Sq  Mean Sq F value  Pr(>F)  
# Groups      3 0.3485 0.116157  3.1987 0.02421 *
# Residuals 227 8.2433 0.036314    
#Significant difference in dispersion 

TukeyHSD(Phys.Geno.disp.AC)
#                   diff          lwr         upr     p adj
# AC10-AC8   0.05658980 -0.033837970 0.147017568 0.3695972
# AC12-AC8  -0.02538566 -0.117400710 0.066629390 0.8914578
# AC15-AC8   0.06923508 -0.022779967 0.161250133 0.2113642
# AC12-AC10 -0.08197546 -0.173616354 0.009665436 0.0975844
# AC15-AC10  0.01264528 -0.078995612 0.104286179 0.9843578
# AC15-AC12  0.09462074  0.001413223 0.187828262 0.0451267 *

##Pairwise comparisons between Genotypes
pairwise.adonis(x=CoralData.log.full.AC[,10:15],factors=CoralData.log.full.AC$Genotype,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')
#           pairs Df SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1 AC10 vs AC12  1 1.7209900 6.9395292 0.057380157   0.001      0.006   *
# 2 AC10 vs AC15  1 1.6063810 5.3384455 0.044733661   0.005      0.030   .
# 3  AC10 vs AC8  1 0.1784113 0.6890473 0.005854812   0.550      1.000    
# 4 AC12 vs AC15  1 0.3951816 1.5677374 0.014051888   0.201      1.000    
# 5  AC12 vs AC8  1 1.6654474 7.9345460 0.065610251   0.001      0.006   *
# 6  AC15 vs AC8  1 1.2468637 4.7359494 0.040225177   0.006      0.036   .

####Site Comparisons

##Check dispersion by Site
Phys.Site.disp.AC<-betadisper(vegdist(CoralData.log.full.AC[,10:15], "euclidean"), CoralData.log.full.AC$Site)
anova(Phys.Site.disp.AC)
#             Df  Sum Sq Mean Sq F value   Pr(>F)   
# Groups      2 0.0170 0.0084973  0.3001  0.741
# Residuals 228 6.4551 0.0283118        
#No significant difference in dispersion 

##Pairwise comparisons between Site
pairwise.adonis(x=CoralData.log.full.AC[,10:15],factors=CoralData.log.full.AC$Site,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')

#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 BD vs KL  1 9.2244303 40.988562 0.21574226   0.001      0.003   *
# 2 BD vs SS  1 0.3679373  1.676848 0.01112877   0.153      0.459    
# 3 KL vs SS  1 6.5952131 29.815714 0.15874984   0.001      0.003   *

####Season Comparisons

##Check dispersion by Season
Phys.Season.disp.AC<-betadisper(vegdist(CoralData.log.full.AC[,10:15], "euclidean"), CoralData.log.full.AC$Season)
anova(Phys.Season.disp.AC)
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups      3 0.2169 0.072295  1.7861 0.1506
# Residuals 227 9.1883 0.040477       
#No significant difference in dispersion 

##Pairwise comparisons between Seasons
pairwise.adonis(x=CoralData.log.full.AC[,10:15],factors=CoralData.log.full.AC$Season,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')

#               pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1   Fall vs Winter  1 0.3079583  1.352679 0.01225779   0.228      1.000    
# 2   Fall vs Spring  1 8.7894796 41.590910 0.26225280   0.001      0.006   *
# 3   Fall vs Summer  1 5.3732719 24.318449 0.17208262   0.001      0.006   *
# 4 Winter vs Spring  1 9.3325295 51.598902 0.31930231   0.001      0.006   *
# 5 Winter vs Summer  1 6.1912282 32.397509 0.22751458   0.001      0.006   *
# 6 Spring vs Summer  1 1.3640987  7.690019 0.06118241   0.001      0.006   *


####Physiology Permanova Analysis- A. palmata####
CoralData.log.full.AP<-subset(CoralData.log.full, Species=="Apal")
rownames(CoralData.log.full.AP)<-paste("A", CoralData.log.full.AP$RandN, sep="")

names(CoralData.log.full.AP)
adonis(vegdist(CoralData.log.full.AP[,10:15], "euclidean")~CoralData.log.full.AP$Genotype+CoralData.log.full.AP$Site+CoralData.log.full.AP$Season+
         CoralData.log.full.AP$Genotype:CoralData.log.full.AP$Site+ CoralData.log.full.AP$Genotype:CoralData.log.full.AP$Season, data=CoralData.log.full.AP, method="euclidean")
#                                                               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# CoralData.log.full.AP$Genotype                                2     0.877  0.4384   2.786 0.02470  0.020 *  
# CoralData.log.full.AP$Site                                    2     1.514  0.7570   4.810 0.04265  0.002 ** 
# CoralData.log.full.AP$Season                                  3    16.751  5.5838  35.479 0.47189  0.001 ***
# CoralData.log.full.AP$Genotype:CoralData.log.full.AP$Site     4     0.875  0.2188   1.390 0.02465  0.182    
# CoralData.log.full.AP$Genotype:CoralData.log.full.AP$Season   6     1.474  0.2457   1.561 0.04152  0.095 .  
# Residuals                                                    89    14.007  0.1574         0.39458           
# Total                                                       106    35.499                 1.00000       

####Genotype Comparisons

##Check dispersion by Genotype
Phys.Geno.disp.AP<-betadisper(vegdist(CoralData.log.full.AP[,10:15], "euclidean"), CoralData.log.full.AP$Genotype)
anova(Phys.Geno.disp.AP)
#             Df Sum Sq  Mean Sq F value  Pr(>F)  
# Groups      2 0.0943 0.047128  1.0324 0.3598
# Residuals 104 4.7476 0.045650 
#No significant difference in dispersion 

##Pairwise comparisons between Genotypes
pairwise.adonis(x=CoralData.log.full.AP[,10:15],factors=CoralData.log.full.AP$Genotype,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')
#           pairs Df SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1 AP1 vs AP10  1 0.4921708 1.3849116 0.019676256   0.230      0.690    
# 2  AP1 vs AP8  1 0.6247566 1.8881854 0.026636108   0.151      0.453    
# 3 AP10 vs AP8  1 0.2017438 0.6450936 0.009131471   0.564      1.000  

####Site Comparisons

##Check dispersion by Site
Phys.Site.disp.AP<-betadisper(vegdist(CoralData.log.full.AP[,10:15], "euclidean"), CoralData.log.full.AP$Site)
anova(Phys.Site.disp.AP)
#             Df Sum Sq  Mean Sq F value Pr(>F)
# Groups      2 0.1463 0.073137  1.5043  0.227
# Residuals 104 5.0564 0.048619     
#No significant difference in dispersion 

##Pairwise comparisons between Site
pairwise.adonis(x=CoralData.log.full.AP[,10:15],factors=CoralData.log.full.AP$Site,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')
#       pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 BD vs OL  1 0.5100391 1.477613 0.02096571   0.216      0.648    
# 2 BD vs SS  1 1.0633550 3.652290 0.05027081   0.028      0.084    
# 3 OL vs SS  1 0.7343067 2.140450 0.02967060   0.119      0.357  

####Season Comparisons

##Check dispersion by Season
Phys.Season.disp.AP<-betadisper(vegdist(CoralData.log.full.AP[,10:15], "euclidean"), CoralData.log.full.AP$Season)
anova(Phys.Season.disp.AP)
#           Df Sum Sq  Mean Sq F value Pr(>F)
# Groups      3 0.0154 0.0051327  0.1745 0.9135
# Residuals 103 3.0301 0.0294189        
#No significant difference in dispersion 

##Pairwise comparisons between Seasons
pairwise.adonis(x=CoralData.log.full.AP[,10:15],factors=CoralData.log.full.AP$Season,sim.function="vegdist", sim.method='euclidean',p.adjust.m='bonferroni')

#               pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1   Fall vs Summer  1  5.8025966 30.546127 0.37458709   0.001      0.006   *
# 2   Fall vs Winter  1  1.1634353  6.265667 0.10753618   0.002      0.012   .
# 3   Fall vs Spring  1  5.6240223 32.955425 0.38791431   0.001      0.006   *
# 4 Summer vs Winter  1 10.6678947 55.043789 0.51906660   0.001      0.006   *
# 5 Summer vs Spring  1  0.2196331  1.230548 0.02355993   0.300      1.000    
# 6 Winter vs Spring  1 10.0020096 57.342056 0.52442819   0.001      0.006   *



####**Physiology Part 4: Variance Partitioning####


####Variance Partitioning of Physiology by Genotype and Environment- A. cervicornis####
Phys.vp.AC<-varpart(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Genotype, CoralData.log.full.AC$Site, CoralData.log.full.AC$Season)
Phys.vp.AC$part
# Individual fractions  Df           Adj.R.square    Testable                                  
# [a] = X1 | X2+X3       3               0.05987     TRUE
# [b] = X2 | X1+X3       2               0.18817     TRUE
# [c] = X3 | X1+X2       3               0.26405     TRUE
# [h] = Residuals                        0.51595    FALSE

##Check model significance
anova(rda(CoralData.log.full.AC[,10:15]~ CoralData.log.full.AC$Genotype + CoralData.log.full.AC$Site + CoralData.log.full.AC$Season))
#The model is significant p=0.001

##Check variance explained by model
RsquareAdj(rda(CoralData.log.full.AC[,10:15]~ CoralData.log.full.AC$Genotype + CoralData.log.full.AC$Site + CoralData.log.full.AC$Season))$adj.r.squared*100
#The model explains 48.40% of the variation in coral physiology in A cerv

##Check effect of Genotype 
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Genotype))
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Genotype, CoralData.log.full.AC[,c(5:6)]))
#Genotype with overlap and Genotype only (controlling for other variables) both significant p=0.001 

##Check effect of Site 
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Site))
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Site,  CoralData.log.full.AC[,c(4,6)]))
#Site with overlap and Site only (controlling for other variables) both significant p=0.001 

##Check effect of Season
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Season))
anova(rda(CoralData.log.full.AC[,10:15], CoralData.log.full.AC$Season,  CoralData.log.full.AC[,c(4:5)]))
#Season with overlap and Season only (controlling for other variables) both significant p=0.001 


####Plot Physiology Variance Partitioning- A. cervicornis####
Phys.vp.AC.Fracs<-data.frame(Fractions=c("Genotype", "Site", "Season", "Geno&Site", "Site&Seas", "Geno&Seas", "All", "Residuals"), 
                             Variance=c(Phys.vp.AC$part$indfract$Adj.R.square))
Phys.vp.AC.Fracs$Percent<-round(Phys.vp.AC.Fracs$Variance*100,2)
Phys.vp.AC.Fracs.cut<-Phys.vp.AC.Fracs[which(Phys.vp.AC.Fracs$Fractions=="Genotype" |
                                               Phys.vp.AC.Fracs$Fractions=="Site" |
                                               Phys.vp.AC.Fracs$Fractions=="Season" |
                                               Phys.vp.AC.Fracs$Fractions=="Residuals"),]
Phys.vp.AC.Fracs.cut$Fractions<-factor(Phys.vp.AC.Fracs.cut$Fractions, levels=c("Genotype", "Site", "Season", "Residuals"), ordered=TRUE)

Phys.vp.AC.plot<-ggplot(Phys.vp.AC.Fracs.cut, aes(x=Fractions, y=Percent, fill=Fractions)) + 
  geom_bar(stat="identity", show.legend=FALSE)+
  scale_fill_manual(values = c(GE.colors.o, "darkgray"))+
  theme_classic()+
  theme(axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.title.sz, colour="black"), 
        axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position = "none")+
  ylab("Percent Variance Explained")+
  xlab(NULL)+
  ylim(0,60)+
  geom_text(aes(label=paste0(sprintf("%1.1f",Percent),"%"), fontface = "bold"),
            vjust=-1.75, size=8)+
  geom_text(aes(label=c("p = 0.001", "p = 0.001", "p = 0.001", ""), fontface = "italic"),
            vjust=-.5, size=sig.sz)

Phys.vp.AC.plot


####Variance Partitioning of Physiology by Genotype and Environment- A. palmata####
Phys.vp.AP<-varpart(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Genotype, CoralData.log.full.AP$Site, CoralData.log.full.AP$Season)
Phys.vp.AP$part
# Individual fractions Df           Adj.R.square    Testable                                   
# [a] = X1 | X2+X3       2               0.01509     TRUE
# [b] = X2 | X1+X3       2               0.03565     TRUE
# [c] = X3 | X1+X2       3               0.47588     TRUE
# [h] = Residuals                        0.49334    FALSE

##Check model significance
anova(rda(CoralData.log.full.AP[,10:15]~ CoralData.log.full.AP$Genotype + CoralData.log.full.AP$Site + CoralData.log.full.AP$Season))
#The model is significant p=0.001

##Check variance explained by model
RsquareAdj(rda(CoralData.log.full.AP[,10:15]~ CoralData.log.full.AP$Genotype + CoralData.log.full.AP$Site + CoralData.log.full.AP$Season))$adj.r.squared*100
#The model explains 50.67% of the variation in coral physiology in A pal

##Check effect of Genotype 
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Genotype))
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Genotype, CoralData.log.full.AP[,c(5:6)]))
#Genotype with overlap not significant p=0.267
#Genotype only (controlling for other variables) significant p=0.031

##Check effect of Site
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Site))
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Site,  CoralData.log.full.AP[,c(4,6)]))
#Site with overlap and Site only (controlling for other variables) both significant p=0.045, 0.001, resp.

##Check effect of Season
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Season))
anova(rda(CoralData.log.full.AP[,10:15], CoralData.log.full.AP$Season,  CoralData.log.full.AP[,c(4:5)]))
#Season with overlap and Season only (controlling for other variables) both significant p=0.001 


####Plot Physiology Variance Partitioning- A. palmata####
Phys.vp.AP.Fracs<-data.frame(Fractions=c("Genotype", "Site", "Season", "Geno&Site", "Site&Seas", "Geno&Seas", "All", "Residuals"), 
                             Variance=c(Phys.vp.AP$part$indfract$Adj.R.square))
Phys.vp.AP.Fracs$Percent<-round(Phys.vp.AP.Fracs$Variance*100,2)
Phys.vp.AP.Fracs.cut<-Phys.vp.AP.Fracs[which(Phys.vp.AP.Fracs$Fractions=="Genotype" |
                                               Phys.vp.AP.Fracs$Fractions=="Site" |
                                               Phys.vp.AP.Fracs$Fractions=="Season" |
                                               Phys.vp.AP.Fracs$Fractions=="Residuals"),]
Phys.vp.AP.Fracs.cut$Fractions<-factor(Phys.vp.AP.Fracs.cut$Fractions, levels=c("Genotype", "Site", "Season", "Residuals"), ordered=TRUE)

Phys.vp.AP.plot<-ggplot(Phys.vp.AP.Fracs.cut, aes(x=Fractions, y=Percent, fill=Fractions)) + 
  geom_bar(stat="identity", show.legend=FALSE)+
  scale_fill_manual(values = c(GE.colors.o, "darkgray"))+
  theme_classic()+
  theme(axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.title.sz, colour="black"), 
        axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position = "none")+
  ylab("Percent Variance Explained")+
  xlab(NULL)+
  ylim(0,60)+
  geom_text(aes(label=paste0(sprintf("%1.1f",Percent),"%"), fontface = "bold"),
            vjust=-1.75, size=8)+
  geom_text(aes(label=c("p = 0.031", "p = 0.001", "p = 0.001", ""), fontface = "italic"),
            vjust=-.5, size=sig.sz)

Phys.vp.AP.plot


####***Figure 2: Plot Physiology Variance Partitioning Panels- A. cervicornis ####
Phys.VarPart_fig<-plot_grid(Phys.vp.AC.plot, Phys.vp.AP.plot,
                            rel_widths=c(1, 1), rel_heights=c(1, 1), 
                            ncol=2, nrow=1, byrow=T, labels = c('A', 'B'), 
                            label_size=20, align="v")
ggsave(filename="Fig2_Acer_and_Apal_Physiology_VarPart.png", plot=Phys.VarPart_fig, dpi=300, width=12, height=6, units="in")
ggsave(filename="Fig2_Acer_and_Apal_Physiology_VarPart.pdf", plot=Phys.VarPart_fig, dpi=300, width=12, height=6, units="in")



#-----*DNA Methylation------------------------------------------------------------


####**Methylation Part 1: DNA Methylation States (MSAP)####

####Load MSAP Fragment Analysis Data####

##Load all data
AC_C1Frag<-read.csv("ACER_C1_FAM_Blue.csv", header=TRUE)
AC_C2Frag<-read.csv("ACER_C2_HEX_Green.csv", header=TRUE)
AC_C4Frag<-read.csv("ACER_C4_HEX_Green.csv", header=TRUE)

AP_C1Frag<-read.csv("APAL_C1_FAM_Blue.csv", header=TRUE)
AP_C2Frag<-read.csv("APAL_C2_HEX_Green.csv", header=TRUE)
AP_C4Frag<-read.csv("APAL_C4_HEX_Green.csv", header=TRUE)

FragMeta<-read.csv("FragMeta.csv", header=TRUE)
PlateMeta<-read.csv("PlateMeta.csv", header=TRUE)


####Organize and Filter Fragment Data-Acropora cervicornis####

##Merge each Primer file with Fragment Meta Data by FileName
#Adds Well ID
AC_C1<-merge(FragMeta, AC_C1Frag, all.x=FALSE)
AC_C2<-merge(FragMeta, AC_C2Frag, all.x=FALSE)
AC_C4<-merge(FragMeta, AC_C4Frag, all.x=FALSE)

##Merge each with Plate Meta Data by WellID
#Adds Random Number
AC_C1<-merge(PlateMeta, AC_C1, all.x=FALSE)
AC_C2<-merge(PlateMeta, AC_C2, all.x=FALSE)
AC_C4<-merge(PlateMeta, AC_C4, all.x=FALSE)

##Merge each with Coral Meta Data
AC_C1<-merge(CoralMeta, AC_C1, all.x=FALSE)
AC_C2<-merge(CoralMeta, AC_C2, all.x=FALSE)
AC_C4<-merge(CoralMeta, AC_C4, all.x=FALSE)

####Filter each Primer for Loci in at least 15% of samples 
nrow(AC_C1) #480 samples in A. cervicornis 

####Primer C1

##Create a new dataframe by transposing columns containing loci
names(AC_C1)
ncol(AC_C1)
AC_C1.t<-data.frame(t(AC_C1[, c(21:247)]))
rownames(AC_C1.t)
names(AC_C1.t)<-AC_C1$Input

##Sum number of peaks for each loci
AC_C1.t$L.Peaks<-apply(AC_C1.t, 1, sum)
boxplot(AC_C1.t$L.Peaks)

#Identify loci present in at least 15% of samples
C1.AC.Loci.15<-rownames(AC_C1.t)[which(AC_C1.t$L.Peaks>(0.15*480))]
length(C1.AC.Loci.15)
#82 Loci

##Create a dataframe retaining only filtered loci
names(AC_C1)
AC_C1_F15<- cbind(AC_C1[,c(1:20)], AC_C1[, names(AC_C1) %in% C1.AC.Loci.15])


####Primer C2

##Create a new dataframe by transposing columns containing loci
names(AC_C2)
ncol(AC_C2)
AC_C2.t<-data.frame(t(AC_C2[, c(21:91)]))
rownames(AC_C2.t)
names(AC_C2.t)<-AC_C2$Input

##Sum number of peaks for each loci
AC_C2.t$L.Peaks<-apply(AC_C2.t, 1, sum)
boxplot(AC_C2.t$L.Peaks)

#Identify loci present in at least 15% of samples
C2.AC.Loci.15<-rownames(AC_C2.t)[which(AC_C2.t$L.Peaks>(0.15*480))]
length(C2.AC.Loci.15)
#9 Loci

##Create a dataframe retaining only filtered loci
names(AC_C2)
AC_C2_F15<- cbind(AC_C2[,c(1:20)], AC_C2[, names(AC_C2) %in% C2.AC.Loci.15])


####Primer C4

##Create a new dataframe by transposing columns containing loci
names(AC_C4)
ncol(AC_C4)
AC_C4.t<-data.frame(t(AC_C4[, c(21:68)]))
rownames(AC_C4.t)
names(AC_C4.t)<-AC_C4$Input

##Sum number of peaks for each loci
AC_C4.t$L.Peaks<-apply(AC_C4.t, 1, sum)
boxplot(AC_C4.t$L.Peaks)

#Identify loci present in at least 15% of samples
C4.AC.Loci.15<-rownames(AC_C4.t)[which(AC_C4.t$L.Peaks>(0.15*480))]
length(C4.AC.Loci.15)
#8 Loci

##Create a dataframe retaining only filtered loci
names(AC_C4)
AC_C4_F15<- cbind(AC_C4[,c(1:20)], AC_C4[, names(AC_C4) %in% C4.AC.Loci.15])

##Merge Primer Datasets together and Organize for msap function 
#First column: Genotype (initial grouping variable)
#Second column: Sample ID
#Third column: Enzyme
names(AC_C1_F15)
ncol(AC_C1_F15)
ncol(AC_C2_F15)
AC_Frags<-merge(AC_C1_F15[,c(4, 2, 17, 21:102)], AC_C2_F15[,c(4, 2, 17, 21:29)])
ncol(AC_C4_F15)
AC_Frags<-merge(AC_Frags, AC_C4_F15[,c(4, 2, 17, 21:28)])

##Prepare organized .csv file for MSAP function
write.csv(AC_Frags, "AC_Frags.csv", row.names = FALSE)


####Methylation State Data- A. cervicornis####

##Check number of loci per primer
length(C1.AC.Loci.15)
length(C2.AC.Loci.15)
length(C4.AC.Loci.15)

##Run MSAP function
AC_MSAP<-msap("AC_Frags.csv", no.bands="h", loci.per.primer=c(82, 9, 8))

# Primer:  1 
# --Number of Methylation-Susceptible Loci (MSL):  82 
# --Number of No Methylated Loci (NML):  0 
# 
# Primer:  2 
# --Number of Methylation-Susceptible Loci (MSL):  9 
# --Number of No Methylated Loci (NML):  0 
# 
# Primer:  3 
# --Number of Methylation-Susceptible Loci (MSL):  8 
# --Number of No Methylated Loci (NML):  0 
# 
# All combinations: 
#   Number of Methylation-Susceptible Loci (MSL):  99 
# Number of No Methylated Loci (NML):  0 
# 
# Number of polymorphic MSL:  91  ( 92 % of total MSL)

# Analysis of MSL
# Report of methylation levels 
#                                                     AC10   AC12   AC15    AC8
# HPA+/MSP+ (Unmethylated)                          0.2195 0.2264 0.2397 0.2286
# HPA+/MSP- (Hemimethylated)                        0.1887 0.1739 0.1785 0.1739
# HPA-/MSP+ (Internal cytosine methylation)         0.1343 0.1172 0.1234 0.1160
# HPA-/MSP- (Full methylation or absence of target) 0.4574 0.4825 0.4584 0.4815
# 
# 
# Performing AMOVA
# AMOVA TABLE 	d.f. 	SSD 		MSD 		Variance
# among groups	 3 	 102 	 34 	 0.3426 
# within groups	 236 	 3172 	 13.44 	 13.44 
# Total        	 239 	 3274 	 13.7  
# 
# Phi_ST =  0.02486   (P<0.0001)

##Extract binary polymorphic MSL 
AC_MSAP.b<-AC_MSAP$transformed.MSL
names(AC_MSAP.b)[c(1:2)]<-c("Genotype" , "ID")
#1: Unmethylated
#2: Methylated

##Create dataframe with methylation patterns by individual for MSL (Methylation-Susceptible Loci)
names(AC_MSAP$patterns)
AC_MetPat<-data.frame(rbind(AC_MSAP$patterns$AC10, AC_MSAP$patterns$AC12, 
                            AC_MSAP$patterns$AC15, AC_MSAP$patterns$AC8))

##Add a column with Genotype and Sample ID
AC_MetPat<-cbind(AC_MSAP.b[,c(1:2)], AC_MetPat)
str(AC_MetPat)

##Change loci character variables to factors 
#Levels= u: Unmethylated, h: Hemimethylated, i: Internal C, f: Full methylation or Lack of target
names(AC_MetPat)
ncol(AC_MetPat)
AC_MetPat[,3:101] <- lapply(AC_MetPat[,3:101], function(x) {factor(x, levels=c("u", "h", "i", "f"))})
str(AC_MetPat)

##Retain only polymorphic MSL for analysis 
poly.ac<-c(names(AC_MSAP.b))
AC_MetPat.p<-AC_MetPat[,names(AC_MetPat) %in% poly.ac]

##Merge Metadata with methylation pattern data
MetData.AC<-merge(CoralMeta, AC_MetPat.p)

##Add rownames of Sample ID
rownames(MetData.AC)<- MetData.AC$ID


####Organize and Filter Fragment Data-Acropora palmata####

##Merge each Primer file with Fragment Meta Data by FileName
#Adds Well ID
AP_C1<-merge(FragMeta, AP_C1Frag, all.x=FALSE)
AP_C2<-merge(FragMeta, AP_C2Frag, all.x=FALSE)
AP_C4<-merge(FragMeta, AP_C4Frag, all.x=FALSE)

##Merge each with Plate Meta Data by WellID
#Adds Random Number
AP_C1<-merge(PlateMeta, AP_C1, all.x=FALSE)
AP_C2<-merge(PlateMeta, AP_C2, all.x=FALSE)
AP_C4<-merge(PlateMeta, AP_C4, all.x=FALSE)

##Merge each with Coral Meta Data
AP_C1<-merge(CoralMeta, AP_C1, all.x=FALSE)
AP_C2<-merge(CoralMeta, AP_C2, all.x=FALSE)
AP_C4<-merge(CoralMeta, AP_C4, all.x=FALSE)


####Filter each Primer for Loci in at least 15% of samples 
nrow(AP_C1) #214 samples in A. palmata 

####Primer C1

##Create a new dataframe by transposing columns containing loci
names(AP_C1)
ncol(AP_C1)
AP_C1.t<-data.frame(t(AP_C1[, c(21:111)]))
rownames(AP_C1.t)
names(AP_C1.t)<-AP_C1$Input

##Sum number of peaks for each loci
AP_C1.t$L.Peaks<-apply(AP_C1.t, 1, sum)
boxplot(AP_C1.t$L.Peaks)

#Identify loci present in at least 15% of samples
C1.AP.Loci.15<-rownames(AP_C1.t)[which(AP_C1.t$L.Peaks>(0.15*214))]
length(C1.AP.Loci.15)
#38 Loci

##Create a dataframe retaining only filtered loci
names(AP_C1)
AP_C1_F15<- cbind(AP_C1[,c(1:20)], AP_C1[, names(AP_C1) %in% C1.AP.Loci.15])


####Primer C2

##Create a new dataframe by transposing columns containing loci
names(AP_C2)
ncol(AP_C2)
AP_C2.t<-data.frame(t(AP_C2[, c(21:70)]))
rownames(AP_C2.t)
names(AP_C2.t)<-AP_C2$Input

##Sum number of peaks for each loci
AP_C2.t$L.Peaks<-apply(AP_C2.t, 1, sum)
boxplot(AP_C2.t$L.Peaks)

#Identify loci present in at least 15% of samples
C2.AP.Loci.15<-rownames(AP_C2.t)[which(AP_C2.t$L.Peaks>(0.15*214))]
length(C2.AP.Loci.15)
#13 Loci

##Create a dataframe retaining only filtered loci
names(AP_C2)
AP_C2_F15<- cbind(AP_C2[,c(1:20)], AP_C2[, names(AP_C2) %in% C2.AP.Loci.15])


####Primer C4

##Create a new dataframe by transposing columns containing loci
names(AP_C4)
ncol(AP_C4)
AP_C4.t<-data.frame(t(AP_C4[, c(21:104)]))
rownames(AP_C4.t)
names(AP_C4.t)<-AP_C4$Input

##Sum number of peaks for each loci
AP_C4.t$L.Peaks<-apply(AP_C4.t, 1, sum)
boxplot(AP_C4.t$L.Peaks)

#Identify loci present in at least 15% of samples
C4.AP.Loci.15<-rownames(AP_C4.t)[which(AP_C4.t$L.Peaks>(0.15*214))]
length(C4.AP.Loci.15)
#24 Loci

##Create a dataframe retaining only filtered loci
names(AP_C4)
AP_C4_F15<- cbind(AP_C4[,c(1:20)], AP_C4[, names(AP_C4) %in% C4.AP.Loci.15])


##Merge Primer Datasets together and Organize for msap function 
#First column: Genotype (initial grouping variable)
#Second column: Sample ID
#Third column: Enzyme
names(AP_C1_F15)
ncol(AP_C1_F15)
ncol(AP_C2_F15)
AP_Frags<-merge(AP_C1_F15[,c(4, 2, 17, 21:58)], AP_C2_F15[,c(4, 2, 17, 21:33)])
ncol(AP_C4_F15)
AP_Frags<-merge(AP_Frags, AP_C4_F15[,c(4, 2, 17, 21:44)])

##Prepare organized .csv file for MSAP function
write.csv(AP_Frags, "AP_Frags.csv", row.names = FALSE)


####Methylation State Data- A. palmata####

##Check number of loci per primer
length(C1.AP.Loci.15)
length(C2.AP.Loci.15)
length(C4.AP.Loci.15)

##Run MSAP function
AP_MSAP<-msap("AP_Frags.csv", no.bands="h", loci.per.primer=c(38, 13, 24))

# Primer:  1 
# --Number of Methylation-Susceptible Loci (MSL):  38 
# --Number of No Methylated Loci (NML):  0 
# 
# Primer:  2 
# --Number of Methylation-Susceptible Loci (MSL):  12 
# --Number of No Methylated Loci (NML):  1 
# 
# Primer:  3 
# --Number of Methylation-Susceptible Loci (MSL):  23 
# --Number of No Methylated Loci (NML):  1 
# 
# All combinations: 
#   Number of Methylation-Susceptible Loci (MSL):  73 
# Number of No Methylated Loci (NML):  2 
# 
# Number of polymorphic MSL:  68  ( 93 % of total MSL)
# Number of polymorphic NML:  2  ( 100 % of total NML)
# 
# Analysis of MSL
# Report of methylation levels 
#                                                      AP1   AP10    AP8
# HPA+/MSP+ (Unmethylated)                          0.2481 0.2622 0.2785
# HPA+/MSP- (Hemimethylated)                        0.2121 0.2036 0.2230
# HPA-/MSP+ (Internal cytosine methylation)         0.1264 0.1008 0.1244
# HPA-/MSP- (Full methylation or absence of target) 0.4133 0.4334 0.3740
# 
# Performing AMOVA
# AMOVA TABLE 	d.f. 	SSD 		MSD 		Variance
# among groups	 2 	 111.7 	 55.85 	 1.314 
# within groups	 104 	 934.4 	 8.985 	 8.985 
# Total        	 106 	 1046 	 9.869  
# 
# Phi_ST =  0.1276   (P<0.0001) 


##Extract binary polymorphic MSL 
AP_MSAP.b<-AP_MSAP$transformed.MSL
names(AP_MSAP.b)[c(1:2)]<-c("Genotype" , "ID")
#1: Unmethylated
#2: Methylated

##Create dataframe with methylation patterns by individual for MSL (Methylation-Susceptible Loci)
names(AP_MSAP$patterns)
AP_MetPat<-data.frame(rbind(AP_MSAP$patterns$AP1, AP_MSAP$patterns$AP10, AP_MSAP$patterns$AP8))

##Add a column with Genotype and Sample ID
AP_MetPat<-cbind(AP_MSAP.b[,c(1:2)], AP_MetPat)
str(AP_MetPat)

##Change loci character variables to factors 
#Levels= u: Unmethylated, h: Hemimethylated, i: Internal C, f: Full methylation or Lack of target
names(AP_MetPat)
ncol(AP_MetPat)
AP_MetPat[,3:75] <- lapply(AP_MetPat[,3:75], function(x) {factor(x, levels=c("u", "h", "i", "f"))})
str(AP_MetPat)

##Retain only polymorphic MSL for analysis 
poly.ap<-c(names(AP_MSAP.b))
AP_MetPat.p<-AP_MetPat[,names(AP_MetPat) %in% poly.ap]

##Merge Metadata with methylation pattern data
MetData.AP<-merge(CoralMeta, AP_MetPat.p)

##Add rownames of Sample ID
rownames(MetData.AP)<- MetData.AP$ID



####Percent of Methylation States- A. cervicornis####

##Create a new dataframe for summary of Methylation Types
names(MetData.AC)
Met.Types.AC<-MetData.AC[,c(1:7)]

##Sum of each methylation state across each sample
#NMT: Nonmethylated
#HMM: Hemimethylated
#ICM: Internal C Methylation
#HPM: Hypermethylation 
Met.Types.AC$NMT<-rowSums(MetData.AC[,-c(1:7)]=="u")
Met.Types.AC$HMM<-rowSums(MetData.AC[,-c(1:7)]=="h")
Met.Types.AC$ICM<-rowSums(MetData.AC[,-c(1:7)]=="i")
Met.Types.AC$HPM<-rowSums(MetData.AC[,-c(1:7)]=="f")

##Check number of polymorphic MSL
length(poly.ac[-c(1:2)])
#91

##Calculate percent of each methylation state (out of 91 loci)
Met.Types.AC$NMT.p<-(Met.Types.AC$NMT/91)*100
Met.Types.AC$HMM.p<-(Met.Types.AC$HMM/91)*100
Met.Types.AC$ICM.p<-(Met.Types.AC$ICM/91)*100
Met.Types.AC$HPM.p<-(Met.Types.AC$HPM/91)*100

##Average percent of each methylation state by Genotype
names(Met.Types.AC)
Met.Types.AC.pa<-aggregate(Met.Types.AC[12:15], list(Met.Types.AC$Genotype), FUN="mean")
names(Met.Types.AC.pa)<-c("Genotype", "NMT", "HMM", "ICM", "HPM")
Met.Types.AC.pa$Genotype<-factor(Met.Types.AC.pa$Genotype, levels=c("AC8", "AC10", "AC12", "AC15"), ordered=TRUE)

##Add cumulative sums for stacked error bars
Met.Types.AC.pa$NMT_stack<-Met.Types.AC.pa$HPM+Met.Types.AC.pa$ICM+Met.Types.AC.pa$HMM+Met.Types.AC.pa$NMT
Met.Types.AC.pa$HMM_stack<-Met.Types.AC.pa$HPM+Met.Types.AC.pa$ICM+Met.Types.AC.pa$HMM
Met.Types.AC.pa$ICM_stack<-Met.Types.AC.pa$HPM+Met.Types.AC.pa$ICM
Met.Types.AC.pa$HPM_stack<-Met.Types.AC.pa$HPM

##Add standard error
names(Met.Types.AC)
Met.Types.AC.pse<-aggregate(Met.Types.AC[12:15], list(Met.Types.AC$Genotype), FUN="std.error")
names(Met.Types.AC.pse)<-c("Genotype", "NMT", "HMM", "ICM", "HPM")
Met.Types.AC.pse$Genotype<-factor(Met.Types.AC.pse$Genotype, levels=c("AC8", "AC10", "AC12", "AC15"), ordered=TRUE)

##Change to long format
Met.Types.AC_long.p <- melt(Met.Types.AC.pa, 
                            id.vars=c("Genotype"), 
                            measure.vars=c("NMT", "HMM", "ICM", "HPM"),
                            variable.name="Type", value.name="Percent")

Met.Types.AC_long.se <- melt(Met.Types.AC.pse, 
                             id.vars=c("Genotype"), 
                             measure.vars=c("NMT", "HMM", "ICM", "HPM"),
                             variable.name="Type.se", value.name="SE")

Met.Types.AC_long.p$SE<-Met.Types.AC_long.se$SE

Met.Types.AC_long.p.stack <- melt(Met.Types.AC.pa, 
                                  id.vars=c("Genotype"), 
                                  measure.vars=c("NMT_stack", "HMM_stack", "ICM_stack", "HPM_stack"),
                                  variable.name="Type", value.name="Percent.Stack")

Met.Types.AC_long.p$Percent.Stack<-Met.Types.AC_long.p.stack$Percent.Stack


####Plot Percent of Methylation States by Genotype- A. cervicornis####
Met.Type.AC.plot<-ggplot(Met.Types.AC_long.p, aes(x=Genotype, y=Percent, fill=Type)) + 
  geom_bar(position="stack", stat="identity")+
  geom_errorbar(aes(ymin=Percent.Stack-SE, ymax=Percent.Stack+SE), width=0.2, size=0.6) +
  scale_fill_manual(values =Met.colors.o)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz))+
  geom_text(aes(label=sprintf("%1.1f",Percent), fontface = "bold"),
            position=position_stack(vjust=0.5), size=levels.sz)+
  theme(plot.title=element_text(hjust=0.5))

Met.Type.AC.plot


####Percent of Methylation States- A. palmata####

##Create a new dataframe for summary of Methylation Types
names(MetData.AP)
Met.Types.AP<-MetData.AP[,c(1:7)]

#Sum of each methylation state across each sample
#NMT: Nonmethylated
#HMM: Hemimethylated
#ICM: Internal C Methylation
#HPM: Hypermethylation 
Met.Types.AP$NMT<-rowSums(MetData.AP[,-c(1:7)]=="u")
Met.Types.AP$HMM<-rowSums(MetData.AP[,-c(1:7)]=="h")
Met.Types.AP$ICM<-rowSums(MetData.AP[,-c(1:7)]=="i")
Met.Types.AP$HPM<-rowSums(MetData.AP[,-c(1:7)]=="f")

##Check number of polymorphic MSL
length(poly.ap[-c(1:2)])
#68

##Calculate percent of each methylation state (out of 68 loci)
Met.Types.AP$NMT.p<-(Met.Types.AP$NMT/68)*100
Met.Types.AP$HMM.p<-(Met.Types.AP$HMM/68)*100
Met.Types.AP$ICM.p<-(Met.Types.AP$ICM/68)*100
Met.Types.AP$HPM.p<-(Met.Types.AP$HPM/68)*100

##Average percent of each methylation state by Genotype
names(Met.Types.AP)
Met.Types.AP.pa<-aggregate(Met.Types.AP[12:15], list(Met.Types.AP$Genotype), FUN="mean")
names(Met.Types.AP.pa)<-c("Genotype", "NMT", "HMM", "ICM", "HPM")
Met.Types.AP.pa$Genotype<-factor(Met.Types.AP.pa$Genotype, levels=c("AP1", "AP8", "AP10"), ordered=TRUE)

##Add cumulative sums for stacked error bars
Met.Types.AP.pa$NMT_stack<-Met.Types.AP.pa$HPM+Met.Types.AP.pa$ICM+Met.Types.AP.pa$HMM+Met.Types.AP.pa$NMT
Met.Types.AP.pa$HMM_stack<-Met.Types.AP.pa$HPM+Met.Types.AP.pa$ICM+Met.Types.AP.pa$HMM
Met.Types.AP.pa$ICM_stack<-Met.Types.AP.pa$HPM+Met.Types.AP.pa$ICM
Met.Types.AP.pa$HPM_stack<-Met.Types.AP.pa$HPM

##Add standard error
names(Met.Types.AP)
Met.Types.AP.pse<-aggregate(Met.Types.AP[12:15], list(Met.Types.AP$Genotype), FUN="std.error")
names(Met.Types.AP.pse)<-c("Genotype", "NMT", "HMM", "ICM", "HPM")
Met.Types.AP.pse$Genotype<-factor(Met.Types.AP.pse$Genotype, levels=c("AP1", "AP8", "AP10"), ordered=TRUE)

##Change to long format
Met.Types.AP_long.p <- melt(Met.Types.AP.pa, 
                            id.vars=c("Genotype"), 
                            measure.vars=c("NMT", "HMM", "ICM", "HPM"),
                            variable.name="Type", value.name="Percent")

Met.Types.AP_long.se <- melt(Met.Types.AP.pse, 
                             id.vars=c("Genotype"), 
                             measure.vars=c("NMT", "HMM", "ICM", "HPM"),
                             variable.name="Type.se", value.name="SE")

Met.Types.AP_long.p$SE<-Met.Types.AP_long.se$SE

Met.Types.AP_long.p.stack <- melt(Met.Types.AP.pa, 
                                  id.vars=c("Genotype"), 
                                  measure.vars=c("NMT_stack", "HMM_stack", "ICM_stack", "HPM_stack"),
                                  variable.name="Type", value.name="Percent.Stack")

Met.Types.AP_long.p$Percent.Stack<-Met.Types.AP_long.p.stack$Percent.Stack


####Plot Percent of Methylation States by Genotype- A. palmata####
Met.Type.AP.plot<-ggplot(Met.Types.AP_long.p, aes(x=Genotype, y=Percent, fill=Type)) + 
  geom_bar(position="stack", stat="identity")+
  geom_errorbar(aes(ymin=Percent.Stack-SE, ymax=Percent.Stack+SE), width=0.2, size=0.6) +
  scale_fill_manual(values =Met.colors.o)+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz))+
  geom_text(aes(label=sprintf("%1.1f",Percent), fontface = "bold"),
            position=position_stack(vjust=0.5), size=levels.sz)+
  theme(plot.title=element_text(hjust=0.5))

Met.Type.AP.plot


####**Methylation Part 2: DNA Methylation Permanova####


####DNA Methylation Permanova Analysis- A. cervicornis####
names(MetData.AC)
adonis(daisy(MetData.AC[, -c(1:7)], "gower")~MetData.AC$Genotype+MetData.AC$Site+MetData.AC$Season+
         MetData.AC$Genotype:MetData.AC$Site+ MetData.AC$Genotype:MetData.AC$Season, data=MetData.AC, method="gower")
#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# MetData.AC$Genotype                     3     1.738 0.57919  3.8884 0.04272  0.001 ***
# MetData.AC$Site                         2     0.276 0.13812  0.9272 0.00679  0.498    
# MetData.AC$Season                       3     3.514 1.17135  7.8638 0.08639  0.001 ***
# MetData.AC$Genotype:MetData.AC$Site     6     1.268 0.21126  1.4183 0.03116  0.039 *  
# MetData.AC$Genotype:MetData.AC$Season   9     1.706 0.18960  1.2729 0.04195  0.071 .  
# Residuals                             216    32.174 0.14895         0.79099           
# Total                                 239    40.676                 1.00000     


####Genotype Comparisons

##Check dispersion by Genotype
Met.Geno.disp.AC<-betadisper(daisy(MetData.AC[, -c(1:7)], "gower"), MetData.AC$Genotype)
anova(Met.Geno.disp.AC)
#             Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      3 0.02831 0.0094365  2.0637 0.1057
# Residuals 236 1.07916 0.0045727  
#No significant difference in dispersion

##Pairwise comparisons between Genotypes
pairwise.adonis(x=MetData.AC[, -c(1:7)],factors=MetData.AC$Genotype,sim.function="daisy", sim.method='gower',p.adjust.m='bonferroni')
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 AC10 vs AC12  1 0.7543433 4.718347 0.03844859   0.001      0.006   *
# 2 AC10 vs AC15  1 0.4965846 2.939908 0.02430883   0.007      0.042   .
# 3  AC10 vs AC8  1 0.7470746 4.703954 0.03833580   0.002      0.012   .
# 4 AC12 vs AC15  1 0.4527583 2.645108 0.02192470   0.015      0.090    
# 5  AC12 vs AC8  1 0.3903534 2.423427 0.02012421   0.029      0.174    
# 6  AC15 vs AC8  1 0.6340539 3.727270 0.03061984   0.002      0.012   .


####Season Comparisons

##Check dispersion by Season
Met.Season.disp.AC<-betadisper(daisy(MetData.AC[, -c(1:7)], "gower"), MetData.AC$Season)
anova(Met.Season.disp.AC)
# Df  Sum Sq  Mean Sq F value   Pr(>F)    
# Groups      3 0.13540 0.045132  11.145 7.23e-07 ***
# Residuals 236 0.95571 0.004050       
#Significant difference in dispersion 

TukeyHSD(Met.Season.disp.AC)
#                       diff         lwr          upr     p adj
# Winter-Fall    0.005995457 -0.02406627  0.036057188 0.9551986
# Spring-Fall   -0.051604273 -0.08166600 -0.021542542 0.0000806 *
# Summer-Fall   -0.033498660 -0.06356039 -0.003436930 0.0221890 *
# Spring-Winter -0.057599730 -0.08766146 -0.027538000 0.0000081 *
# Summer-Winter -0.039494118 -0.06955585 -0.009432387 0.0043827 *
# Summer-Spring  0.018105612 -0.01195612  0.048167343 0.4045382

##Pairwise comparisons between Seasons
pairwise.adonis(x=MetData.AC[, -c(1:7)],factors=MetData.AC$Season,sim.function="daisy", sim.method='gower',p.adjust.m='bonferroni')
# pairs Df SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1   Fall vs Winter  1 0.2643552  1.521254 0.012727893   0.130      0.780    
# 2   Fall vs Spring  1 1.4964084  9.714286 0.076062640   0.001      0.006   *
# 3   Fall vs Summer  1 1.5623687  9.804852 0.076717372   0.001      0.006   *
# 4 Winter vs Spring  1 1.8025571 11.585657 0.089405395   0.001      0.006   *
# 5 Winter vs Summer  1 1.7457151 10.850387 0.084209194   0.001      0.006   *
# 6 Spring vs Summer  1 0.1567202  1.110253 0.009321222   0.310      1.000    


####DNA Methylation Permanova Analysis- A. palmata####
names(MetData.AP)
adonis(daisy(MetData.AP[, -c(1:7)], "gower")~MetData.AP$Genotype+MetData.AP$Site+MetData.AP$Season+
         MetData.AP$Genotype:MetData.AP$Site+MetData.AP$Genotype:MetData.AP$Season, data=MetData.AP, method="gower")
#                                         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# MetData.AP$Genotype                     2    2.6462 1.32308 11.6011 0.16471  0.001 ***
# MetData.AP$Site                         2    0.3877 0.19384  1.6996 0.02413  0.055 .  
# MetData.AP$Season                       3    1.8527 0.61758  5.4151 0.11533  0.001 ***
# MetData.AP$Genotype:MetData.AP$Site     4    0.5644 0.14111  1.2373 0.03513  0.171    
# MetData.AP$Genotype:MetData.AP$Season   6    0.4640 0.07734  0.6781 0.02889  0.964    
# Residuals                              89   10.1502 0.11405         0.63181           
# Total                                 106   16.0653                 1.00000        


####Genotype Comparisons

##Check dispersion by Genotype
Met.Geno.disp.AP<-betadisper(daisy(MetData.AP[, -c(1:7)], "gower"), MetData.AP$Genotype)
anova(Met.Geno.disp.AP)
#             Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      2 0.01255 0.0062746   1.593 0.2082
# Residuals 104 0.40963 0.0039387 
#No significant difference in dispersion 

##Pairwise comparisons between Genotypes
pairwise.adonis(x=MetData.AP[, -c(1:7)],factors=MetData.AP$Genotype,sim.function="daisy", sim.method='gower',p.adjust.m='bonferroni')
#         pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 AP1 vs AP10  1 0.8118454  6.059056 0.08072385   0.001      0.003   *
# 2  AP1 vs AP8  1 1.5408252 12.100885 0.14920781   0.001      0.003   *
# 3 AP10 vs AP8  1 1.6124838 12.816245 0.15475520   0.001      0.003   *


####Season Comparisons

##Check dispersion by Season
Met.Season.disp.AP<-betadisper(daisy(MetData.AP[, -c(1:7)], "gower"), MetData.AP$Season)
anova(Met.Season.disp.AP)
#             Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups      3 0.00402 0.0013389  0.4213 0.7381
# Residuals 103 0.32735 0.0031781 
#No significant difference in dispersion

##Pairwise comparisons between Seasons
pairwise.adonis(x=MetData.AP[, -c(1:7)],factors=MetData.AP$Season,sim.function="daisy", sim.method='gower',p.adjust.m='bonferroni')
#               pairs Df  SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1   Fall vs Summer  1 0.99275467 7.3087467 0.125345632   0.001      0.006   *
# 2   Fall vs Winter  1 0.04611207 0.3416025 0.006526404   0.959      1.000    
# 3   Fall vs Spring  1 0.86260893 6.1618208 0.105942708   0.001      0.006   *
# 4 Summer vs Winter  1 0.92672873 6.8199138 0.117950950   0.001      0.006   *
# 5 Summer vs Spring  1 0.11085602 0.7862751 0.015183079   0.620      1.000    
# 6 Winter vs Spring  1 0.77511854 5.5347368 0.096198178   0.001      0.006   *



####**Methylation Part 3: Variance Partitioning####


####Variance Partitioning of Methylation by Genotype and Season- A. cervicornis####
Met.vp.AC<-varpart(daisy(MetData.AC[, -c(1:7)], "gower"), MetData.AC$Genotype, MetData.AC$Season)
Met.vp.AC$part
# Individual fractions  Df R.squared Adj.R.squared Testable                                 
# [a] = X1|X2           3                 0.03190     TRUE
# [b]                   0                -0.00136    FALSE
# [c] = X2|X1           3                 0.07613     TRUE
# [d] = Residuals                         0.89332    FALSE

##Check model significance
met.rda.AC<-dbrda(daisy(MetData.AC[, -c(1:7)]) ~ Genotype + Season, data=MetData.AC, dist="gower")
anova(met.rda.AC)
#The model is significant p=0.001

##Check variance explained by model
RsquareAdj(met.rda.AC)$adj.r.squared*100
#The model explains 10.67% of the variation in coral DNA methylation in A cerv

##Check effect of Genotype 
anova(dbrda(daisy(MetData.AC[, -c(1:7)])~Genotype, data=MetData.AC, dist="gower"))
anova(dbrda(daisy(MetData.AC[, -c(1:7)])~Genotype+Condition(Season), data=MetData.AC, dist="gower"))
#Genotype with overlap and Genotype only (controlling for Season) both significant p=0.001 

##Check effect of Season
anova(dbrda(daisy(MetData.AC[, -c(1:7)])~Season, data=MetData.AC, dist="gower"))
anova(dbrda(daisy(MetData.AC[, -c(1:7)])~Season+Condition(Genotype), data=MetData.AC, dist="gower"))
#Season with overlap and Season only (controlling for Genotype) both significant p=0.001 


####Variance Partitioning of Methylation by Genotype and Season- A. palmata####
Met.vp.AP<-varpart(daisy(MetData.AP[, -c(1:7)], "gower"), MetData.AP$Genotype, MetData.AP$Season)
Met.vp.AP$part
# Individual fractions Df R.squared Adj.R.squared Testable
# [a] = X1|X2           2                 0.15501     TRUE
# [b]                   0                -0.00636    FALSE
# [c] = X2|X1           3                 0.09609     TRUE
# [d] = Residuals                         0.75526    FALSE

##Check model significance
met.rda.AP<-dbrda(daisy(MetData.AP[, -c(1:7)]) ~ Genotype + Season, data=MetData.AP, dist="gower")
anova.cca(met.rda.AP)
#The model is significant p=0.001

##Check variance explained by model
RsquareAdj(met.rda.AP)$adj.r.squared*100
#The model explains 24.47% of the variation in coral DNA methylation in A pal

##Check effect of Genotype 
anova.cca(dbrda(daisy(MetData.AP[, -c(1:7)])~Genotype, data=MetData.AP, dist="gower"))
anova.cca(dbrda(daisy(MetData.AP[, -c(1:7)])~Genotype+Condition(Season), data=MetData.AP, dist="gower"))
#Genotype with overlap and Genotype only (controlling for Season) both significant p=0.001 

##Check effect of Season
anova.cca(dbrda(daisy(MetData.AP[, -c(1:7)])~Season, data=MetData.AP, dist="gower"))
anova.cca(dbrda(daisy(MetData.AP[, -c(1:7)])~Season+Condition(Genotype), data=MetData.AP, dist="gower"))
#Season with overlap and Season only (controlling for Genotype) both significant p=0.001 



####**Methylation Part 4: Discriminant Analysis of Principle Components####


####Methylation DAPC- A. cervicornis####

##Create new dataframe to group by Significant Factors (Genotype and Season)
MetData.AC.set<-MetData.AC

##Create a new variable of Seasonal Period
#Cool: "Cooling" time points of Fall and Winter 
#Warm: "Warming" time points of Spring and Summer 
MetData.AC.set$Period<-"Cooling"
MetData.AC.set$Period[which(MetData.AC.set$Season=="Spring" | MetData.AC.set$Season=="Summer")]<-"Warming"
MetData.AC.set$Period<-factor(MetData.AC.set$Period, levels=c("Cooling", "Warming"), ordered=TRUE)

##Create a new variable of Set by Genotype and Period 
MetData.AC.set$Set<-paste(MetData.AC.set$Genotype, MetData.AC.set$Period, sep=".")
MetData.AC.set$Set<-factor(MetData.AC.set$Set, levels=c("AC8.Cooling", "AC8.Warming", "AC10.Cooling", "AC10.Warming",
                                                        "AC12.Cooling", "AC12.Warming", "AC15.Cooling", "AC15.Warming"), ordered=TRUE)

##Prepare Input Data
names(MetData.AC.set)
Gen_Data.AC <- df2genind(MetData.AC.set[,c(8:98)], ploidy=1, ncode=1)
pop.ord.AC<-factor(MetData.AC.set$Set, ordered=FALSE)

##Number of PCA's to retain as k-1
n.pca.AC<-length(levels(pop.ord.AC))-1

##DAPC Analysis with Parameters Obtained
MetDAPC.AC <- dapc(Gen_Data.AC, pop=pop.ord.AC, n.pca=n.pca.AC, n.da=2, scale=FALSE,
                   var.contrib=TRUE)
summary(MetDAPC.AC)

##Calculate percent of variance explained by each axis
#LD1
(MetDAPC.AC$eig[1]/sum(MetDAPC.AC$eig))*100
#63.83%

#LD2
(MetDAPC.AC$eig[2]/sum(MetDAPC.AC$eig))*100
#30.25%

##Prepare for Plotting
met.dapc.plot.data.AC<-data.frame(MetDAPC.AC$ind.coord)
met.dapc.plot.data.AC$ID<-rownames(MetDAPC.AC$ind.coord)
met.dapc.plot.data.AC<-merge(met.dapc.plot.data.AC, MetData.AC.set)  


####Plot DAPC Scatter Plot of DNA Methylation- A. cervicornis####
Met.DAPC.AC.plot<-ggplot(data = met.dapc.plot.data.AC, aes(x = LD1, y = LD2)) + 
  geom_point(data = met.dapc.plot.data.AC, size = 3,  aes(colour=Set, alpha = Set)) + 
  theme_classic()+
  scale_x_continuous(limits = c(-4, 4))+
  scale_y_continuous(limits = c(-4.7, 3.3))+
  scale_alpha_manual(values=c(rep(c(0.45, 0.8), 4)), name=NULL)+
  stat_conf_ellipse(data = met.dapc.plot.data.AC, geom = "polygon", 
                    aes(x = LD1, y = LD2, color = Set, fill=Set, alpha = Set))+
  scale_fill_manual(values =c(rep(AC.genotype.colors.o[1],2), rep(AC.genotype.colors.o[2],2), 
                              rep(AC.genotype.colors.o[3],2), rep(AC.genotype.colors.o[4],2)), name=NULL)+
  scale_colour_manual(values=c(rep(AC.genotype.colors.o[1],2), rep(AC.genotype.colors.o[2],2), 
                               rep(AC.genotype.colors.o[3],2), rep(AC.genotype.colors.o[4],2)), name=NULL)+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .075), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), 
        legend.box.background = element_rect(color = "black"))+
  xlab("LD 1 (63.83%)")+
  ylab("LD 2 (30.25%)")+
  guides(fill=guide_legend(nrow=2))

Met.DAPC.AC.plot


####DAPC Loci Contributions- A. cervicornis####

#Identifying loci within the top 90th percentile of loading scores
loadingplot(MetDAPC.AC$var.contr, axis=1, thres=as.numeric(quantile(MetDAPC.AC$var.contr[,1], .9)),  
            srt=90, adj=-0.1,  ylab="Discriminant Axis 1 Loadings", xlab="Loci", main=NULL)
AC.DF1.90<-c("AC_A23", "AC_A25", "AC_A32", "AC_A33", "AC_A39", "AC_A44", "AC_A47", "AC_A62", 
             "AC_A74", "AC_A80", "AC_A93", "AC_A103", "AC_A108", "AC_A109", "AC_A115", "AC_A119", 
             "AC_A120", "AC_A125", "AC_A126", "AC_A128", "AC_A134", "AC_A138", "AC_A176", 
             "AC_B18", "AC_D6", "AC_D23", "AC_D40")

loadingplot(MetDAPC.AC$var.contr, axis=2, thres=as.numeric(quantile(MetDAPC.AC$var.contr[,2], .9)),  
            srt=90, adj=-0.1,  ylab="Discriminant Axis 2 Loadings", xlab="Loci", main=NULL)
AC.DF2.90<-c("AC_A40", "AC_A79", "AC_A80", "AC_A103", "AC_A108", "AC_A115", "AC_A125", "AC_A127", 
             "AC_A134", "AC_A136", "AC_A138", "AC_A139", "AC_A164", "AC_A166", "AC_A170", 
             "AC_A177", "AC_A187", "AC_A189", "AC_D6", "AC_D12", "AC_D13", "AC_D22", "AC_D23")

##Prepare for panel plots
DAPC.AC.df1<- as_grob(function() {loadingplot(MetDAPC.AC$var.contr, axis=1, 
                                              thres=as.numeric(quantile(MetDAPC.AC$var.contr[,1], .9)),  
                                              srt=90, adj=-0.1,  ylab="Discriminant Axis 1 Loadings", xlab="Loci", 
                                              main=NULL)})

DAPC.AC.df2<- as_grob(function() {loadingplot(MetDAPC.AC$var.contr, axis=2, 
                                              thres=as.numeric(quantile(MetDAPC.AC$var.contr[,2], .9)),  
                                              srt=90, adj=-0.1,  ylab="Discriminant Axis 2 Loadings", xlab="Loci", 
                                              main=NULL)})

####Methylation DAPC- A. palmata####

##Create new dataframe to group by Significant Factors (Genotype and Season)
MetData.AP.set<-MetData.AP

##Create a new variable of Seasonal Period
#Cool: "Cooling" time points of Fall and Winter 
#Warm: "Warming" time points of Spring and Summer 
MetData.AP.set$Period<-"Cooling"
MetData.AP.set$Period[which(MetData.AP.set$Season=="Spring" | MetData.AP.set$Season=="Summer")]<-"Warming"
MetData.AP.set$Period<-factor(MetData.AP.set$Period, levels=c("Cooling", "Warming"), ordered=TRUE)

##Create a new variable of Set by Genotype and Period 
MetData.AP.set$Set<-paste(MetData.AP.set$Genotype, MetData.AP.set$Period, sep=".")
MetData.AP.set$Set<-factor(MetData.AP.set$Set, levels=c("AP1.Cooling", "AP1.Warming", 
                                                        "AP8.Cooling", "AP8.Warming", "AP10.Cooling", "AP10.Warming"), ordered=TRUE)

##Prepare Input Data
names(MetData.AP.set)
Gen_Data.AP <- df2genind(MetData.AP.set[,c(8:75)], ploidy=1, ncode=1)
pop.ord.AP<-factor(MetData.AP.set$Set, ordered=FALSE)

##Number of PCA's to retain as k-1
n.pca.AP<-length(levels(pop.ord.AP))-1

##DAPC Analysis with Parameters Obtained
MetDAPC.AP <- dapc(Gen_Data.AP, pop=pop.ord.AP, n.pca=n.pca.AP, n.da=2,  scale=FALSE,
                   var.contrib=TRUE)
summary(MetDAPC.AP)

##Calculate percent of variance explained by each axis
#LD1
(MetDAPC.AP$eig[1]/sum(MetDAPC.AP$eig))*100
#80.11%

#LD2
(MetDAPC.AP$eig[2]/sum(MetDAPC.AP$eig))*100
#12.92%

##Prepare for Plotting
met.dapc.plot.data.AP<-data.frame(MetDAPC.AP$ind.coord)
met.dapc.plot.data.AP$ID<-rownames(MetDAPC.AP$ind.coord)
met.dapc.plot.data.AP<-merge(met.dapc.plot.data.AP, MetData.AP.set)  


####Plot DAPC Scatter Plot of DNA Methylation- A. palmata####
Met.DAPC.AP.plot<-ggplot(data = met.dapc.plot.data.AP, aes(x = LD1, y = LD2)) + 
  geom_point(data = met.dapc.plot.data.AP, size = 3,  aes(colour=Set, alpha = Set)) + 
  theme_classic()+
  scale_x_continuous(limits = c(-5.2, 4.8))+
  scale_y_continuous(limits = c(-5, 5))+
  scale_alpha_manual(values=c(rep(c(0.45, 0.8), 3)), name=NULL)+
  stat_conf_ellipse(data = met.dapc.plot.data.AP, geom = "polygon", 
                    aes(x = LD1, y = LD2, color = Set, fill=Set, alpha = Set))+
  scale_fill_manual(values =c(rep(AP.genotype.colors.o[1],2), rep(AP.genotype.colors.o[2],2), 
                              rep(AP.genotype.colors.o[3],2)), name=NULL)+
  scale_colour_manual(values=c(rep(AP.genotype.colors.o[1],2), rep(AP.genotype.colors.o[2],2), 
                               rep(AP.genotype.colors.o[3],2)), name=NULL)+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, .075), legend.direction = "horizontal", legend.text=element_text(size=leg.txt.sz), 
        legend.box.background = element_rect(color = "black"))+
  xlab("LD 1 (80.11%)")+
  ylab("LD 2 (12.92%)")+
  guides(fill=guide_legend(nrow=2))

Met.DAPC.AP.plot


####DAPC Loci Contributions- A. palmata####

#Identifying loci within the top 90th percentile of loading scores
loadingplot(MetDAPC.AP$var.contr, axis=1, thres=as.numeric(quantile(MetDAPC.AP$var.contr[,1], .9)),  
            srt=90, adj=-0.1,  ylab="Discriminant Axis 1 Loadings", xlab="Loci", main=NULL)
AP.DF1.90<-c("AP_A24", "AP_A27", "AP_A33", "AP_A34", "AP_A36", "AP_A38", "AP_A39", "AP_A43", 
             "AP_A50", "AP_A51", "AP_A56", "AP_A60", "AP_A62", "AP_A91", "AP_B20", "AP_B35", 
             "AP_B37", "AC_B38", "AP_D23", "AP_D32")

loadingplot(MetDAPC.AP$var.contr, axis=2, thres=as.numeric(quantile(MetDAPC.AP$var.contr[,2], .9)),  
            srt=90, adj=-0.1,  ylab="Discriminant Axis 2 Loadings", xlab="Loci", main=NULL)
AP.DF2.90<-c("AP_A14", "AP_A35", "AP_A78", "AP_B20", "AP_B37", "AP_D8", "AP_D9", "AP_D15", 
             "AP_D17", "AP_D22", "AP_D32", "AP_D35", "AP_D44", "AP_D46", "AP_D51", "AP_D52", 
             "AP_D53")

##Prepare for panel plots
DAPC.AP.df1<- as_grob(function() {loadingplot(MetDAPC.AP$var.contr, axis=1, 
                                              thres=as.numeric(quantile(MetDAPC.AP$var.contr[,1], .9)),  
                                              srt=90, adj=-0.1,  ylab="Discriminant Axis 1 Loadings", xlab="Loci", 
                                              main=NULL)})

DAPC.AP.df2<- as_grob(function() {loadingplot(MetDAPC.AP$var.contr, axis=2, 
                                              thres=as.numeric(quantile(MetDAPC.AP$var.contr[,2], .9)),  
                                              srt=90, adj=-0.1,  ylab="Discriminant Axis 2 Loadings", xlab="Loci", 
                                              main=NULL)})


####***Figure 3, S9: Plot Methylation Panels- A. cervicornis and A. palmata ####

##Figure 3: Methylation Types and DAPC 
Met.AC.AP_fig<-plot_grid(Met.Type.AC.plot, Met.DAPC.AC.plot,
                         Met.Type.AP.plot,Met.DAPC.AP.plot,
                         rel_widths=c(0.75, 1, 0.75, 1), rel_heights=c(1, 1, 1, 1), 
                         ncol=2, nrow=2, byrow=T, labels = c('A', 'C', 'B', 'D'), 
                         label_size=panel.lab.sz, align="hV")
ggsave(filename="Fig3_Acer_and_Apal_Methylation.png", plot=Met.AC.AP_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="Fig3_Acer_and_Apal_Methylation.pdf", plot=Met.AC.AP_fig, dpi=300, width=12, height=12, units="in")

##Figure S9: Top Loci Loading Plots
Met.Loci_fig<-plot_grid(DAPC.AC.df1, DAPC.AC.df2,
                        DAPC.AP.df1, DAPC.AP.df2, 
                        rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1),
                        ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                        label_size=panel.lab.sz, align="V")

ggsave(filename="FigS9_Acer_and_Apal_Top_Loci.png", plot=Met.Loci_fig, dpi=300, width=16, height=11, units="in", bg="white")
ggsave(filename="FigS9_Acer_and_Apal_Top_Loci.pdf", plot=Met.Loci_fig, dpi=300, width=16, height=11, units="in", bg="white")


####**Methylation Part 5: Patterns of Influential Loci####

####Most Influential Loci from DAPC####

####A. cervicornis

##AC DF1 - Seasonal
length(AC.DF1.90) #27

##AC DF2 - Genotype
length(AC.DF2.90) #23

Top.Loci.AC.90<-unique(c(AC.DF1.90, AC.DF2.90))
length(Top.Loci.AC.90) #41
(length(Top.Loci.AC.90)/91)*100 #45.05%


####A. palmata

##AP DF1 - Genotype
length(AP.DF1.90) #20

##AP DF2 - Seasonal
length(AP.DF2.90) #17

Top.Loci.AP.90<-unique(c(AP.DF1.90, AP.DF2.90))
length(Top.Loci.AP.90) #34
(length(Top.Loci.AP.90)/68)*100 #50%


####Patterns of Most Influential Loci- A. cervicornis####

##Organize Metdata for heatmap annotation
Genotype.AC <- data.frame(MetData.AC[,2])
Site.AC <- data.frame(MetData.AC[,4])
Season.AC <- data.frame(MetData.AC[,5])
Metadata.AC <- cbind(Genotype.AC, Site.AC, Season.AC)
colnames(Metadata.AC)<-c("Genotype", "Site", "Season")
rownames(Metadata.AC)<- MetData.AC$ID
Metadata.AC$Period <- ifelse(Metadata.AC$Season %in% c("Fall", "Winter"), "Cooling", ifelse(Metadata.AC$Season %in% c("Spring", "Summer"), "Warming", NA))

##Set colors as names vectors for heatmap annotation 
AC.genotype.colors.o.ha <- c(AC8 = "#F5E926FF", AC10 = "#FDB130FF", AC12 = "#EB7556FF", AC15 = "#D7566CFF")
Period.colors.o.ha <- c(Cooling = "#455ED2FF", Warming = "#A51301FF")
Met.colors.o.ha <- c("4" = "#DEF5E5FF", "2" = "#5ACCADFF", "3" = "#3482A4FF", "1" = "#3C5397FF")

##Heatmap annotation
ha.AC <- HeatmapAnnotation(df = Metadata.AC[, c("Period", "Genotype")], which = "row", col = list(Period = Period.colors.o.ha, Genotype = AC.genotype.colors.o.ha))

##Organize Methylation State Data for heatmap
#Filtering for Loci in 90th percentile of loading scores on both DAPC DF1 and DF2 axes
patt.AC.90 <- subset(MetData.AC, select = names(MetData.AC) %in% Top.Loci.AC.90)

num_mat.AC.90 <- matrix(NA, nrow = nrow(patt.AC.90), ncol = length(patt.AC.90))
num_mat.AC.90[patt.AC.90 == "f"] <- 1
num_mat.AC.90[patt.AC.90 == "h"] <- 2
num_mat.AC.90[patt.AC.90 == "i"] <- 3
num_mat.AC.90[patt.AC.90 == "u"] <- 4

##Gower distance function to call in heatmap
gowerDist <- function(x){daisy(as.matrix(x), metric="gower")}


####Plot Methylation Heatmap- A. cervicornis####
Met.heatmap.AC.90.plot <- Heatmap(num_mat.AC.90, 
                                  clustering_distance_rows = gowerDist, 
                                  clustering_distance_columns = gowerDist, 
                                  column_title = "MSAP Loci", 
                                  column_title_gp = gpar(fontface = "bold"), 
                                  row_title = "Specimen", 
                                  row_title_gp = gpar(fontface = "bold"), 
                                  column_dend_height = unit(1, "cm"), 
                                  show_column_names = FALSE, 
                                  show_row_names = FALSE, col = Met.colors.o.ha, 
                                  border = TRUE, 
                                  row_km = 2, 
                                  row_km_repeats = 1000, 
                                  heatmap_legend_param = list(title = "Methylation", 
                                                              labels = c("NMT", "HMM", "ICM", "HPM")), 
                                  right_annotation = ha.AC)

Met.heatmap.AC.90.plot


####Patterns of Most Influential Loci- A. palmata####

##Organize Metdata for heatmap annotation
Genotype.AP <- data.frame(MetData.AP[,2])
Site.AP <- data.frame(MetData.AP[,4])
Season.AP <- data.frame(MetData.AP[,5])
Metadata.AP <- cbind(Genotype.AP, Site.AP, Season.AP)
colnames(Metadata.AP)<-c("Genotype", "Site", "Season")
rownames(Metadata.AP)<- MetData.AP$ID
Metadata.AP$Period <- ifelse(Metadata.AP$Season %in% c("Fall", "Winter"), "Cooling", ifelse(Metadata.AP$Season %in% c("Spring", "Summer"), "Warming", NA))

##Set colors as names vectors for heatmap annotation 
AP.genotype.colors.o.ha <- c(AP1 = "#B52F8CFF", AP8 = "#6E00A8FF", AP10 = "#240691FF")
Period.colors.o.ha <- c(Cooling = "#455ED2FF", Warming = "#A51301FF")
Met.colors.o.ha <- c("4" = "#DEF5E5FF", "2" = "#5ACCADFF", "3" = "#3482A4FF", "1" = "#3C5397FF")

##Heatmap annotation
ha.AP <- HeatmapAnnotation(df = Metadata.AP[, c("Period", "Genotype")], which = "row", col = list(Period = Period.colors.o.ha, Genotype = AP.genotype.colors.o.ha))

##Organize Methylation State Data for heatmap
#Filtering for Loci in 90th percentile of loading scores on both DAPC DF1 and DF2 axes
patt.AP.90 <- subset(MetData.AP, select = names(MetData.AP) %in% Top.Loci.AP.90)

num_mat.AP.90 <- matrix(NA, nrow = nrow(patt.AP.90), ncol = length(patt.AP.90))
num_mat.AP.90[patt.AP.90 == "f"] <- 1
num_mat.AP.90[patt.AP.90 == "h"] <- 2
num_mat.AP.90[patt.AP.90 == "i"] <- 3
num_mat.AP.90[patt.AP.90 == "u"] <- 4


Met.heatmap.AP.90.plot <- Heatmap(num_mat.AP.90,  clustering_distance_rows = gowerDist, 
                                  clustering_distance_columns = gowerDist, 
                                  column_title = "MSAP Loci", 
                                  column_title_gp = gpar(fontface = "bold"), 
                                  row_title = "Specimen",                       
                                  row_title_gp = gpar(fontface = "bold"), 
                                  column_dend_height = unit(1, "cm"),
                                  row_dend_width = unit(1, "cm"), 
                                  show_column_names = FALSE, 
                                  show_row_names = FALSE, col = Met.colors.o.ha, 
                                  border = TRUE, 
                                  row_km = 2, 
                                  row_km_repeats = 1000, 
                                  heatmap_legend_param = list(title = "Methylation", 
                                                              labels = c("NMT", "HMM", "ICM", "HPM")), 
                                  right_annotation = ha.AP)

Met.heatmap.AP.90.plot


####***Figure S10: Plot Methylation Heatmaps- A. cervicornis and A. palmata ####  

Met.heatmap.AC.90.grid<-grid.grabExpr(Met.heatmap.AC.90.plot)

Met.heatmap.AP.90.grid<-grid.grabExpr(Met.heatmap.AP.90.plot) 

Met.Heatmap_fig<-plot_grid(Met.heatmap.AC.90.grid, Met.heatmap.AP.90.grid,
                           rel_widths=c(1, 1), rel_heights=c(1, 1),
                           ncol=2, nrow=1, byrow=T, labels = c('A', 'B'), 
                           label_size=panel.lab.sz, align="h")

ggsave(filename="FigS10_Acer_and_Apal_Methylation_Heatmap.png", plot=Met.Heatmap_fig, dpi=300, width=12, height=8, units="in")
ggsave(filename="FigS10_Acer_and_Apal_Methylation_Heatmap.pdf", plot=Met.Heatmap_fig, dpi=300, width=12, height=8, units="in")



#-----*Physiology and Epigenetic Plasticity Comparisons------------------------------------------------------------


####**Plasticity Part 1: RDA Relating DNA Methylation to Physiological Metrics####


####Distance Based RDA of Methylation and Physiology- Acropora cervicornis####

##Filtering Methylation Data for Loci in 90th percentile of loading scores on both DAPC DF1 and DF2 axes
names(MetData.AC.set)
AC_Data<-cbind(MetData.AC.set[,c(1:7,99:100)], MetData.AC.set[,names(MetData.AC.set) %in% Top.Loci.AC.90])

##Merge Methylation and Physiology Data
names(CoralData.log.full.AC)
AC_Data<-merge(AC_Data, CoralData.log.full.AC, all.x=FALSE, all.y=FALSE)
rownames(AC_Data)<-AC_Data$ID

##Distance-based redundancy analysis of Methylation patterns explained by Physiological metrics
names(AC_Data)
AC.met.phys.rda<-dbrda(daisy(AC_Data[, c(10:50)], "gower") ~ TP_ug.cm2_C + TP_ug.cm2_S + 
                         Chl.a_ug.cm2 + Chl.c2_ug.cm2 + AFDW_mg.cm2_C + AFDW_mg.cm2_S, 
                       data=AC_Data, dist="gower")

##Check model significance
anova(AC.met.phys.rda)
#           Df SumOfSqs      F Pr(>F)    
# Model      6    5.328 5.6595  0.001 ***
# Residual 224   35.148 
#The model is significant 

##Check variance explained by model
(summary(AC.met.phys.rda)$constr.chi/summary(AC.met.phys.rda)$tot.chi)*100
#13.16% of the variance in Methylation patterns is constrained by Physiological metrics 

RsquareAdj(AC.met.phys.rda)$adj.r.squared*100
#Adjusted: The model explains 10.84% of the variation in coral DNA methylation in A cerv

##Check variance explained by each axis
anova(AC.met.phys.rda, by="axis")
#           Df SumOfSqs       F Pr(>F)    
# dbRDA1     1    4.120 26.2548  0.001 ***
# dbRDA2     1    0.571  3.6366  0.025 *  
# dbRDA3     1    0.318  2.0240  0.383    
# dbRDA4     1    0.236  1.5045  0.583    
# dbRDA5     1    0.044  0.2781  1.000    
# dbRDA6     1    0.041  0.2593  0.979    
# Residual 224   35.148 

summary(AC.met.phys.rda)$cont
#dbRDA1 explains 10.18% of total variance 
#dbRDA2 explains 1.41% of total variance 

summary(AC.met.phys.rda)$concont
#dbRDA1 explains 77.32% of constrained variance 
#dbRDA2 explains 10.71% of constrained variance 

##Check variance explained by each Physiological metric
anova.cca(AC.met.phys.rda, by="terms")
#                 Df SumOfSqs       F Pr(>F)    
# TP_ug.cm2_C     1    1.029  6.5592  0.001 ***
# TP_ug.cm2_S     1    2.955 18.8312  0.001 ***
# Chl.a_ug.cm2    1    0.347  2.2090  0.020 *  
# Chl.c2_ug.cm2   1    0.387  2.4653  0.018 *  
# AFDW_mg.cm2_C   1    0.101  0.6460  0.724    
# AFDW_mg.cm2_S   1    0.509  3.2466  0.002 ** 
# Residual      224   35.148     

##Prepare for Plotting
AC.met.phys.rda.scores<- data.frame(scores(AC.met.phys.rda)$sites)
AC.met.phys.rda.scores$ID<- rownames(AC.met.phys.rda.scores)
names(AC_Data)
AC.met.phys.rda.scores<-merge(AC_Data[,c(1:9)], AC.met.phys.rda.scores)
names(AC.met.phys.rda.scores)

AC.met.phys.rda.phys<-data.frame(summary(AC.met.phys.rda)$biplot[,1:2]) * ordiArrowMul(AC.met.phys.rda, display="bp")
AC.met.phys.rda.phys$Covariate<-rownames(AC.met.phys.rda.phys)
AC.met.phys.rda.phys$Covariate
AC.met.phys.rda.phys$Metric<-c("Protein_Host", "Protein_Symb", "Chlorophyll-a", "Chlorophyll-c2", "Biomass_Host", "Biomass_Symb")
AC.met.phys.rda.phys

####Plot dbRDA Biplot of DNA Methylation vs Physiological Metrics- A. cervicornis####
Met.Phys.RDA.AC.plot<-ggplot(data = AC.met.phys.rda.scores, aes(x = dbRDA1, y = dbRDA2)) + 
  geom_point(data = AC.met.phys.rda.scores, size = 3, alpha = 0.8, aes(colour = Period, fill = Period)) + 
  theme_classic()+
  xlim(-5.25, 5.25)+
  ylim(-5.25, 5.25)+
  scale_colour_manual(values=Period.colors.o)+  
  scale_fill_manual(values=Period.colors.o)+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, 0.05), legend.direction = "horizontal", legend.text=element_text(size=leg.title.sz), 
        legend.title=element_text(size=leg.title.sz), legend.box.background = element_rect(color = "black"))+
  xlab("RDA 1 (77.32% Fitted, 10.18% Total)")+
  ylab("RDA 2 (10.71% Fitted, 1.41% Total)")+
  geom_segment(aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               data = AC.met.phys.rda.phys, size =1, alpha = 0.5, colour = "grey20", 
               arrow=arrow(length=unit(0.01,"npc")))+
  geom_text(data = AC.met.phys.rda.phys[-5,], 
            aes(x=c(dbRDA1[1], dbRDA1[2]-1.5, dbRDA1[3]-0.25, dbRDA1[4]+0.75, dbRDA1[5]+0.15), 
                y = c(dbRDA2[1],dbRDA2[2],dbRDA2[3]-0.15,dbRDA2[4]+0.15, dbRDA2[5]+0.5),
                label = Metric, fontface = "bold",
                hjust=0.5*(1-sign(dbRDA1)),vjust=0.6*(1-sign(dbRDA2))), 
            colour = "grey20", size=levels.sz-1)+
  geom_text(data = AC.met.phys.rda.phys[5,], aes(x=c(dbRDA1[1]+0.2), y = c(dbRDA2[1]+0.25), label = Metric, fontface = "italic", 
                                                 hjust=0.5*(1-sign(dbRDA1)),vjust=0.6*(1-sign(dbRDA2))), 
            colour = "grey20", size=levels.sz-1)+
  annotate("text", x = 0.25, y = -4.25, label=expression(bolditalic(paste("Adj. ",R^2, " = 0.108, p = 0.001"))), size=sig.sz, hjust = 0)

Met.Phys.RDA.AC.plot


####Distance Based RDA of Methylation and Physiology- Acropora palmata####

##Filtering Methylation Data for Loci in 90th percentile of loading scores on both DAPC DF1 and DF2 axes
names(MetData.AP.set)
AP_Data<-cbind(MetData.AP.set[,c(1:7,76:77)], MetData.AP.set[,names(MetData.AP.set) %in% Top.Loci.AP.90])

##Merge Methylation and Physiology Data
names(CoralData.log.full.AP)
AP_Data<-merge(AP_Data, CoralData.log.full.AP, all.x=FALSE, all.y=FALSE)
rownames(AP_Data)<-AP_Data$ID

##Distance-based redundancy analysis of Methylation patterns explained by Physiological metrics
names(AP_Data)
AP.met.phys.rda<-dbrda(daisy(AP_Data[, c(10:42)], "gower") ~ TP_ug.cm2_C + TP_ug.cm2_S + 
                         Chl.a_ug.cm2 + Chl.c2_ug.cm2 + AFDW_mg.cm2_C + AFDW_mg.cm2_S, 
                       data=AP_Data, dist="gower")



##Check model significance
anova(AP.met.phys.rda)
#           Df SumOfSqs      F Pr(>F)    
# Model      6   4.1423 4.6794  0.001 ***
# Residual 100  14.7537 
#The model is significant 

##Check variance explained by model
(summary(AP.met.phys.rda)$constr.chi/summary(AP.met.phys.rda)$tot.chi)*100
#21.92% of the variance in Methylation patterns is constrained by Physiological metrics 

RsquareAdj(AP.met.phys.rda)$adj.r.squared*100
#Adjusted: The model explains 17.24% of the variation in coral DNA methylation in A pal

##Check variance explained by each axis
anova(AP.met.phys.rda, by="axis")
#           Df SumOfSqs       F Pr(>F)    
# dbRDA1     1   2.5596 17.3491  0.001 ***
# dbRDA2     1   0.7716  5.2299  0.008 ** 
# dbRDA3     1   0.3911  2.6506  0.222    
# dbRDA4     1   0.3149  2.1342  0.289    
# dbRDA5     1   0.0724  0.4907  0.987    
# dbRDA6     1   0.0328  0.2220  0.950    
# Residual 100  14.7537   

summary(AP.met.phys.rda)$cont
#dbRDA1 explains 13.55% of total variance 
#dbRDA2 explains 4.08% of total variance 

summary(AP.met.phys.rda)$concont
#dbRDA1 explains 61.79% of constrained variance 
#dbRDA2 explains 18.63% of constrained variance 

##Check variance explained by each Physiological metric
anova(AP.met.phys.rda, by="terms")
#                 Df SumOfSqs       F Pr(>F)    
# TP_ug.cm2_C     1   1.5757 10.6798  0.001 ***
# TP_ug.cm2_S     1   1.1414  7.7362  0.001 ***
# Chl.a_ug.cm2    1   0.3857  2.6142  0.018 *  
# Chl.c2_ug.cm2   1   0.3923  2.6587  0.013 *  
# AFDW_mg.cm2_C   1   0.4425  2.9990  0.003 ** 
# AFDW_mg.cm2_S   1   0.2049  1.3886  0.221    
# Residual      100  14.7537       


##Prepare for Plotting
AP.met.phys.rda.scores<- data.frame(scores(AP.met.phys.rda)$sites)
AP.met.phys.rda.scores$ID<- rownames(AP.met.phys.rda.scores)
names(AP_Data)
AP.met.phys.rda.scores<-merge(AP_Data[,c(1:9)], AP.met.phys.rda.scores)
names(AP.met.phys.rda.scores)

AP.met.phys.rda.phys<-data.frame(summary(AP.met.phys.rda)$biplot[,1:2]) * ordiArrowMul(AP.met.phys.rda, display="bp")
AP.met.phys.rda.phys$Covariate<-rownames(AP.met.phys.rda.phys)
AP.met.phys.rda.phys$Covariate
AP.met.phys.rda.phys$Metric<-c("Protein_Host", "Protein_Symb", "Chlorophyll-a", "Chlorophyll-c2", "Biomass_Host", "Biomass_Symb")
AP.met.phys.rda.phys


####Plot dbRDA Biplot of DNA Methylation vs Physiological Metrics- A. palmata####
Met.Phys.RDA.AP.plot<-ggplot(data = AP.met.phys.rda.scores, aes(x = dbRDA1, y = dbRDA2)) + 
  geom_point(data = AP.met.phys.rda.scores, size = 3, alpha = 0.8, aes(colour = Period, fill = Period)) + 
  theme_classic()+
  xlim(-3, 4)+
  ylim(-3.5, 3.5)+
  scale_colour_manual(values=Period.colors.o)+  
  scale_fill_manual(values=Period.colors.o)+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.position=c(.5, 0.05), legend.direction = "horizontal", legend.text=element_text(size=leg.title.sz), 
        legend.title=element_text(size=leg.title.sz), legend.box.background = element_rect(color = "black"))+
  xlab("RDA 1 (61.79% Fitted, 13.55% Total)")+
  ylab("RDA 2 (18.63% Fitted, 4.08% Total)")+
  geom_segment(aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               data = AP.met.phys.rda.phys, size =1, alpha = 0.5, colour = "grey20", 
               arrow=arrow(length=unit(0.01,"npc")))+
  geom_text(data = AP.met.phys.rda.phys[-6,], 
            aes(x=c(dbRDA1[1]-0.5, dbRDA1[2]+1, dbRDA1[3]-0.6, dbRDA1[4]-0.75, dbRDA1[5]), 
                y = c(dbRDA2[1]+0.15,dbRDA2[2]+0.15, dbRDA2[3]+0.2, dbRDA2[4]-0.3, dbRDA2[5]),
                label = Metric, fontface = "bold",
                hjust=0.5*(1-sign(dbRDA1)),vjust=0.6*(1-sign(dbRDA2))), 
            colour = "grey20", size=levels.sz-1)+
  geom_text(data = AP.met.phys.rda.phys[6,], aes(x=c(dbRDA1-1.75), y = c(dbRDA2+0.2), label = Metric, fontface = "italic", 
                                                 hjust=0.5*(1-sign(dbRDA1)),vjust=0.6*(1-sign(dbRDA2))), 
            colour = "grey20", size=levels.sz-1)+
  annotate("text", x = 0.5, y = -2.8, label=expression(bolditalic(paste("Adj. ",R^2, " = 0.172, p = 0.001"))), size=sig.sz, hjust = 0)

Met.Phys.RDA.AP.plot


####**Plasticity Part 2: Plasticity as Shifts in Centroid Distances####


####Quantify Physiological Plasticity using Centroid Distances- A. cervicornis####

##Obtain centroid coordinates of each Sample Set (Replicates of the same Genotype, Site, Season) with betadisper type="centroid"
Phys.disp.AC<-betadisper(vegdist(CoralData.log.full.AC[,10:15], "euclidean"), CoralData.log.full.AC$SampSet, type="centroid")

##Calculate the variance of Physiology explained by the first two PCoA axes
Phys.disp.AC$eig[1]/sum(Phys.disp.AC$eig[which(Phys.disp.AC$eig>0)])*100
Phys.disp.AC$eig[2]/sum(Phys.disp.AC$eig[which(Phys.disp.AC$eig>0)])*100
sum(Phys.disp.AC$eig[1:2])/sum(Phys.disp.AC$eig[which(Phys.disp.AC$eig>0)])*100
#82.13% explained by first two axes 

##Calculate the Euclidean distance between centroids using coordinates of the first 2 PCoA axes
Phys.centdist.AC<-as.matrix(dist(Phys.disp.AC$centroids[,c(1:2)], method = "euclidean"))

##Create a dataframe of distances between centroids of Sample Sets
Phys.cents.AC<-data.frame(SetA=colnames(Phys.centdist.AC)[col(Phys.centdist.AC)], SetB=rownames(Phys.centdist.AC)[row(Phys.centdist.AC)], Cent.Dist.Phys=c(Phys.centdist.AC))

##Remove rows comparing self-pairs of Sample Sets 
Phys.cents.AC<-Phys.cents.AC[-c(which(Phys.cents.AC$SetA==Phys.cents.AC$SetB)),]

##Remove repeats (Ex: Set 1 vs Set 2 and Set 2 vs Set 1)
Phys.cents.AC<-Phys.cents.AC[!duplicated(Phys.cents.AC$Cent.Dist.Phys),]

##Split Sets column to add variables of interest
Phys.cents.AC<-separate(Phys.cents.AC, "SetA", into=c("Site.A", "Genotype.A", "Season.A"), remove=FALSE)
Phys.cents.AC<-separate(Phys.cents.AC, "SetB", into=c("Site.B", "Genotype.B", "Season.B"), remove=FALSE)

##Filter to only keep rows comparing the same Genotype
Phys.cents.AC<-Phys.cents.AC[c(which(Phys.cents.AC$Genotype.A==Phys.cents.AC$Genotype.B)),]

##Add variable of Pair for merging
Phys.cents.AC$Pair<-paste(Phys.cents.AC$SetA, Phys.cents.AC$SetB, sep="_")

##Dataframe with distances between Seasons by only keeping rows comparing the same Site
Phys.cents.AC.Seas<-Phys.cents.AC[c(which(Phys.cents.AC$Site.A==Phys.cents.AC$Site.B)),]


####Quantify Physiological Plasticity using Centroid Distances- A. palmata####

##Obtain centroid coordinates of each Sample Set (Replicates of the same Genotype, Site, Season) with betadisper type="centroid"
Phys.disp.AP<-betadisper(vegdist(CoralData.log.full.AP[,10:15], "euclidean"), CoralData.log.full.AP$SampSet, type="centroid")

##Calculate the variance of Physiology explained by the first two PCoA axes
Phys.disp.AP$eig[1]/sum(Phys.disp.AP$eig[which(Phys.disp.AP$eig>0)])*100
Phys.disp.AP$eig[2]/sum(Phys.disp.AP$eig[which(Phys.disp.AP$eig>0)])*100
sum(Phys.disp.AP$eig[1:2])/sum(Phys.disp.AP$eig[which(Phys.disp.AP$eig>0)])*100
#84.61% explained by first two axes 

##Calculate the Euclidean distance between centroids using coordinates of the first 2 PCoA axes
Phys.centdist.AP<-as.matrix(dist(Phys.disp.AP$centroids[,c(1:2)], method = "euclidean"))

##Create a dataframe of distances between centroids of Sample Sets
Phys.cents.AP<-data.frame(SetA=colnames(Phys.centdist.AP)[col(Phys.centdist.AP)], SetB=rownames(Phys.centdist.AP)[row(Phys.centdist.AP)], Cent.Dist.Phys=c(Phys.centdist.AP))

##Remove rows comparing self-pairs of Sample Sets 
Phys.cents.AP<-Phys.cents.AP[-c(which(Phys.cents.AP$SetA==Phys.cents.AP$SetB)),]

##Remove repeats (Ex: Set 1 vs Set 2 and Set 2 vs Set 1)
Phys.cents.AP<-Phys.cents.AP[!duplicated(Phys.cents.AP$Cent.Dist.Phys),]

##Split Sets column to add variables of interest
Phys.cents.AP<-separate(Phys.cents.AP, "SetA", into=c("Site.A", "Genotype.A", "Season.A"), remove=FALSE)
Phys.cents.AP<-separate(Phys.cents.AP, "SetB", into=c("Site.B", "Genotype.B", "Season.B"), remove=FALSE)

##Filter to only keep rows comparing the same Genotype
Phys.cents.AP<-Phys.cents.AP[c(which(Phys.cents.AP$Genotype.A==Phys.cents.AP$Genotype.B)),]

##Add variable of Pair for merging
Phys.cents.AP$Pair<-paste(Phys.cents.AP$SetA, Phys.cents.AP$SetB, sep="_")

##Dataframe with distances between Seasons by only keeping rows comparing the same Site
Phys.cents.AP.Seas<-Phys.cents.AP[c(which(Phys.cents.AP$Site.A==Phys.cents.AP$Site.B)),]


####Quantify Methylation Plasticity using Centroid Distances- A. cervicornis####

##Obtain centroid coordinates of each Sample Set (Replicates of the same Genotype, Site, Season) with betadisper type="centroid"
#Filtering for Most Influential Loci (90th percentile of loading scores on DAPC axes)
Met.disp.AC<-betadisper(daisy(MetData.AC[,names(MetData.AC) %in% Top.Loci.AC.90], "gower"), MetData.AC$SampSet, type="centroid")

##Calculate the variance of Methylation explained by the first two PCoA axes
Met.disp.AC$eig[1]/sum(Met.disp.AC$eig[which(Met.disp.AC$eig>0)])*100
Met.disp.AC$eig[2]/sum(Met.disp.AC$eig[which(Met.disp.AC$eig>0)])*100
sum(Met.disp.AC$eig[1:2])/sum(Met.disp.AC$eig[which(Met.disp.AC$eig>0)])*100
#27.12% explained by first two axes 

##Calculate the Euclidean distance between centroids using coordinates of the first 2 PCoA axes
Met.centdist.AC<-as.matrix(dist(Met.disp.AC$centroids[,c(1:2)], method = "euclidean"))

##Create a dataframe of distances between centroids of Sample Sets
Met.cents.AC<-data.frame(SetA=colnames(Met.centdist.AC)[col(Met.centdist.AC)], 
                         SetB=rownames(Met.centdist.AC)[row(Met.centdist.AC)], 
                         Cent.Dist.Met=c(Met.centdist.AC))

##Remove rows comparing self-pairs of Sample Sets 
Met.cents.AC<-Met.cents.AC[-c(which(Met.cents.AC$SetA==Met.cents.AC$SetB)),]

##Remove repeats (Ex: Set 1 vs Set 2 and Set 2 vs Set 1)
Met.cents.AC<-Met.cents.AC[!duplicated(Met.cents.AC$Cent.Dist.Met),]

##Split Sets column to add variables of interest
Met.cents.AC<-separate(Met.cents.AC, "SetA", into=c("Site.A", "Genotype.A", "Season.A"), remove=FALSE)
Met.cents.AC<-separate(Met.cents.AC, "SetB", into=c("Site.B", "Genotype.B", "Season.B"), remove=FALSE)

##Filter to only keep rows comparing the same Genotype
Met.cents.AC<-Met.cents.AC[c(which(Met.cents.AC$Genotype.A==Met.cents.AC$Genotype.B)),]

##Add variable of Pair for merging
Met.cents.AC$Pair<-paste(Met.cents.AC$SetA, Met.cents.AC$SetB, sep="_")

##Dataframe with distances between Seasons by only keeping rows comparing the same Site
Met.cents.AC.Seas<-Met.cents.AC[c(which(Met.cents.AC$Site.A==Met.cents.AC$Site.B)),]


####Quantify Methylation Plasticity using Centroid Distances- A. palmata####

##Obtain centroid coordinates of each Sample Set (Replicates of the same Genotype, Site, Season) with betadisper type="centroid"
#Filtering for Most Influential Loci (90th percentile of loading scores on DAPC axes)
Met.disp.AP<-betadisper(daisy(MetData.AP[,names(MetData.AP) %in% Top.Loci.AP.90], "gower"), MetData.AP$SampSet, type="centroid")

##Calculate the variance of Methylation explained by the first two PCoA axes
Met.disp.AP$eig[1]/sum(Met.disp.AP$eig[which(Met.disp.AP$eig>0)])*100
Met.disp.AP$eig[2]/sum(Met.disp.AP$eig[which(Met.disp.AP$eig>0)])*100
sum(Met.disp.AP$eig[1:2])/sum(Met.disp.AP$eig[which(Met.disp.AP$eig>0)])*100
#34.95% explained by first two axes 

##Calculate the Euclidean distance between centroids using coordinates of the first 2 PCoA axes
Met.centdist.AP<-as.matrix(dist(Met.disp.AP$centroids[,c(1:2)], method = "euclidean"))

##Create a dataframe of distances between centroids of Sample Sets
Met.cents.AP<-data.frame(SetA=colnames(Met.centdist.AP)[col(Met.centdist.AP)], SetB=rownames(Met.centdist.AP)[row(Met.centdist.AP)], Cent.Dist.Met=c(Met.centdist.AP))

##Remove rows comparing self-pairs of Sample Sets 
Met.cents.AP<-Met.cents.AP[-c(which(Met.cents.AP$SetA==Met.cents.AP$SetB)),]

##Remove repeats (Ex: Set 1 vs Set 2 and Set 2 vs Set 1)
Met.cents.AP<-Met.cents.AP[!duplicated(Met.cents.AP$Cent.Dist.Met),]

##Split Sets column to add variables of interest
Met.cents.AP<-separate(Met.cents.AP, "SetA", into=c("Site.A", "Genotype.A", "Season.A"), remove=FALSE)
Met.cents.AP<-separate(Met.cents.AP, "SetB", into=c("Site.B", "Genotype.B", "Season.B"), remove=FALSE)

##Filter to only keep rows comparing the same Genotype
Met.cents.AP<-Met.cents.AP[c(which(Met.cents.AP$Genotype.A==Met.cents.AP$Genotype.B)),]

##Add variable of Pair for merging
Met.cents.AP$Pair<-paste(Met.cents.AP$SetA, Met.cents.AP$SetB, sep="_")

##Dataframe with distances between Seasons by only keeping rows comparing the same Site
Met.cents.AP.Seas<-Met.cents.AP[c(which(Met.cents.AP$Site.A==Met.cents.AP$Site.B)),]



####**Plasticity Part 3: Phenotypic vs Epigenetic Seasonal Plasticity####


####Phenotypic vs Epigenetic Seasonal Plasticity- A. cervicornis####

##Merge Physiology and Methylation Seasonal plasticity dataframes
Cent.dists.AC.Seas<-merge(Phys.cents.AC.Seas, Met.cents.AC.Seas)

##Add Genotype as ordered factor
Cent.dists.AC.Seas$Genotype<-Cent.dists.AC.Seas$Genotype.A
Cent.dists.AC.Seas$Genotype<-factor(Cent.dists.AC.Seas$Genotype, levels=c("AC8", "AC10", "AC12", "AC15"), ordered=TRUE)


##Correlation between Shifts in Physiology and Methylation across Seasons
cor.test(Cent.dists.AC.Seas$Cent.Dist.Phys, Cent.dists.AC.Seas$Cent.Dist.Met, method="spearman")
# S = 34648, p-value = 0.0001149
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.4429224 

####Plot Methylation vs Physiology Seasonal Plasticity- A. cervicornis####
Phys.Met.AC.plot<- ggplot(data=Cent.dists.AC.Seas, aes(x=Cent.Dist.Met,y=Cent.Dist.Phys)) +
  geom_point(aes(colour = Genotype), alpha=0.8,  size=3) +
  scale_fill_manual(values =AC.genotype.colors.o)+
  scale_colour_manual(values =AC.genotype.colors.o, name="Genotype")+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz))+
  labs(x= "DNA Methylation Centroid Distance", y="Physiology Centroid Distance")+
  geom_smooth(method=lm, color=Species.colors.o[1], aes(fill="#FAEBDDFF"), alpha=0.2, show.legend = FALSE)+
  annotate("text", x = 0.25, y = 0.05, label=expression(bolditalic(paste(r[S], " = 0.443, p < 0.001"))), size=sig.sz, hjust = 0)

Phys.Met.AC.plot

####Phenotypic vs Epigenetic Seasonal Plasticity- A. palmata####

##Merge Physiology and Methylation Seasonal plasticity dataframes
Cent.dists.AP.Seas<-merge(Phys.cents.AP.Seas, Met.cents.AP.Seas)

##Add Genotype as ordered factor
Cent.dists.AP.Seas$Genotype<-Cent.dists.AP.Seas$Genotype.A
Cent.dists.AP.Seas$Genotype<-factor(Cent.dists.AP.Seas$Genotype, levels=c("AP1", "AP8", "AP10"), ordered=TRUE)

##Correlation between Shifts in Physiology and Methylation across Seasons
cor.test(Cent.dists.AP.Seas$Cent.Dist.Phys, Cent.dists.AP.Seas$Cent.Dist.Met, method="spearman")
# S = 11178, p-value = 8.507e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.573928 


####Plot Methylation vs Physiology Seasonal Plasticity- A. palmata####
Phys.Met.AP.plot<-ggplot(data=Cent.dists.AP.Seas, aes(x=Cent.Dist.Met,y=Cent.Dist.Phys)) +
  geom_point(aes(colour = Genotype), alpha=0.8,  size=3) +
  scale_fill_manual(values =AP.genotype.colors.o)+
  scale_colour_manual(values =AP.genotype.colors.o, name="Genotype")+
  theme_classic()+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz))+
  labs(x= "DNA Methylation Centroid Distance", y="Physiology Centroid Distance")+
  geom_smooth(method=lm, color=Species.colors.o[2], aes(fill="#FAEBDDFF"), alpha=0.2, show.legend = FALSE)+
  annotate("text", x = 0.25, y = 0.05, label=expression(bolditalic(paste(r[S], " =  0.574, p < 0.001"))), size=sig.sz, hjust = 0)

Phys.Met.AP.plot

####***Figure 4: Plot Physiology vs Methylation Panels- A. cervicornis and A. palmata ####
Phys.Met_fig<-plot_grid(Met.Phys.RDA.AC.plot, Met.Phys.RDA.AP.plot, 
                        Phys.Met.AC.plot, Phys.Met.AP.plot,
                        rel_widths=c(1, 1, 1, 1), rel_heights=c(1, 1, 1, 1), 
                        ncol=2, nrow=2, byrow=T, labels = c('A', 'B', 'C', 'D'), 
                        label_size=20, align="h")
ggsave(filename="Fig4_Acer_and_Apal_Physiology_Methylation.png", plot=Phys.Met_fig, dpi=300, width=12, height=12, units="in")
ggsave(filename="Fig4_Acer_and_Apal_Physiology_Methylation.pdf", plot=Phys.Met_fig, dpi=300, width=12, height=12, units="in")




