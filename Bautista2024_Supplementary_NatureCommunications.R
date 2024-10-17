##Supplementary Figures Bautista_2024###
#In order to place the figures in the same folder as the data tables,
#you need to specify the name of this directory in the 'Set Directory' section below.

####Packages####
rm(list=ls())
#install.packages("agricolae")
library(agricolae)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(BiocManager)
#install.packages("broom")
library(broom)
#install.packages("car",type="binary")
library(car)
#install.packages("cowplot")
library(cowplot)
#install.packages("data.table")
library(data.table)
#install.packages("dplyr")
library(dplyr)
#install.packages("drc")
library(drc)
#install.packages("egg")
library(egg)
#install.packages("factoextra")
library(factoextra)
#install.packages("foreign")
library(foreign)
#install.packages("flowCore")
library(flowCore)
#install.packages("flowViz")
library(flowViz)
#install.packages("gdata")
library(gdata)
#install.packages("ggimage")
library(ggimage)
#install.packages("ggpattern")
library(ggpattern)
#install.packages("ggsignif")
library(ggsignif)
#install.packages("ggtext")
library(ggtext)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggprism")
library(ggprism)
#install.packages("grid")
library(ggpubr)
#install.packages("ggthemes")
library(ggthemes)
#BiocManager::install("ggtree")
library(ggtree)
#install.packages("glue")
library(glue)
#install.packages("grid")
library(grid)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("growthcurver") 
library(growthcurver)
#install.packages("gtable")
library(gtable)
#install.packages("janitor")
library(janitor)
#install.packages("lme4",type="binary")
library(lme4)
#install.packages("magick")
library(magick)
#install.packages("magrittr")
library(magrittr)
#install.packages("Matrix") 
library(Matrix)
#install.packages("modelr")
library(modelr)
#install.packages("multcomp")
library(multcomp)
#install.packages("multcompView")
library(multcompView)
#install.packages("nlme") 
library(nlme)
#install.packages("pdftools")
library(pdftools)
#install.packages("pegas")
library(pegas)
#install.packages("plyr")  
library(plyr)
#install.packages("png")
library(png)
#install.packages("quantreg", type="binary")
library(quantreg)
#install.packages("readr")
library(readr)
#install.packages("readxl")
library(readxl)
#install.packages("reshape2")
library(reshape2)
#install.packages("rstatix")
library(rstatix)
#install.packages("stats")
library(stats)
#install.packages("stringr")
library(stringr)
#install.packages("svglite")
library(svglite)
#BiocManager::install("trackViewer")
library(trackViewer)
#install.packages("tidyr")
library(tidyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("UpSetR")
library(UpSetR)
#BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
#install.packages("vcfR")
library(vcfR)
#################

####Set directory####
#Define the directory where you save all the data from Bautista_2024
setwd("")
#################

####Get Data####
Supp_coverage <- read_csv("1Supp_coverage.csv")
Supp_genomics <- read_csv("2Supp_growth.csv")
Supp_ploidy <- read.csv("4Supp_ploidy.csv",stringsAsFactors = F)
Supp_ploidy_heatmap <- read.csv("4Supp_ploidy_heatmap.csv")
Supp_Table_counts_aneuploidies <- read.csv("5Supp_Table_counts_aneuploidies.csv")
Supp_Aneuploidy_chrom <- read_excel("6Supp_Aneuploidy_chrom.xlsx")
Supp_genome <- read_csv("7Supp_new_genome3c_new.csv")
Supp_SNPs<- read_csv("9Supp_SNP.csv")
Supp_LOH_chrom<-read_csv("8Supp_Table_LOH.csv")
Supp_counts_LOH <- read_csv("8Supp_Table_counts_LOH.csv")
Supp10_GO <- image_read_pdf("10Supp_GO.pdf")
Supp_PDR1 <- readPNG("11Supp_PDR1.png")
Supp_candida <-image_read_pdf("12Supp_candida.pdf")
Supp_raw_cov<-read.csv("13Supp_raw_cov.csv")
Supp_growth<-read_csv("14Supp_growth.csv")
Supp_expression<-read_csv("14Supp_expression.csv")
Supp_growth_plasmids<-read.csv("15Supp_growth_plasmids.csv")
Supp_growth_curves_plasmids<-read.csv("15Supp_growth_curves_plasmids.csv")
Supp_growth_plasmids2<-read.csv("15Supp_growth_plasmids2.csv")
Supp_growth_curves_plasmids2<-read.csv("15Supp_growth_curves_plasmids2.csv")
Supp16_boxplot <- read_csv("16Supp_boxplot.csv")
Supp16_curves <- read_csv("16Supp_curves.csv")
Supp_expevol<-read_csv("17Supp_expevol.csv")
Supp17_Sanger <-image_read_pdf("17Supp_Sanger.pdf")
Supp17_LOH_table_parents <- read_excel("17Supp_LOH_parents_table.xlsx")
Supp17_LOH_parents_table_proportion<-read_excel("17Supp_LOH_parents_table_proportion.xlsx")
Supp17_LOH_parents_table_proportion_hybrid<-read_csv("17Supp_LOH_parents_table_proportion_hybrid.csv")
Supp17_LOH_line13_final<-read_excel("17Supp_LOH_line13_final.xlsx")
Supp17_LOH_line28_final<- read_excel("17Supp_LOH_line28_final.xlsx")
Supp17_LOH_line30_final<-read_excel("17Supp_LOH_line30_final.xlsx")
Supp_loh_table_day2<-read_csv("18Supp_loh_table_day2.csv")
Supp_fdata<-read_csv("18Supp_fdata.csv")
##########

############Figure Supplementary 1############
#Define some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Supp_coverage$library <- "Sample"

hybrid <-Supp_coverage %>% filter(Strain=="3Hybrid")
hybrid_noNQO <- hybrid %>% filter(Type!="Evolved_NQO")
hybrid_NQO <- hybrid %>% filter(Type=="Evolved_NQO")

hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=10)
hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=25)

nohybrid <-Supp_coverage %>% filter(Strain!="3Hybrid")

joined <- rbind(nohybrid,hybrid_noNQO,hybrid_NQO)
joined$Replicate<-as.numeric(joined$Replicate)

Supp_coverage<-joined 

mean(Supp_coverage$Mean_Read_Depth,na.rm = TRUE)

FigA <- Supp_coverage  %>% 
  ggplot(aes(y=Mean_Read_Depth, x=library))+
  theme_prism()+
  geom_boxplot(outlier.size=0.1,col="black",alpha=0.2)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, col="red4",fill="red4") +
  annotate("text",
           y = 125,
           x = 1,
           label = "105.5",
           family = "", fontface = 1, size=3,
           col="red4")+
  ylab("Mean read depth")+
  scale_x_discrete(breaks=c("Sample"),labels=c("All the \n samples")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text= element_text(size=9)) 

Supp_coverage$Type <- as.factor(Supp_coverage$Type)
Supp_coverage$Strain <- as.factor(Supp_coverage$Strain)

amod <- aov(Mean_Read_Depth~Type*Strain, data=Supp_coverage)
summary(amod)

legend_title <- " "
FigA1 <- Supp_coverage %>% 
  ggplot(aes(x=Type,y=Mean_Read_Depth, col=Strain, group=interaction(Strain,Type)))+
  geom_boxplot(outlier.size=0.1,col="black",aes(fill=Strain),alpha=0.2)+ 
  theme_prism() + 
  theme_bw(base_size=24) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  ylab("Mean read depth")+
  theme(legend.position="bottom",
        axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust=1,size=10)) +
  scale_x_discrete(limit =c("Ancestor","Evolved_control", "Evolved_NQO"),
                   label= c("Ancestor","Evolved in control","Evolved in UV mimetic"))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10))

leg1 <- get_legend(FigA1)

FigB <- Supp_coverage %>% 
  ggplot(aes(x=Type,y=Mean_Read_Depth, col=Strain, group=interaction(Strain,Type)))+
  geom_boxplot(outlier.size=0.1,col="black",aes(fill=Strain),alpha=0.2)+ 
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_prism() + 
  theme(axis.text.x = element_text()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Mean read depth")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(size=13),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust=1,size=13)) +
  scale_x_discrete(limit =c("Ancestor","Evolved_control", "Evolved_NQO"),
                   label= c("Ancestor","Evolved\n in control","Evolved in \n UV mimetic"))+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=11, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        axis.text= element_text(size=9))
#################
############Assemble and save Supplementary Figure 1############
FigA_label <- plot_grid(FigA, labels = c("a"),label_size=20)
FigA2_label <- plot_grid(FigB, labels = c("b"),label_size=20)

Fig1 <- plot_grid(FigA_label,FigA2_label,nrow = 1)
Fig1leg <- plot_grid(leg1,Fig1, nrow = 2, rel_heights = c(2,12))

ggsave (plot = Fig1leg, filename = "Supplementary_Fig1_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =12, dpi = 300,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.png", units = "cm", device = "png",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.jpg", units = "cm", device = "jpg",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.svg", units = "cm", device = "svg",width = 20, height =12, dpi = 1000,bg = "white")
ggsave (plot = Fig1leg, filename = "Supplementary_Fig1.pdf", units = "cm", device = "pdf",width = 20, height =12, dpi = 1000,bg = "white")
#################

############Figure Supplementary 2############
#NQO
AV<- dplyr::filter(Supp_genomics, condition=="NQO")
nq <- dplyr::filter(AV, Evolved=="Evolved_NQO")
con <- dplyr::filter(AV, Evolved=="Evolved_control")
anc <- dplyr::filter(AV, Evolved=="Ancestor")
AVB<-rbind(nq,con,anc)

#Define some functions
my_comparisons <- list( c("1.1Scer", "2.1Scer"), c("1.3Hybrid","2.3Hybrid"),
                        c("1.2Spar", "2.2Spar"))

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor"
)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

my_comparisons <- list( c("1Scer", "2Spar"), c("1Scer","3Hybrid"),
                        c("2Spar", "3Hybrid"))

legend_title <- " "

hybrid <-AVB %>% filter(Specie=="3Hybrid")
hybrid_noNQO <- hybrid %>% filter(Evolved!="Evolved_NQO")
hybrid_NQO <- hybrid %>% filter(Evolved=="Evolved_NQO")

hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=10)
hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=25)

nohybrid <-AVB %>% filter(Specie!="3Hybrid")

AVB <- rbind(nohybrid,hybrid_noNQO,hybrid_NQO)

Fig_ext1 <- AVB %>% dplyr::filter(day==2)%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  theme_prism() + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.position="top",legend.text.align = 0) +
  guides(color = guide_legend(title = "Genotype")) +
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=15),
        legend.direction="horizontal")+ylim(0,1)

legfitness <- get_legend(Fig_ext1)
daytwo <- AVB %>% dplyr::filter(day==2)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_NQO")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_A <- daytwo %>% dplyr::filter(Evolved=="Evolved_NQO")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.05", "p < 0.001", "p < 0.001"),
           family = "", fontface = 3, size=3)+
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_control")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_B <- daytwo %>% dplyr::filter(Evolved=="Evolved_control")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                                 labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p > 0.05", "p < 0.05"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Ancestor")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_C <- daytwo %>% dplyr::filter(Evolved=="Ancestor")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in UV mimetic conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                                 labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.67, yend=0.67)) + 
  geom_segment(aes(x=2, xend=3, y=0.70, yend=0.70)) + 
  geom_segment(aes(x=1, xend=3, y=0.74, yend=0.74)) +
  annotate("text", x = 1, y = 0.78, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.69, 0.72, 0.76),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p < 0.0001", "p < 0.001"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +ylim(0,0.8)+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

#CONTROL
AV<- dplyr::filter(Supp_genomics, condition=="Control")
nq <- dplyr::filter(AV, Evolved=="Evolved_NQO")
con <- dplyr::filter(AV, Evolved=="Evolved_control")
anc <- dplyr::filter(AV, Evolved=="Ancestor")
AVB<-rbind(nq,con,anc)

hybrid <-AVB %>% filter(Specie=="3Hybrid")
hybrid_noNQO <- hybrid %>% filter(Evolved!="Evolved_NQO")
hybrid_NQO <- hybrid %>% filter(Evolved=="Evolved_NQO")

hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=10)
hybrid_NQO <- hybrid_NQO %>% filter(Replicate!=25)

nohybrid <-AVB %>% filter(Specie!="3Hybrid")

AVB <- rbind(nohybrid,hybrid_noNQO,hybrid_NQO)

#Define some functions
my_comparisons <- list( c("1.1Scer", "2.1Scer"), c("1.3Hybrid","2.3Hybrid"),
                        c("1.2Spar", "2.2Spar"))

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor"
)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

my_comparisons <- list( c("1Scer", "2Spar"), c("1Scer","3Hybrid"),
                        c("2Spar", "3Hybrid"))

legend_title <- " "

Fig_extD <- AVB %>% dplyr::filter(day==2)%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  theme_prism() + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.position="top",legend.text.align = 0) +
  guides(color = guide_legend(title = "Genotype")) +
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=15),
        legend.direction="horizontal")+ylim(0,1)

legfitness <- get_legend(Fig_ext1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_NQO")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

daytwo <- AVB %>% dplyr::filter(day==2)

Fig_ext1_D <- daytwo %>% dplyr::filter(Evolved=="Evolved_NQO")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.77, yend=0.77)) + 
  geom_segment(aes(x=2, xend=3, y=0.80, yend=0.80)) + 
  geom_segment(aes(x=1, xend=3, y=0.84, yend=0.84)) +
  annotate("text", x = 1, y = 0.88, size = 3,
           label = c("p < 0.001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.79, 0.82, 0.86),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.01", "p < 0.01", "p > 0.05"),
           family = "", fontface = 3, size=3)+
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Evolved_control")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_E <- daytwo %>% dplyr::filter(Evolved=="Evolved_control")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.87, yend=0.87)) + 
  geom_segment(aes(x=2, xend=3, y=0.90, yend=0.90)) + 
  geom_segment(aes(x=1, xend=3, y=0.94, yend=0.94)) +
  annotate("text", x = 1, y = 0.98, size = 3,
           label = c("p < 0.0001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.89, 0.92, 0.96),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.0001", "p > 0.05", "p < 0.01"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ylim(0,1)

AVB1<-AVB %>%dplyr::filter(Evolved=="Ancestor")
AVB1<-AVB1 %>%dplyr::filter(day==2)

AVB1$Specie<-as.factor(AVB1$Specie)
amod <- aov(rval ~ Specie, data = AVB1)
summary(amod)
inter.test1 <- glht(amod,  mcp(Specie = "Tukey"))
summary(inter.test1)
cld(inter.test1)

Fig_ext1_F <- daytwo %>% dplyr::filter(Evolved=="Ancestor")%>%
  ggplot(aes(x=Specie,y=rval, col=Specie, group=Specie))+
  geom_boxplot(outlier.shape=NA,aes(fill=Specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  xlab("Group") + ylab("Growth rate (OD/hour) \n in control conditions") + scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                                                                                              labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=0.87, yend=0.87)) + 
  geom_segment(aes(x=2, xend=3, y=0.90, yend=0.90)) + 
  geom_segment(aes(x=1, xend=3, y=0.94, yend=0.94)) +
  annotate("text", x = 1, y = 0.98, size = 3,
           label = c("p < 0.001"),
           family = "", fontface = 3)+
  annotate("text",
           y = c(0.89, 0.92, 0.96),
           x = c(1.5, 2.5, 2),
           label = c("p < 0.001", "p < 0.01", "p > 0.5"),
           family = "", fontface = 3, size=3) + 
  guides(color = guide_legend(title = "Genotype")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_grid(.~Evolved, labeller = as_labeller(cond.labs))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +ylim(0,1)
#################
############Assemble and save Figure Supplementary 2############
Fig_Supp2<-plot_grid(Fig_ext1_C,Fig_ext1_B,Fig_ext1_A,
                    Fig_ext1_F,Fig_ext1_E,Fig_ext1_D,
                    nrow=2)
Fig_Supp2leg <- plot_grid(legfitness, Fig_Supp2, nrow=2,rel_heights = c(2,20))

ggsave (plot = Fig_Supp2leg, filename = "Supplementary_Fig2_low_quality.jpg", units = "cm", device = "jpg",width =40, height =25, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp2leg, filename = "Supplementary_Fig2.png", units = "cm", device = "png",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp2leg, filename = "Supplementary_Fig2.jpg", units = "cm", device = "jpg",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp2leg, filename = "Supplementary_Fig2.svg", units = "cm", device = "svg",width =40, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp2leg, filename = "Supplementary_Fig2.pdf", units = "cm", device = "pdf",width =40, height =25, dpi = 1000, bg = "white")
#################

############Figure Figure Supplementary 3############
#Make some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

data <- data.frame(
  specie = c(rep("1Scer", 17), rep("2Spar", 18), rep("3Hybrid", 12)),
  evolved = c(rep("Ancestor", 5), rep("Evolved_control", 5), rep("Evolved_UV", 7),
              rep("Ancestor", 4), rep("Evolved_control", 5), rep("Evolved_NQO", 9),
              rep("Ancestor", 3), rep("Evolved_control", 5), rep("Evolved_NQO", 4)),
  replicate = c(3, 5, 11, 12, 19, 1, 3, 5, 11, 18, 3, 5, 11, 18, 21, 13,
                1, 2, 4, 23, 1, 1, 4, 23, 26, 2, 2,5, 10, 11, 16, 17, 23, 24, 28,
                23, 24, 27, 20, 23, 26, 27, 28, 1, 26, 27, 28)
)

data <- data.frame(
  specie = rep("1Scer", 25),
  evolved = rep("Ancestor", 25),
  ploidy = rep("Diploid", 25),
  line = c(1, 2, 4, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
)

data <- data.frame(
  ploidy = c(rep("Polyploid", 5), rep("Diploid", 25),
             rep("Polyploid", 5), rep("Diploid", 25),
             rep("Polyploid", 6), rep("Diploid", 24),
             
             rep("Polyploid", 4), rep("Diploid", 26),
             rep("Polyploid", 5), rep("Diploid", 25),
             rep("Polyploid", 6), rep("Diploid", 24),
             
             rep("Polyploid", 3), rep("Diploid", 27),
             rep("Polyploid", 5), rep("Diploid", 25),
             rep("Polyploid", 4), rep("Diploid", 26)),
  specie = rep(c("1Scer", "2Spar", "3Hybrid"), each = 90),
  evolved = c(rep("Ancestor", 30), rep("Evolved_control", 30), rep("Evolved_NQO", 30),
              rep("Ancestor", 30), rep("Evolved_control", 30), rep("Evolved_NQO", 30),
              rep("Ancestor", 30), rep("Evolved_control", 30), rep("Evolved_NQO", 30)),
  line = c(3, 5, 11, 12, 19,
           1, 2, 4,  6, 7, 8, 9, 10, 13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,
           1, 3, 5, 11, 18,
           2, 4,  6, 7, 8, 9, 10, 12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30,
           3, 5, 11,13,18, 21,
           1,2, 4,6, 7, 8, 9, 10, 12,14,15,16,17,19,20,22,23,24,25,26,27,28,29,30,
           1, 2, 4, 23,
           3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30,
           1, 2, 4, 23, 26,
           3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,27,28,29,30,
           2, 5, 16, 23, 24, 28,
           1,3,4,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,25,26,27,29,30,
           23, 24, 27,
           1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,29,21,22,25,26,28,29,30,
           20, 23, 26, 27, 28, 
           1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,24,25,29,30,
           1, 26, 27, 28,
           2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,29,30))

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

data_nohyb <- filter(data, specie!="3Hybrid")
data_hyb <- filter(data, specie=="3Hybrid")
data_hyb_nonq <- filter(data_hyb, evolved!="Evolved_NQO")
data_hyb_nq <- filter(data_hyb, evolved=="Evolved_NQO")
data_hyb_nq <- filter(data_hyb_nq, line!=10)
data_hyb_nq <- filter(data_hyb_nq, line!=25)

data<-rbind(data_nohyb,data_hyb_nonq,data_hyb_nq)

count <- data %>%
  dplyr::group_by(specie, evolved,ploidy) %>%
  dplyr::summarise(n_lines = n())

Fig3A <- count %>%
  filter(evolved != "Ancestor") %>%
  ggplot(aes(x = reorder(interaction(specie, evolved, ploidy)), y = n_lines, fill = specie)) +
  geom_bar_pattern(
    stat = 'identity',
    width = 0.9,
    color = "black",
    alpha = 0.7,
    aes(pattern = evolved),
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.05,
    pattern_alpha = 0.25
  ) +
  scale_pattern_manual(" ", 
                       values = c("Evolved_control" = "none", "Evolved_NQO" = "stripe"),
                       labels = c("Evolved_control" = "Evolved in control", "Evolved_NQO" = "Evolved in UV mimetic")) +
  scale_fill_manual(" ", values = c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  xlab("Ploidy") +
  ylab("Count") +
  theme_light() +
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85"
  ) +
  scale_x_discrete(labels = c("Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid")) +
  geom_text(aes(label = n_lines), position = position_stack(vjust = 0.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.box = "vertical",
    legend.text = element_text(size=12),
    legend.key.width = unit(1, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = NA))) +  # Asegurarse de que la leyenda de `pattern` no interfiera con el `fill`
  geom_vline(xintercept = c(4.5, 8.5))

legheatmap0 <- get_legend(Fig3A)

Fig3A <- count %>%
  filter(evolved != "Ancestor") %>%
  ggplot(aes(x = reorder(interaction(specie, evolved, ploidy)), y = n_lines, fill = specie)) +
  geom_bar_pattern(
    stat = 'identity',
    width = 0.9,
    color = "black",
    alpha = 0.7,
    aes(pattern = evolved),
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.05,
    pattern_alpha = 0.25
  ) +
  scale_pattern_manual(" ", 
                       values = c("Evolved_control" = "none", "Evolved_NQO" = "stripe"),
                       labels = c("Evolved_control" = "Evolved in control", "Evolved_NQO" = "Evolved in UV mimetic")) +
  scale_fill_manual(" ", values = c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  xlab(" ") +
  ylab("Count") +
  theme_light() +
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85"
  ) +
  scale_x_discrete(labels = c("Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid")) +
  geom_text(aes(label = n_lines), position = position_stack(vjust = 0.5)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none") +
  guides(fill = guide_legend(override.aes = list(pattern = "none")),
         pattern = guide_legend(override.aes = list(fill = NA))) +  # Asegurarse de que la leyenda de `pattern` no interfiera con el `fill`
  geom_vline(xintercept = c(4.5, 8.5))

Ext2_A<-plot_grid(legheatmap0,Fig3A,nrow=2,rel_heights = c(0.5,2))

#related to fitness
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

data1<- data %>% filter(evolved=="Evolved_NQO")

data2<- full_join(data1,genomics4,by=c("line"="Replicate","specie"="Specie"))
data2 <- na.omit(data2)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Fig_Fitness <- data2 %>%
  ggplot(aes(x = interaction(ploidy,specie), y = percentage, fill = specie)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(" ", values = c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  xlab("Species and Ploidy") +
  ylab("Fitness") +
  theme_light() +
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85"
  ) +
  scale_x_discrete(labels = c("Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) +
  geom_vline(xintercept = c(2.5, 4.5))

#The interaction is not significant
data2$specie<-as.factor(data2$specie)
data2$ploidy<-as.factor(data2$ploidy)
amod <- aov(percentage ~ specie*ploidy, data = data2)
summary(amod)

inter.test1 <- glht(amod,  mcp(ploidy = "Tukey"))
summary(inter.test1)
cld(inter.test1)

data2$specie_ploidy <- interaction(data2$specie, data2$ploidy)
amod2 <- aov(percentage ~ specie_ploidy, data = data2)

inter.test <- glht(amod2, linfct = mcp(specie_ploidy = "Tukey"))
summary(inter.test)

Ext2_B<- data2 %>% 
  ggplot(aes(x=interaction(ploidy,specie),y=percentage, col=specie, group=interaction(ploidy,specie)))+
  geom_boxplot(aes(fill=specie),col="black",alpha=0.2)+ 
  geom_point(pch=21, aes(fill=specie), col="black", show.legend = F, alpha=0.5,size=2)+
  scale_x_discrete(labels = c("Diploid", "Polyploid",
                              "Diploid", "Polyploid",
                              "Diploid", "Polyploid")) +
  xlab("Group") + ylab("Increase in groth rate (%) \n in UV mimetic conditions") + 
  scale_color_manual(values=c("green4", "dodgerblue1", "black"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_fill_manual(" ",values=c("green4", "dodgerblue1", "#FF9999","grey"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_bw() +
  geom_segment(aes(x=1, xend=2, y=200, yend=200), show.legend = FALSE) + 
  geom_segment(aes(x=3, xend=4, y=200, yend=200),show.legend = FALSE) + 
  geom_segment(aes(x=5, xend=6, y=200, yend=200),show.legend = FALSE) +
  annotate("text", x = 1.5, y = 206, size = 3,
           label = c("p > 0.05"),
           family = "", fontface = 3)+
  annotate("text", x = 3.5, y = 206, size = 3,
           label = c("p > 0.05"),
           family = "", fontface = 3)+
  annotate("text", x = 5.5, y = 206, size = 3,
           label = c("p > 0.05"),
           family = "", fontface = 3)+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "top",
        legend.direction = "horizontal",
  )
#################
############Assemble and save Figure Supplementary 3############
Fig_Extended_2A_label<- plot_grid(Ext2_A,labels="a",label_size=15)
Fig_Extended_2B_label<- plot_grid(Ext2_B,labels="b",label_size=15)

Fig_Supp3_ploidy<- plot_grid(Fig_Extended_2A_label,Fig_Extended_2B_label,nrow=2,rel_heights = c(1,1))

ggsave (plot = Fig_Supp3_ploidy, filename = "Supplementary_Fig3_low_quality.jpg", units = "cm", device = "jpg",width =15, height =22, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp3_ploidy, filename = "Supplementary_Fig3.png", units = "cm", device = "png",width =15, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp3_ploidy, filename = "Supplementary_Fig3.jpg", units = "cm", device = "jpg",width =15, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp3_ploidy, filename = "Supplementary_Fig3.svg", units = "cm", device = "svg",width =15, height =22, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp3_ploidy, filename = "Supplementary_Fig3.pdf", units = "cm", device = "pdf",width =15, height =22, dpi = 1000, bg = "white")
#################

############Figure Supplementary 4############
ploidy <- Supp_ploidy

x<-ploidy%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))

fdata <- x

newdat53n <- filter(newdat5, LMH>184)

triploides <- newdat53n %>% group_by(evolved,same_Replicate,specie) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[5]  <- "counts"
control3n<-filter(x, evolved=="3n")
control2n<-filter(x, evolved=="2n")
controles <- filter(x, evolved=="2n" | evolved=="3n")

fdata$log_FL1A = log1p(fdata$FL1_A)
control2n$log_FL1A = log1p(control2n$FL1_A)
control3n$log_FL1A = log1p(control3n$FL1_A)

legend_title <- "Sample"

conto <- Supp_ploidy_heatmap %>% filter(replicate=="Control")

neworder <- c("n","2n","3n","4n")
iris2 <- arrange(transform(conto,
                           evolved=factor(evolved,levels=neworder)),evolved)

legend_title <- "Cell count (Density)"

#Graphs control
con2<- iris2 %>% filter(evolved=="2n") %>% ggplot(aes(LMH, evolved))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue',limits=c(0,800)) +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line=element_blank())+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  ylab("Controls")+
  xlab("DNA content (Fluorescence)")+
  theme(axis.text.y=element_text(size=5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

con3<- iris2 %>%filter(evolved=="2n") %>% ggplot(aes(LMH, evolved))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue',limits=c(0,800)) +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.line=element_blank())+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  xlab("DNA content (Fluorescence)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

newdat$replicate<-as.factor(newdat$replicate)

legend_title <- "Cell count (Density)"

cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor",
  `1Scer` = "S. cerevisiae",
  `2Spar` = "S. paradoxus",
  `3Hybrid` = "Hybrid"
)

Supp_ploidy_heatmap<-Supp_ploidy_heatmap %>% filter(replicate!="Control")
Supp_ploidy_heatmap$specie <- factor(Supp_ploidy_heatmap$specie,
                                     labels = c("bolditalic(S.cerevisiae)",
                                                "bolditalic(S.paradoxus)",
                                                "bold(Hybrid)"))

Supp_ploidy_heatmap$evolved <- factor(Supp_ploidy_heatmap$evolved,
                                      labels = c("bold(Ancestor)",
                                                 "bold(Control)",
                                                 "bold(UV_mimetic)"))

Supp_ploidy_heatmap <- Supp_ploidy_heatmap %>%
  filter(!(evolved == "bold(UV_mimetic)" & same_rep == "3Hybrid_25"))

Supp_ploidy_heatmap <- Supp_ploidy_heatmap %>%
  filter(!(evolved == "bold(UV_mimetic)" & same_rep == "3Hybrid_10"))

Samples_heatmap_leg <- Supp_ploidy_heatmap %>% 
  ggplot(aes(LMH, replicate))+
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue') +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.title = element_text(angle = 0, vjust= 0.3,hjust= 0.01),legend.title.align=0.5,
        legend.box.just = "center",
        legend.direction = "vertical",
        legend.position = "right")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.8, size=12, face = "bold"))+
  facet_grid(evolved~specie, scales= "free", space = "free",labeller=(label_parsed))+
  xlim(300, 500)

legheatmap <- get_legend(Samples_heatmap_leg)

Samples_heatmap <- Supp_ploidy_heatmap %>% filter(replicate!="Control") %>% ggplot(aes(LMH, replicate))+#,label = interaction(specie,info,glycerol), position = "identity")) + 
  geom_tile(aes(fill = counts)) + scale_fill_gradient(legend_title,low = 'white', high = 'darkblue') +
  ylab("Samples") +
  xlab("Fluorescence")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  )+
  theme(axis.text = element_text(face="bold"))+
  theme(strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  xlab("DNA content (Arb. units fluorescence)")+
  theme(axis.text.y=element_text(size=5))+
  theme(legend.position="none")+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=3)+
  guides(fill = guide_colorbar(title.position = "top"))+
  theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  scale_y_discrete(name ="Line", 
                   breaks=c("1","2","3","4","5","6",
                            "7","8","9","10","11",
                            "12","13","14","15","16",
                            "17","18","19","20","21",
                            "22","23","24","25","26",
                            "27","28","29","30"),
                   limits=c("1","2","3","4","5","6",
                            "7","8","9","10","11",
                            "12","13","14","15","16",
                            "17","18","19","20","21",
                            "22","23","24","25","26",
                            "27","28","29","30"),
                   label=c("1","2","3","4","5","6",
                           "7","8","9","10","11",
                           "12","13","14","15","16",
                           "17","18","19","20","21",
                           "22","23","24","25","26",
                           "27","28","29","30"))+
  facet_grid(evolved~specie, scales= "free", space = "free",labeller=label_parsed)+
  xlim(300, 500)+
  geom_vline(xintercept = c(392, 420), linetype = "dashed", color = "black",size=0.5)

Samples_heatmap2 <- Samples_heatmap +  theme(panel.spacing.x = unit(0.52, "cm"))
#################
############Assemble and save Supplementary Figure 4############
leg <- plot_grid(con2,con3,con3, ncol = 4,rel_widths = c(13.6,12.2,12,1.6))
Fig4SuppA <- plot_grid(leg,Samples_heatmap2, nrow = 2, rel_heights = c(0.8,10))
Fig4SuppAleg <- plot_grid(Fig4SuppA,legheatmap, nrow = 1, rel_widths  = c(10,2))

ggsave (plot = Fig4SuppAleg, filename = "Supplementary_Fig4_low_quality.jpg", units = "cm", device = "jpg",width = 25, height = 25, dpi = 300,bg = "white")
ggsave (plot = Fig4SuppAleg, filename = "Supplementary_Fig4.png", units = "cm", device = "png",width = 25, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig4SuppAleg, filename = "Supplementary_Fig4.jpg", units = "cm", device = "jpg",width = 25, height = 25, dpi = 1000,bg = "white")
ggsave (plot = Fig4SuppAleg, filename = "Supplementary_Fig4.svg", units = "cm", device = "svg",width = 25, height =25, dpi = 1000,bg = "white")
ggsave (plot = Fig4SuppAleg, filename = "Supplementary_Fig4.pdf", units = "cm", device = "pdf",width = 25, height =25, dpi = 1000,bg = "white")
#################

############Figure Supplementary 5############
###Boxplot of aneuploidies
#Define some elements for the graph
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Fisher test
scer <- dplyr:::filter(Supp_Table_counts_aneuploidies,Strain=="1Scer")
spar <- dplyr:::filter(Supp_Table_counts_aneuploidies,Strain=="2Spar")
hyb <- dplyr:::filter(Supp_Table_counts_aneuploidies,Strain=="3Hybrid")

#For S.cerevisiae
#        With aneu    Without aneu
#control    0               29
#NQO        3               26
m_scer <- matrix(c(0,3,29,26),nrow = 2,ncol=2,byrow=FALSE)
fisher.test(m_scer)
#p-value: 0.2368

#For S.paradoxus
#        With aneu    Without aneu
#control   0               30
#NQO       6               24
m_spar <- matrix(c(0,6,30,24),nrow = 2,ncol=2,byrow=FALSE)
fisher.test(m_spar)
#p-value: 0.02372

#For Hybrid
#        With aneu    Without aneu
#control   1               25
#NQO       4               22
m_hyb <- matrix(c(1,4,25,22),nrow = 2,ncol=2,byrow=FALSE)
fisher.test(m_hyb)
#p-value: 0.3497

B_Table_counts_aneuploidies_filtered <- Supp_Table_counts_aneuploidies %>%
  dplyr:::filter(!(Strain=="1Scer" & Replicate==16))
B_Table_counts_aneuploidies_filtered <- B_Table_counts_aneuploidies_filtered %>%
  dplyr:::filter(!(Strain=="3Hybrid" & Replicate==1))
B_Table_counts_aneuploidies_filtered <- B_Table_counts_aneuploidies_filtered %>%
  dplyr:::filter(!(Strain=="3Hybrid" & Replicate==10))
B_Table_counts_aneuploidies_filtered <- B_Table_counts_aneuploidies_filtered %>%
  dplyr:::filter(!(Strain=="3Hybrid" & Replicate==21))
B_Table_counts_aneuploidies_filtered <- B_Table_counts_aneuploidies_filtered %>%
  dplyr:::filter(!(Strain=="3Hybrid" & Replicate==25))

FigSupp5 <- ggplot(B_Table_counts_aneuploidies_filtered, aes(x = interaction(Type,Strain), y = counts,fill=Strain))+
  geom_boxplot(outlier.shape = NA,aes(fill=as.factor(Strain)))+
  geom_jitter(colour="black",pch=21,height = 0, size=3,alpha=0.7, aes(fill=as.factor(Strain)))+
  ylab("Number of \n aneuploidies / line")+
  xlab("Genotype")+
  scale_fill_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_x_discrete("", labels=c("Evolved in \n control","Evolved in \n UV mimetic",
                                "Evolved in \n control","Evolved in \n UV mimetic",
                                "Evolved in \n control","Evolved in \n UV mimetic"))+
  theme_bw(base_size=24) +
  geom_segment(aes(x=3, xend=4, y=1.3, yend=1.3)) +
  annotate("text",
           y = c(1.45),
           x = c(3.5),
           label = c("p < 0.05"),
           family = "", fontface = 3, size=4) +
  scale_y_continuous(breaks = c(0, 1, 2))+
  theme_bw(base_size=24) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(size=15, hjust=0, face = "bold", color="black"),
        axis.text.x = element_text(size = 11, color="black",face = "bold"),
        axis.text.y = element_text(size = 13, color="black"),
        axis.title = element_text(size = 15,face="bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.margin = margin(-10, -10,-6,-2))+
  guides(fill = guide_legend(title = " "))

#################
############Assemble and save Supplementary Figure 5############
ggsave (plot = FigSupp5, filename = "Supplementary_Fig5_low_quality.jpg", units = "cm", device = "jpg",width = 25, height = 10, dpi = 300,bg = "white")
ggsave (plot = FigSupp5, filename = "Supplementary_Fig5.png", units = "cm", device = "png",width = 25, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp5, filename = "Supplementary_Fig5.jpg", units = "cm", device = "jpg",width = 25, height = 10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp5, filename = "Supplementary_Fig5.svg", units = "cm", device = "svg",width = 25, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp5, filename = "Supplementary_Fig5.pdf", units = "cm", device = "pdf",width = 25, height =10, dpi = 1000,bg = "white")
#################

############Figure Supplementary 6############
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Supp_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

Supp_Table_counts_aneuploidies_NQO <- Supp_Table_counts_aneuploidies %>% 
  dplyr:::filter(Type=="Evolved_NQO")

every_aneu_NQO<- right_join(Supp_Aneuploidy_chrom,Supp_Table_counts_aneuploidies_NQO, by=c("Specie"="Strain",
                                                                                           "Replicate"="Replicate"))
FITNESS_LOH <-left_join(every_aneu_NQO,genomics4,by=c("Specie"="Specie","Replicate"="Replicate"))

#Define some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

NQO <- dplyr:::filter(FITNESS_LOH, Type=="Evolved_NQO")

#Get statistics
NQO <- dplyr:::filter(FITNESS_LOH, Type=="Evolved_NQO")
NQO_scer <- dplyr:::filter(NQO, Specie=="1Scer")
ggscatter(NQO_scer, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

NQO_spar <- dplyr:::filter(NQO, Specie=="2Spar")
ggscatter(NQO_spar, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

NQO_hyb <- dplyr:::filter(NQO, Specie=="3Hybrid")
ggscatter(NQO_hyb, x = "counts", y = "percentage",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

legend_title <- " "

FigSupp6 <- FITNESS_LOH %>% filter(Type=="Evolved_NQO")%>%
  ggplot(aes(x=as.factor(counts),y = percentage))+
  stat_smooth(method="lm", alpha=0.1,aes(group=Specie,col=Specie))+
  geom_jitter(pch=21, aes(fill=Specie), col="black", show.legend = F, alpha=0.5,size=2,width = 0.1, height = 0.1)+
  xlab("Number of aneuploidies / line")+ylab("Increase in growth rate (%)")+
  border()  +
  scale_fill_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  scale_color_manual(legend_title,values=c("green4","dodgerblue1","#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus","Hybrid"))) +
  theme(legend.position="none")+
  stat_smooth(method="lm", alpha=0.1)+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(-15, -15,-15,-15))+
  annotate("text",
           y = 190,x = 3,
           label = c("rs = -0.18, p > 0.05"),
           family = "", fontface = 3, size=3, col="darkgreen") +
  annotate("text",
           y = 180,x = 3,
           label = c("rs = -0.22, p > 0.05"),
           family = "", fontface = 3, size=3, col="dodgerblue3") +
  annotate("text",
           y = 170,x = 3,
           label = c("rs = -0.27, p > 0.05"),
           family = "", fontface = 3, size=3, col="lightpink4")+
  stat_smooth(method="lm", alpha=0.2) 
#################
############Assemble and save Supplementary Figure 6############
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6_low_quality.jpg", units = "cm", device = "jpg",width = 15, height =10, dpi = 300,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.png", units = "cm", device = "png",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.jpg", units = "cm", device = "jpg",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.svg", units = "cm", device = "svg",width = 15, height =10, dpi = 1000,bg = "white")
ggsave (plot = FigSupp6, filename = "Supplementary_Fig6.pdf", units = "cm", device = "pdf",width = 15, height =10, dpi = 1000,bg = "white")
#################

############Figure Figure Supplementary 7############
hybrid_nq <-  dplyr:::filter(Supp_genome, Strain=="3Hybrid")
hybrid_nq <-  dplyr:::filter(hybrid_nq, Type=="Evolved_NQO")

exampleHybrid <- dplyr:::filter(hybrid_nq, Strain=="3Hybrid")
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=1)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=21)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=10)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=25)

nq_snps_HYBRID <-  exampleHybrid

#Normalize max and min values
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read >1.5, 1.5,Normalized_read))
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read < -1.5, -1.5,Normalized_read))

nq_genome<- nq_snps_HYBRID
legend_title <- "Relative Read Depth"

Fig_Supp7 <-nq_genome %>% dplyr:::filter(Type=="Evolved_NQO") %>% ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Normalized_read)) + 
  theme_prism() + 
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_tile(aes(fill =Normalized_read),height=0.65,size=0.2) +
  scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "lavender", high = "orange3", midpoint = 0, limits = c(-2.9, 2.9),breaks= c(-2,-1,0,1,2))+
  scale_y_discrete(limits=c("32","31","30","29","28","27","26","25","24","23","22","21","20","19","18","17",
                            "16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"),
                   labels=c("Spar XVI","Scer chrXVI","Spar chrXV","Scer chrXV",
                            "Spar chrXIV","Scer chrXIV","Spar chrXIII","Scer chrXIII",
                            "Spar chrXII","Scer chrXII","Spar chrXI","Scer chrXI",
                            "Spar chrX","Scer chrX","Spar chrIX","Scer chrIX",
                            "Spar chrVIII","Scer chrVIII","Spar chrVII","Scer chrVII",
                            "Spar chrVI","Scer chrVI","Spar chrV","Scer chrV",
                            "Spar chrIV","Scer chrIV","Spar chrIII","Scer chrIII",
                            "Spar chrII","Scer chrII","Spar chrI","Scer chrI"))+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  facet_wrap(.~ Replicate, ncol=13)+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome") +
  ggtitle("Hybrids evolved in UV mimetic conditions")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 40, face = "bold"))+
  theme(legend.background = element_rect(fill="snow2",size=0.4, linetype="solid",colour ="grey"))+
  theme(legend.key.size = unit(2, 'cm'),
        legend.key.height = unit(2, 'cm'),
        legend.key.width = unit(2, 'cm'),
        legend.title = element_text(size=25,angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1,
        legend.text = element_text(size=25))
#################
############Assemble and save Figure Supplementary 7############
ggsave (plot = Fig_Supp7, filename = "Supplementary_Fig7_low_quality.jpg", units = "cm", device = "jpg",width =95, height =50, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp7, filename = "Supplementary_Fig7.png", units = "cm", device = "png",width =95, height =50, dpi = 800, bg = "white")
ggsave (plot = Fig_Supp7, filename = "Supplementary_Fig7.jpg", units = "cm", device = "jpg",width =95, height =50, dpi = 800, bg = "white")
ggsave (plot = Fig_Supp7, filename = "Supplementary_Fig7.svg", units = "cm", device = "svg",width =95, height =50, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp7, filename = "Supplementary_Fig7.pdf", units = "cm", device = "pdf",width =95, height =50, dpi = 1000, bg = "white")
#################

############Figure Supplementary 8############
#Figure Supp8
size_chrom <- Supp_genome %>% group_by(chrom3, Replicate,Strain,Type) %>% dplyr::summarise(Size_chrom = last(end))
size_chrom <- size_chrom %>% dplyr:::filter(Type!="Ancestor")

LOH_events_counts_correlation<- Supp_LOH_chrom %>% dplyr:::filter(LOH=="LOH") %>% group_by(Type,LOH,chrom3) %>%   
  dplyr:::summarise(n()) %>% ungroup()
colnames(LOH_events_counts_correlation)[4]  <- "counts"

LOH_events_counts_correlation$Strain <- "3Hybrid"
LOH_events_counts_correlation$Strain <- "3Hybrid"

corr_aneu5a<-LOH_events_counts_correlation  %<>% mutate(same_chrom = ifelse(chrom3=="Scer chrI",1,
                                                                            ifelse(chrom3=="Scer chrII",2,
                                                                                   ifelse(chrom3=="Scer chrIII",3,
                                                                                          ifelse(chrom3=="Scer chrIV",4,
                                                                                                 ifelse(chrom3=="Scer chrV",5,
                                                                                                        ifelse(chrom3=="Scer chrVI",6,
                                                                                                               ifelse(chrom3=="Scer chrVII",7,
                                                                                                                      ifelse(chrom3=="Scer chrVIII",8,
                                                                                                                             ifelse(chrom3=="Scer chrIX",9,
                                                                                                                                    ifelse(chrom3=="Scer chrX",10,
                                                                                                                                           ifelse(chrom3=="Scer chrXI",11,
                                                                                                                                                  ifelse(chrom3=="Scer chrXII",12,
                                                                                                                                                         ifelse(chrom3=="Scer chrXIII",13,
                                                                                                                                                                ifelse(chrom3=="Scer chrXIV",14,
                                                                                                                                                                       ifelse(chrom3=="Scer chrXV",15,
                                                                                                                                                                              ifelse(chrom3=="Scer chrXVI",16,
                                                                                                                                                                                     ifelse(chrom3=="Spar chrI",1,
                                                                                                                                                                                            ifelse(chrom3=="Spar chrII",2,
                                                                                                                                                                                                   ifelse(chrom3=="Spar chrIII",3,
                                                                                                                                                                                                          ifelse(chrom3=="Spar chrIV",4,
                                                                                                                                                                                                                 ifelse(chrom3=="Spar chrV",5,
                                                                                                                                                                                                                        ifelse(chrom3=="Spar chrVI",6,
                                                                                                                                                                                                                               ifelse(chrom3=="Spar chrVII",7,
                                                                                                                                                                                                                                      ifelse(chrom3=="Spar chrVIII",8,
                                                                                                                                                                                                                                             ifelse(chrom3=="Spar chrIX",9,
                                                                                                                                                                                                                                                    ifelse(chrom3=="Spar chrX",10,
                                                                                                                                                                                                                                                           ifelse(chrom3=="Spar chrXI",11,
                                                                                                                                                                                                                                                                  ifelse(chrom3=="Spar chrXII",12,
                                                                                                                                                                                                                                                                         ifelse(chrom3=="Spar chrXIII",13,
                                                                                                                                                                                                                                                                                ifelse(chrom3=="Spar chrXIV",14,
                                                                                                                                                                                                                                                                                       ifelse(chrom3=="Spar chrXV",15,
                                                                                                                                                                                                                                                                                              ifelse(chrom3=="Spar chrXVI",16,NA)))))))))))))))))))))))))))))))))
corr_aneu5b <- ddply(corr_aneu5a, .(same_chrom), summarise, sum_counts=sum(counts))

corr_aneu5b<-corr_aneu5b %<>% mutate(same_chrom2 = ifelse(same_chrom==3,"chrIII",
                                                          ifelse(same_chrom==7,"chrVII",
                                                                 ifelse(same_chrom==10,"chrX",
                                                                        ifelse(same_chrom==11,"chrXI",
                                                                               ifelse(same_chrom==12,"chrXII",
                                                                                      ifelse(same_chrom==14,"chrXIV",
                                                                                             ifelse(same_chrom==15,"chrXV",
                                                                                                    ifelse(same_chrom==16,"chrXVI",NA)))))))))
corr_aneu5b <- corr_aneu5b %>% mutate(Size_chrom2=ifelse(same_chrom==15,1048000,
                                                         ifelse(same_chrom==12,1003000,
                                                                ifelse(same_chrom==14,742000,
                                                                       ifelse(same_chrom==11,643000,
                                                                              ifelse(same_chrom==3,295000,
                                                                                     ifelse(same_chrom==7,1058000,
                                                                                            ifelse(same_chrom==10,699000,
                                                                                                   ifelse(same_chrom==16,908000,NA)))))))))
corr_aneu_5 <- corr_aneu5b
cond.labs <- c(
  `Evolved_control` = "Evolved in control",
  `Evolved_NQO` = "Evolved in UV mimetic",
  `Ancestor` = "Ancestor",
  `1Scer` = "Saccharomyces cerevisiae",
  `2Spar` = "Saccharomyces paradoxus",
  `3Hybrid` = "Hybrid"
)
corr_aneu_5$Strain <- "3Hybrid" 

ggscatter(corr_aneu_5, x = "Size_chrom2", y = "sum_counts",
          color = "black", shape = 21, size = 3,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), 
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))


Fig_Supp8_A <-corr_aneu_5 %>%filter(sum_counts>0)%>%
  ggplot(aes(x=Size_chrom2,y = sum_counts))+
  stat_smooth(method="lm", alpha=0.1)+
  annotate("text",
           y = 8,x =300000,
           label = c("chrIII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y =8,x =643000,
           label = c("chrXI"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = -5,x =699000,
           label = c("chrX"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 16,x =742000,
           label = c("chrXIV"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 10,x =908000,
           label = c("chrXVI"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 46,x =1003000,
           label = c("chrXII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 14,x =1058000,
           label = c("chrVII"),
           family = "", fontface = 3, size=2.2)+
  annotate("text",
           y = 34,x =1048000,
           label = c("chrXV"),
           family = "", fontface = 3, size=2.2)+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, aes(fill=Strain), col="black", show.legend = F, alpha=0.5,size=2)+
  xlab("Chromosome Size (kbp)")+ylab("Total number of LOH")+
  border()+
  scale_x_continuous(breaks=c(400000,600000,800000,1000000),
                     labels=c(400,600,800,1000))+
  scale_shape_manual("",values=c(1,19),
                     label=c("Evolved in control","Evolved in UV mimetic"))+theme(strip.text.x = element_blank())+
  theme(panel.spacing = unit(0.72, "cm"))+
  guides(fill=guide_legend(title=" "))+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position="none",
        axis.title = element_text(size=10, face = "bold"),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.background = element_blank())+
  annotate("text",
           y = 40,x = 450000,
           label = c("rs = 0.81, p < 0.05"),
           family = "", fontface = 3, size=3.5) 

#Figure Ext4_B
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics2<-genomics1 %>% dplyr:::filter(Evolved=="Evolved_NQO")
genomics1 <- dplyr:::filter(Extended_genomics,condition=="NQO")
genomics1 <- dplyr:::filter(genomics1,day==2)
genomics1$Replicate<-as.numeric(genomics1$Replicate)
genomics3<-genomics1 %>% dplyr:::filter(Evolved=="Ancestor")

genomics4<- left_join(genomics2,genomics3, by=c("Specie"="Specie","condition"="condition","Replicate"="Replicate"))

genomics4$fitness_gain<- (genomics4$rval.x - genomics4$rval.y) / genomics4$rval.y
genomics4$percentage <- genomics4$fitness_gain * 100

newloh <-Supp_counts_LOH

newloh$Replicate <- as.numeric(newloh$Replicate)

FITNESS_LOH <-left_join(newloh,genomics4,by=c("Type"="Evolved.x","Replicate"="Replicate"))

#Some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

FITNESS_LOH_nq<- FITNESS_LOH %>% dplyr::filter(Type=="Evolved_NQO")

FITNESS_LOH_nq<- FITNESS_LOH_nq %>% dplyr::filter(Replicate!=25)
FITNESS_LOH_nq<- FITNESS_LOH_nq %>% dplyr::filter(Replicate!=10)

FITNESS_LOH_nq_hyb<-dplyr::filter(FITNESS_LOH_nq,Specie=="3Hybrid")
FITNESS_LOH_nq_hyb$new_counts <- FITNESS_LOH_nq_hyb$counts /2

ggscatter(FITNESS_LOH_nq_hyb, x = "new_counts", y = "percentage",
          color = "black", shape = 21, size = 3,
          add = "reg.line",  
          add.params = list(color = "blue", fill = "lightgray"), 
          cor.coef = TRUE,
          cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"))

Fig_Supp8_B <- FITNESS_LOH_nq_hyb %>% dplyr::filter(Type=="Evolved_NQO")%>%
  ggplot(aes(x=new_counts,y = percentage))+
  stat_smooth(method="lm", alpha=0.1)+
  annotate("text",
           y = 180,x = 0.6,
           label = c("rs = 0.33, p > 0.05"),
           family = "", fontface = 3, size=3.5) + 
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, aes(fill=Type), col="black", show.legend = F, alpha=0.5,size=2)+
  xlab("Number of LOH / line")+ylab("Increase in growth rate (%)")+
  border()  +
  scale_fill_manual(values=c("#FF9999"), 
                    labels = toexpr(c("Hybrid"))) +
  scale_color_manual(values=c("#FF9999"), 
                     labels = toexpr(c("Hybrid"))) +
  theme(legend.position="none")+
  stat_smooth(method="lm", alpha=0.1)+
  theme_prism()+
  theme_bw(base_size=24) +
  theme(legend.position="none",
        axis.title = element_text(size=10, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        panel.background = element_blank()) +theme(strip.text.x = element_blank())
#################
############Assemble and save Extended Figure 4############
Fig8A_label<- plot_grid(Fig_Supp8_A,labels="a",label_size=15)
Fig8B_label<- plot_grid(Fig_Supp8_B,labels="b",label_size=15)

Fig8_supp_LOH<-plot_grid(Fig8A_label,Fig8B_label)

ggsave (plot = Fig8_supp_LOH, filename = "Supplementary_Fig8_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =8, dpi = 300, bg = "white")
ggsave (plot = Fig8_supp_LOH, filename = "Supplementary_Fig8.png", units = "cm", device = "png",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig8_supp_LOH, filename = "Supplementary_Fig8.jpg", units = "cm", device = "jpg",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig8_supp_LOH, filename = "Supplementary_Fig8.svg", units = "cm", device = "svg",width = 20, height =8, dpi = 1000, bg = "white")
ggsave (plot = Fig8_supp_LOH, filename = "Supplementary_Fig8.pdf", units = "cm", device = "pdf",width = 20, height =8, dpi = 1000, bg = "white")
#################

############Figure Supplementary 9############
SNP_after_norm2<-Supp_SNPs %<>% mutate(chrom2 = ifelse(CHROM=="utg351_pilon","Spar chrI",
                                                           ifelse(CHROM=="utg1271_pilon","Spar chrII",
                                                                  ifelse(CHROM=="utg584_pilon","Spar chrIII",
                                                                         ifelse(CHROM=="utg988_pilon","Spar chrIV",
                                                                                ifelse(CHROM=="utg69_pilon","Spar chrV",
                                                                                       ifelse(CHROM=="utg639_pilon","Spar chrVI",
                                                                                              ifelse(CHROM=="utg199_pilon","Spar chrVII",
                                                                                                     ifelse(CHROM=="utg176_pilon","Spar chrVIII",
                                                                                                            ifelse(CHROM=="utg1121_pilon","Spar chrIX",
                                                                                                                   ifelse(CHROM=="utg245_pilon","Spar chrX",
                                                                                                                          ifelse(CHROM=="utg298_pilon","Spar chrXI",
                                                                                                                                 ifelse(CHROM=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                        ifelse(CHROM=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                               ifelse(CHROM=="utg210_pilon","Spar chrXIII",
                                                                                                                                                      ifelse(CHROM=="utg48_pilon","Spar chrXIV",
                                                                                                                                                             ifelse(CHROM=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                    ifelse(CHROM=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                           ifelse(CHROM=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                  ifelse(CHROM=="chrI","Scer chrI",
                                                                                                                                                                                         ifelse(CHROM=="chrII","Scer chrII",
                                                                                                                                                                                                ifelse(CHROM=="chrIII","Scer chrIII",
                                                                                                                                                                                                       ifelse(CHROM=="chrIV","Scer chrIV",
                                                                                                                                                                                                              ifelse(CHROM=="chrV","Scer chrV",
                                                                                                                                                                                                                     ifelse(CHROM=="chrVI","Scer chrVI",
                                                                                                                                                                                                                            ifelse(CHROM=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                   ifelse(CHROM=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                          ifelse(CHROM=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                 ifelse(CHROM=="chrX","Scer chrX",
                                                                                                                                                                                                                                                        ifelse(CHROM=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                               ifelse(CHROM=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                      ifelse(CHROM=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                             ifelse(CHROM=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                    ifelse(CHROM=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                           ifelse(CHROM=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))
SNP_after_norm2$sample <- paste(SNP_after_norm2$Strain,SNP_after_norm2$Replicate,SNP_after_norm2$condition, sep="")
SNP_after_norm2_nooutliers <- SNP_after_norm2 %>% dplyr::filter(sample!="1Scer10NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar1NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar3NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar5NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar7NQO")
SNP_after_norm2_nooutliers <- SNP_after_norm2_nooutliers %>% dplyr::filter(sample!="2Spar10NQO")

#Get some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

all2<-SNP_after_norm2_nooutliers

all2$Type <- all2$condition

count_data <- all2 %>%
  group_by(Strain,Type, Replicate) %>%
  mutate(count = n())

rows.per.group  <- aggregate(rep(1, length(paste0(all2$Strain, all2$Type,all2$Replicate))),
                             by=list(all2$Strain, all2$Type,all2$Replicate), sum)

colnames(rows.per.group) <- c("Strain","Type","Replicate","Counts")

rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Mean=mean(Counts), Max=max(Counts), Min=min(Counts), Median=median(Counts), Std=sd(Counts))

p_meds <- rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Median=median(Counts))

#Remove some lines
all2_2_nohyb<-dplyr::filter(rows.per.group,Strain!="3Hybrid")
all2_2_hyb<-dplyr::filter(rows.per.group,Strain=="3Hybrid")
all2_2_hyb_con<-dplyr::filter(all2_2_hyb,Type!="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb,Type=="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=1)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=10)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=21)
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=25)

rows.per.groupB<-rbind(all2_2_nohyb,all2_2_hyb_con,all2_2_hyb_nq)

rows.per.group<-rows.per.groupB

p_meds <- rows.per.group%>%
  group_by(Type,Strain)%>% 
  summarise(Median=median(Counts))

rows.per.group_nqo<- dplyr::filter(rows.per.group, Type=="NQO")

p_meds <- rows.per.group_nqo%>%
  group_by(Type,Strain)%>% 
  dplyr::summarise(Median=median(Counts))

my_comparisons <- list( c("NQO.1Scer", "NQO.2Spar"), c("NQO.1Scer", "NQO.3Hybrid"), c("NQO.2Spar", "NQO.3Hybrid") )

legend_title <- " "
Fig9SuppA_1 <- rows.per.group %>% dplyr:::filter(Type=="NQO") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in UV mimetic")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position="none")+
  ylim(0,450)

rows.per.group_nqo<- dplyr::filter(rows.per.group, Type=="control")

p_meds <- rows.per.group_nqo%>%
  group_by(Type,Strain)%>% 
  dplyr::summarise(Median=median(Counts))

my_comparisons <- list( c("control.1Scer", "control.2Spar"), c("control.1Scer", "control.3Hybrid"), c("control.2Spar", "control.3Hybrid") )

legend_title<- " "
Fig9SuppA_2 <- rows.per.group %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in control")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position = "top",
        legend.direction = "horizontal")+
  ylim(0,450)+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title.align=0.1,
        legend.text = element_text(size=15))

leg <- get_legend(Fig9SuppA_2)

Fig9SuppA_2 <- rows.per.group %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=interaction(Type,Strain),y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = p_meds, aes(x = interaction(Type,Strain), y = Median, label =  Median),fontface = "bold",
             size = 3.5)+
  ggtitle("All variants in control")+
  ylab("Count")+xlab(" ")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  theme(legend.position = "none")+
  ylim(0,450)

#Missense variants
all2_2 <- dplyr::filter(Extended_GO,Consequence=="missense_variant")
all2_2_nohyb<-dplyr::filter(all2_2,Strain!="3Hybrid")
all2_2_hyb<-dplyr::filter(all2_2,Strain=="3Hybrid")
all2_2_hyb_con<-dplyr::filter(all2_2_hyb,condition!="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb,condition=="NQO")
all2_2_hyb_nq<-dplyr::filter(all2_2_hyb_nq,Replicate!=25)
all2_2_b <-rbind(all2_2_nohyb,all2_2_hyb_con,all2_2_hyb_nq)
all2_2<-all2_2_b
rows.per.group  <- aggregate(rep(1, length(paste0(all2_2$Strain, all2_2$condition,all2_2$Replicate))),
                             by=list(all2_2$Strain, all2_2$condition,all2_2$Replicate), sum)
colnames(rows.per.group) <- c("Strain","Type","Replicate","Counts")
medians <- aggregate(Counts ~  Strain*Type,rows.per.group, median)
con <- rows.per.group %>% dplyr::filter(Type=="control")
con2 <- medians %>% dplyr::filter(Type=="control")
nq <- rows.per.group %>% filter(Type=="NQO")
nq2 <- medians %>% dplyr::filter(Type=="NQO")

#Get some functions
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Fig9SuppB_1 <- con %>% dplyr::filter(Type=="control") %>%
  ggplot(aes(x=Strain,y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = con2, aes(x = Strain, y = Counts, label =  Counts),fontface = "bold",
             size = 3.5)+
  theme_bw(base_size=15) + 
  ylab("Count")+xlab(" ")+
  ggtitle("Missense variants in control")+
  theme(legend.position = "none")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  theme(legend.position="none")+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=3)+
  ylim(0,450)

Fig9SuppB_2 <- nq %>% dplyr::filter(Type=="NQO") %>%
  ggplot(aes(x=Strain,y=Counts)) +
  geom_boxplot(aes(fill=Strain))+
  geom_point(colour="black",pch=21,size=2, aes(fill=Strain))+
  background_grid(
    major = c("xy"),
    minor = c("y"),
    size.major = 0.5,
    size.minor = 0.3,
    color.major = "grey85",
    color.minor = "grey85") + ylab("Counts") +
  geom_label(data = nq2, aes(x = Strain, y = Counts, label =  Counts),fontface = "bold",
             size = 3.5)+
  theme_bw(base_size=15) + 
  ylab("Count")+xlab(" ")+
  ggtitle("Missense variants in UV mimetic")+
  theme(legend.position = "none")+
  scale_x_discrete(" ", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  scale_fill_manual(legend_title,values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid")))+
  theme_bw() +
  theme(axis.title = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 9,face = "bold"),
        strip.background = element_blank())+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  annotate("text",
           y = c(440),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=3)+
  ylim(0,450)+
  theme(legend.position = "none")

#################
############Assemble and save Figure Supplementary 9############
Figure9a<-plot_grid(Fig9SuppA_2,Fig9SuppA_1,rel_widths = c(1,1))
Figure9b<-plot_grid(Fig9SuppB_1,Fig9SuppB_2,rel_widths = c(1,1))
Fig9B_label<- plot_grid(Figure9b,labels="b",label_size=15)
Fig9A_label<- plot_grid(Figure9a,labels="a",label_size=15)
Figure9_Supp_SNP <- plot_grid(Fig9A_label,Fig9B_label, nrow = 2)
Figure9_Supp_SNPleg <- plot_grid(leg, Figure9_Supp_SNP, nrow=2,rel_heights = c(2,20))

ggsave (plot = Figure9_Supp_SNPleg, filename = "Supplementary_Fig9_low_quality.jpg", units = "cm", device = "jpg",width = 22, height = 22, dpi = 300, bg = "white")
ggsave (plot = Figure9_Supp_SNPleg, filename = "Supplementary_Fig9.png", units = "cm", device = "png",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure9_Supp_SNPleg, filename = "Supplementary_Fig9.jpg", units = "cm", device = "jpg",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure9_Supp_SNPleg, filename = "Supplementary_Fig9.svg", units = "cm", device = "svg",width = 22, height = 22, dpi = 1000, bg = "white")
ggsave (plot = Figure9_Supp_SNPleg, filename = "Supplementary_Fig9.pdf", units = "cm", device = "pdf",width = 22, height = 22, dpi = 1000, bg = "white")
#################

############Figure Supplementary 10############
#Special Packages of enrichment

#install.packages("AnnotationDbi")
library(AnnotationDbi)
#BiocManager::install("AnnotationHub")
library(AnnotationHub)
#install.packages("biomaRt")
library(biomaRt)
#install.packages("biomartr")
library(biomartr) 
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("cola")
library(cola)
#BiocManager::install("DO.db")
library(DO.db)
#install.packages("DOSE")
library(DOSE)
#install.packages("enrichplot")
library(enrichplot)
#BiocManager::install("gage")
library(gage)
#install.packages("gageData")
library(gageData)
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("gsEasy")
library(gsEasy)
#BiocManager::install("hu6800.db")
library(hu6800.db)
#BiocManager::install("org.Sc.sgd.db")
library(org.Sc.sgd.db)
#BiocManager::install("ReactomePA")
library(ReactomePA)
#BiocManager::install("simplifyEnrichment")
library(simplifyEnrichment) 

m1 <- read.csv("10Supp_GO_scer_NQO_without_outliers.csv")
m1_miss <- filter(m1, Consequence=="missense_variant")
m1 <- m1_miss

m3 <- read.csv("10Supp_GO_spar_NQO_without_outliers.csv")
m3_miss <- filter(m3, Consequence=="missense_variant")
m3 <- m3_miss

m5 <- read.csv("10Supp_GO_hybrid_NQO_without_outliers_hybrid_concatenated.csv")
m5_miss <- filter(m5, Consequence=="missense_variant")
m5 <- m5_miss

#Remove 25
m5 <- filter(m5,Replicate!=25)

m1a<- m1$ensembl_gene_id
m1a<-as.data.frame(m1a)
colnames(m1a) <- c("gene")
m1a<- m1$ensembl_gene_id
hyb_NQO <- enrichGO(m1a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)

x1 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A1<- emapplot(x1,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Scer 4-NQO")

m3a<- m3$ensembl_gene_id
m3a<-as.data.frame(m3a)
colnames(m3a) <- c("gene")
m3a<- m3$ensembl_gene_id
hyb_NQO <- enrichGO(m3a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)
x2 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A2<- emapplot(x2,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Spar 4-NQO")

m5a<- m5$ensembl_gene_id
m5a<-as.data.frame(m5a)
colnames(m5a) <- c("gene")
m5a<- m5$ensembl_gene_id
hyb_NQO <- enrichGO(m5a,"org.Sc.sgd.db", keyType = "ENSEMBL",
                    ont = "ALL",
                    qvalueCutoff = 0.5,pvalueCutoff=0.5)
x3 <- pairwise_termsim(hyb_NQO) 
set.seed(123)
A3<- emapplot(x3,color = "p.adjust",group_category = T, group_legend = T, cex_label_group = 1.5, node_label = "ONTOLOGY",nCluster = 5)+
  ggtitle("Hybrid 4-NQO")

network<- plot_grid(A1,A2,A3, nrow=1)
#################
############Assemble and save raw Supplementary Figure 10############
ggsave (plot =  network, filename = "Supplementary_Fig10_raw_low_quality.jpg", units = "cm", device = "jpg",width = 60, height =25, dpi = 400,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig10_raw.png", units = "cm", device = "png",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig10_raw.jpg", units = "cm", device = "jpg",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig10_raw.svg", units = "cm", device = "svg",width = 60, height =25, dpi = 1000,bg = "white")
ggsave (plot =  network, filename = "Supplementary_Fig10_raw.pdf", units = "cm", device = "pdf",width = 60, height =25, dpi = 1000,bg = "white")
#################

#Calculate ratio of enrichment
x1b<-as.data.frame(x1)
x1b$Specie<-"Scer"
x2b<-as.data.frame(x2)
x2b$Specie<-"Spar"
x3b<-as.data.frame(x3)
x3b$Specie<-"Hybrid"

xb<-rbind(x1b,x2b,x3b)

fraction_a_number <- function(fraction) {
  eval(parse(text = fraction))
}
View(xb)

xb$GeneRatio_number <- sapply(xb$GeneRatio, fraction_a_number)
xb$BgRatio_number <- sapply(xb$BgRatio, fraction_a_number)
xb$enrichment_Ratio <-  xb$GeneRatio_number / xb$BgRatio_number

#So the most 2 enriched GO term (higher value of this division of ratios) per specie:
#S.cerevisiae: replication-born double-strand break repair via sister chromatid exchange and regulation of cell differentiation
#S.paradoxus: ER-associated misfolded protein catabolic process, endocytic vesicle
#Hybrid: trehalose metabolic process, ABC-type transporter activity (normally involved in efflux pumps poner)

#################
############Assemble and save Supplementary Figure 10############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Supp10_GO, interpolate=TRUE)

FigSupp10<-plot_grid(gpp)

ggsave (plot = FigSupp10, filename = "Supplementary_Fig10_low_quality.jpg", units = "cm", device = "jpg",width = 20, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp10, filename = "Supplementary_Fig10.png", units = "cm", device = "png",width = 20, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp10, filename = "Supplementary_Fig10.jpg", units = "cm", device = "jpg",width = 20, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp10, filename = "Supplementary_Fig10.svg", units = "cm", device = "svg",width = 20, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp10, filename = "Supplementary_Fig10.pdf", units = "cm", device = "pdf",width = 20, height =40, dpi = 300,bg = "white")
#################

############Figure Supplementary 11############
#Positions of mutations
SNP <- c(280,280, 282,282, 298,308,308, 308,308,484,516,691,729,762,820,867,1041,1041,1047,1045)

sample.gr <- GRanges("chr7", IRanges(SNP, width=1, names=paste0("", SNP)))
features <- GRanges("chr7", IRanges(c(0),
                                    width=c(1050),
                                    names=paste0("block", 1:1)))
sample.gr$color <- sample.int(6, length(SNP), replace=TRUE)

sample.gr$color <- (c("dodgerblue1", "#FF9999",
                      "green4", "#FF9999",
                      "#FF9999", "dodgerblue1",
                      "dodgerblue1", "#FF9999",
                      "#FF9999", "dodgerblue1",
                      "#FF9999", "dodgerblue1",
                      "green4", "dodgerblue1",
                      "#FF9999", "dodgerblue1",
                      "dodgerblue1", "green4",
                      "#FF9999","#FF9999"))

plot <- lolliplot(sample.gr, features)
#################
############Assemble and save Figure Supplementary 11############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Supp_PDR1, interpolate=TRUE)

Figure11_Supp<-plot_grid(gpp)

ggsave (plot = Figure11_Supp, filename = "Supplementary_Fig11_low_quality.jpg", units = "cm", device = "jpg",width =18, height =3, dpi = 300, bg = "white")
ggsave (plot = Figure11_Supp, filename = "Supplementary_Fig11.png", units = "cm", device = "png",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Figure11_Supp, filename = "Supplementary_Fig11.jpg", units = "cm", device = "jpg",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Figure11_Supp, filename = "Supplementary_Fig11.svg", units = "cm", device = "svg",width =18, height =3, dpi = 1000, bg = "white")
ggsave (plot = Figure11_Supp, filename = "Supplementary_Fig11.pdf", units = "cm", device = "pdf",width =18, height =3, dpi = 1000, bg = "white")
#################

############Figure Supplementary 12############
#Create tree
tree <- read.tree(text = "(Nakaseomyces\nglabratus:0.02,(Saccharomyces\nbayanus:0.2,(Saccharomyces\nkudriavzevii:0.1,(Saccharomyces\nmikatae:0.2,(Saccharomyces\ncerevisiae:0.1,Saccharomyces\nparadoxus:0.1))))C:10);")

ggtree_obj <- ggtree(tree)

a<- ggtree_obj + geom_tiplab(fontface=4) +hexpand(0.5,direction=1)
#################
############Assemble and save Figure Supplementary 12############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Supp_candida, width = unit(1, "npc"), height = unit(1, "npc"))

Figure12_Supp<-plot_grid(gpp)

ggsave (plot = Figure12_Supp, filename = "Supplementary_Fig12_low_quality.jpg", units = "cm", device = "jpg",width =45, height =45, dpi = 300, bg = "white")
ggsave (plot = Figure12_Supp, filename = "Supplementary_Fig12.png", units = "cm", device = "png",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Figure12_Supp, filename = "Supplementary_Fig12.jpg", units = "cm", device = "jpg",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Figure12_Supp, filename = "Supplementary_Fig12.svg", units = "cm", device = "svg",width =45, height =45, dpi = 1000, bg = "white")
ggsave (plot = Figure12_Supp, filename = "Supplementary_Fig12.pdf", units = "cm", device = "pdf",width =45, height =45, dpi = 1000, bg = "white")
#################

############Figure Supplementary 13############
##Panel A
Supp_new_genome3c_new2 <-  dplyr:::filter(Supp_genome, Strain=="3Hybrid")
Supp_new_genome3c_new2 <-  dplyr:::filter(Supp_genome, Type=="Evolved_NQO")

exampleHybrid <- dplyr:::filter(Supp_genome, Strain=="3Hybrid")
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=1)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=21)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=10)
exampleHybrid <-dplyr:::filter(exampleHybrid,Replicate!=25)
exampleHybrid<-filter(exampleHybrid,Type=="Evolved_NQO")

nq_snps_HYBRID <-  dplyr:::filter(exampleHybrid, Strain=="3Hybrid")
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read >1.5, 1.5,Normalized_read))
nq_snps_HYBRID<- nq_snps_HYBRID %>% mutate(Normalized_read=if_else(Normalized_read < -1.5, -1.5,Normalized_read))

nq_snps_HYBRID<- filter(nq_snps_HYBRID,Replicate==28)
nq_genome<- nq_snps_HYBRID

legend_title <- "Relative Read Depth"

FigSupp13A <-nq_genome %>% dplyr:::filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x = as.numeric(start2),y = as.factor(chrom4),fill=Normalized_read)) + 
  theme_prism() + 
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_tile(aes(fill =Normalized_read),height=0.65,size=0.2) +
  scale_fill_gradient2(legend_title,low = "darkslateblue", mid = "lavender", high = "orange3", midpoint = 0, limits = c(-2.9, 2.9),breaks= c(-2,-1,0,1,2))+
  scale_y_discrete(limits=c("32","31","30","29","28","27","26","25","24","23","22","21","20","19","18","17",
                            "16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"),
                   labels=c("Spar XVI","Scer chrXVI","Spar chrXV","Scer chrXV",
                            "Spar chrXIV","Scer chrXIV","Spar chrXIII","Scer chrXIII",
                            "Spar chrXII","Scer chrXII","Spar chrXI","Scer chrXI",
                            "Spar chrX","Scer chrX","Spar chrIX","Scer chrIX",
                            "Spar chrVIII","Scer chrVIII","Spar chrVII","Scer chrVII",
                            "Spar chrVI","Scer chrVI","Spar chrV","Scer chrV",
                            "Spar chrIV","Scer chrIV","Spar chrIII","Scer chrIII",
                            "Spar chrII","Scer chrII","Spar chrI","Scer chrI"))+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),
                     labels=c(0,500,1000,1500))+
  facet_wrap(.~ Replicate, ncol=13)+
  xlab("Coordinate (kbp)")+
  ylab("Chromosome") +
  theme(strip.text = element_text(
    size = 20))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  theme(legend.text = element_text(size = 30))+
  theme(legend.background = element_rect(fill="snow2", 
                                         size=0.4, linetype="solid",colour ="grey"))+
  theme(legend.title = element_text(angle = 270, vjust= 0.2,hjust= 0.9),
        legend.title.align=0.1)+
  theme(legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))+
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text =element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust=1))

#Panel B
aneuploidy<-Supp_raw_cov
aneuploidy$Replicate<-as.numeric(aneuploidy$Replicate)

all_mean_genome_chrom<-aneuploidy %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon","Spar chrI",
                                                              ifelse(chrom=="utg1271_pilon","Spar chrII",
                                                                     ifelse(chrom=="utg584_pilon","Spar chrIII",
                                                                            ifelse(chrom=="utg988_pilon","Spar chrIV",
                                                                                   ifelse(chrom=="utg69_pilon","Spar chrV",
                                                                                          ifelse(chrom=="utg639_pilon","Spar chrVI",
                                                                                                 ifelse(chrom=="utg199_pilon","Spar chrVII",
                                                                                                        ifelse(chrom=="utg176_pilon","Spar chrVIII",
                                                                                                               ifelse(chrom=="utg1121_pilon","Spar chrIX",
                                                                                                                      ifelse(chrom=="utg245_pilon","Spar chrX",
                                                                                                                             ifelse(chrom=="utg298_pilon","Spar chrXI",
                                                                                                                                    ifelse(chrom=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                           ifelse(chrom=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                                  ifelse(chrom=="utg210_pilon","Spar chrXIII",
                                                                                                                                                         ifelse(chrom=="utg48_pilon","Spar chrXIV",
                                                                                                                                                                ifelse(chrom=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                       ifelse(chrom=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                              ifelse(chrom=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                     ifelse(chrom=="chrI","Scer chrI",
                                                                                                                                                                                            ifelse(chrom=="chrII","Scer chrII",
                                                                                                                                                                                                   ifelse(chrom=="chrIII","Scer chrIII",
                                                                                                                                                                                                          ifelse(chrom=="chrIV","Scer chrIV",
                                                                                                                                                                                                                 ifelse(chrom=="chrV","Scer chrV",
                                                                                                                                                                                                                        ifelse(chrom=="chrVI","Scer chrVI",
                                                                                                                                                                                                                               ifelse(chrom=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                      ifelse(chrom=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                             ifelse(chrom=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                    ifelse(chrom=="chrX","Scer chrX",
                                                                                                                                                                                                                                                           ifelse(chrom=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                                  ifelse(chrom=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                         ifelse(chrom=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                                ifelse(chrom=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                       ifelse(chrom=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                              ifelse(chrom=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))
all_pools_samples<- all_mean_genome_chrom 
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples$gainloss <-ifelse(all_pools_samples$log2>7,"Gain",NA)
all_pools_samples$gainloss <-ifelse(all_pools_samples$log2<2,"Loss",all_pools_samples$gainloss)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))
all_pools_samples1 <- filter(all_pools_samples, Unique_value=="Yes")
all_pools_samples1$Replicate <- as.numeric(all_pools_samples1$Replicate) 

NQ_hyb<-all_pools_samples%>% filter(Type=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(Strain=="3Hybrid")
NQ_hyb28 <- NQ_hyb %>% filter(Replicate == 28)

hybnq <-NQ_hyb28

FigSupp13B <-hybnq %>%
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read \n depth (log2)")+
  ylim(0,10)+#+scale_x_discrete(breaks=c("a","b"),labels=c("2","3"))
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Panel C
ploidy <- Supp_ploidy
x<-ploidy%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))

fdata <- x
newdat53n <- filter(newdat, LMH>184)

triploides <- newdat53n %>% group_by(evolved,same_Replicate,specie) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[5]  <- "counts"
control3n<-filter(x, evolved=="3n")
control2n<-filter(x, evolved=="2n")
controles <- filter(x, evolved=="2n" | evolved=="3n")

fdata$log_FL1A = log1p(fdata$FL1_A)
control2n$log_FL1A = log1p(control2n$FL1_A)
control3n$log_FL1A = log1p(control3n$FL1_A)

legend_title <- " "

NQ_hyb<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(specie=="3Hybrid")
NQ_hyb28 <- NQ_hyb %>% filter(replicate == 28)

tri28graph <-rbind(NQ_hyb28,control3n,control2n)

FigSupp13C<- tri28graph %>% ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+
  ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "honeydew4", "#FF9999"),
                    labels = (c("2n", "3n", "Replicate 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(panel.spacing = unit(1, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#################
############Assemble and save Supplementary Figure 13############
FigA_label<- plot_grid(FigSupp13A,labels="a",label_size=25)
FigB_label<- plot_grid(FigSupp13B,labels="b",label_size=25)
FigC_label<- plot_grid(FigSupp13C,labels="c",label_size=25)

Figsupp13_left <- plot_grid(FigB_label,FigC_label,nrow=2)

Figsupp13_all <- plot_grid(FigA_label,Figsupp13_left,nrow=1, rel_widths = c(1,1.7))

ggsave (plot = Figsupp13_all, filename = "Supplementary_Fig13_low_quality.jpg", units = "cm", device = "jpg",width = 40, height =20, dpi = 300,bg = "white")
ggsave (plot = Figsupp13_all, filename = "Supplementary_Fig13.png", units = "cm", device = "png",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp13_all, filename = "Supplementary_Fig13.jpg", units = "cm", device = "jpg",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp13_all, filename = "Supplementary_Fig13.svg", units = "cm", device = "svg",width = 40, height =20, dpi = 1000,bg = "white")
ggsave (plot = Figsupp13_all, filename = "Supplementary_Fig13.pdf", units = "cm", device = "pdf",width = 40, height =20, dpi = 1000,bg = "white")
#################

############Figure Supplementary 14############
#Supp 14A
Supp_14A<- Supp_growth%>% dplyr::filter(Specie=="BY4741")%>%
  ggplot(aes(x = hour_Rounded, y = mean_od, ymin = mean_od - error_estandar, ymax = mean_od + error_estandar, fill = Mutation)) +
  geom_ribbon(alpha = 0.5) +
  scale_shape_identity() +  
  geom_line(aes(col=as.factor(Mutation),alpha=0.4)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  theme_bw()+
  xlim(0,40)+
  theme(legend.position = "none",
        axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_color_manual(values = c("indianred1","tan1","gold4","mediumaquamarine",
                                "violet","steelblue1"),
                     breaks= c("Empty", "Wild type (WT)", "G280R","G280S","M308I","G1041/2W"),
                     labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  scale_fill_manual(values = c("indianred1","tan1","gold4","mediumaquamarine",
                               "violet","steelblue1"),
                    breaks= c("Empty", "Wild type (WT)", "G280R","G280S","M308I","G1041/2W"),
                    labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))
##Supp 14B
#Statistics
Supp_expression$Mutation <- as.factor(Supp_expression$Mutation)
amod <- aov(`GRN-B-HLog`~Mutation, data=Extended_expression)
summary(amod)
inter.test1 <- glht(amod,  mcp(Mutation = "Tukey"))
summary(inter.test1)
cld(inter.test1)

orden_desired <- c("Empty", "WT", "G280R","G280S","M308I","G1041/2W")
Supp_expression$Mutation <- factor(Supp_expression$Mutation, levels = orden_desired)

legend_title <- " "

Supp_14B<-Supp_expression %>% ggplot(aes(x = Mutation,y=`GRN-B-HLog`,fill=Mutation),  group=well) +
  geom_violin(trim=FALSE,alpha=0.7) + xlab("Mutation")+ylab("Pdr5-GFP (U.A. Fluorescence)")+ 
  theme(panel.spacing = unit(0.72, "cm"))+
  scale_color_manual(legend_title,values = c("indianred1","tan1","gold4","mediumaquamarine",
                                             "violet","steelblue1"),
                     breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                     labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  scale_fill_manual(legend_title,values = c("indianred1","tan1","gold4","mediumaquamarine",
                                            "violet","steelblue1"),
                    breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                    labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  guides(fill = guide_legend(nrow = 1))+
  annotate("text",
           y = c(4),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=5) +
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.position = "top")

leg <- get_legend(Supp_14B)

Supp_14B<- Supp_expression %>% ggplot(aes(x = Mutation,y=`GRN-B-HLog`,fill=Mutation),  group=well) +
  geom_violin(trim=FALSE,alpha=0.7) + 
  xlab("Amino acid change")+ylab("Pdr5-GFP \n (Arb. units fluorescence)")+ 
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="none")+
  scale_fill_manual(values = c("indianred1","tan1",
                               "gold4",
                               "mediumaquamarine",
                               "violet","steelblue1"))+
  scale_x_discrete(breaks= c("Empty", "WT", "G280R","G280S","M308I","G1041/2W"),
                   labels = c("Empty", "WT", "G280R","G280S","M308I","G1042W"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank())+
  annotate("text",
           y = c(4),
           x = c(1.3),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=4) 

#################
############Assemble and save Figure Supplementary 14############
FigA_label<- plot_grid(Supp_14A,labels="a",label_size=15)
FigB_label<- plot_grid(Supp_14B,labels="b",label_size=15)

Fig14AB<- plot_grid(FigA_label,FigB_label,nrow=1,rel_widths = c(1,1.05))

Fig14ABleg<- plot_grid(leg,Fig14AB,nrow=2,rel_heights = c(0.2,1))

ggsave (plot = Fig14ABleg, filename = "Supplementary_Fig14_low_quality.jpg", units = "cm", device = "jpg",width =22, height =12, dpi = 300, bg = "white")
ggsave (plot = Fig14ABleg, filename = "Supplementary_Fig14.png", units = "cm", device = "png",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig14ABleg, filename = "Supplementary_Fig14.jpg", units = "cm", device = "jpg",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig14ABleg, filename = "Supplementary_Fig14.svg", units = "cm", device = "svg",width =22, height =12, dpi = 1000, bg = "white")
ggsave (plot = Fig14ABleg, filename = "Supplementary_Fig14.pdf", units = "cm", device = "pdf",width =22, height =12, dpi = 1000, bg = "white")
#################

############Figure Supplementary 15############
a1<- filter(Supp_growth_plasmids, Type=="WT + pmoby")
a2<- filter (Supp_growth_plasmids, Type=="pdr1Δ + pmoby")

all_data2<-rbind(a1,a2)

all_data2_day1<-all_data2%>% filter(Day==1)

fdata<-Supp_growth_curves_plasmids

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr2<- mutate(fdata1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2,Day==2)

fdataapdr2$hour<-as.numeric(fdataapdr2$hour) 
fdataapdr2$hour_Rounded <- round(fdataapdr2$hour, digits = 1)

data_new3<- fdataapdr2%>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(mean_od = mean(od))

error_estandar <- fdataapdr2 %>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(error_estandar = sd(od) / sqrt(n()))

data_new4<-  full_join(data_new3,error_estandar,by=c("Type"="Type","Specie"="Specie","Mutation"="Mutation","hour_Rounded"="hour_Rounded","Condition.y"="Condition.y","Type2"="Type2"))

data_new4<- filter(data_new4,Mutation!="G282V")

data_new4<- data_new4 %<>% mutate(Mutation = ifelse(Mutation=="Control","WT",Mutation))

data_new4<- data_new4 %<>% mutate(Specie = ifelse(Specie=="BY4741 a","BY4741",Specie))

data_new4<- data_new4 %<>% filter(Specie!="BY4741 alpha")

data_new4<- data_new4 %<>% mutate(Type = ifelse(Type=="WT + pmoby","WT + pmoby_scer",Type))

data_new4<- data_new4 %<>% mutate(Type = ifelse(Type=="Δpdr1 + pmoby","Δpdr1 + pmoby_scer",Type))

new_names <- c("BY4741", "Scer", "Spar")
levels(data_new4$Specie) <- c("BY4741", "italic('S.cerevisiae')", "italic('S.paradoxus')")
levels(data_new4$Specie) <- c("italic('p=')*0.001", "italic('p=')*0.01", "italic('p=')*0.05")

data_bold<- data_new4
data_bold$Specie <- factor(data_bold$Specie,        # Change factor labels
                           labels = c("BY4741",
                                      "italic(S.cerevisiae)",
                                      "italic(S.paradoxus)"))

#Make changes for figures purposes
data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmoby", "pdr1Δ + pmoby", data_bold$Type)
data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmobyEMPTY", "pdr1Δ + pmobyEMPTY", data_bold$Type)
data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmoby_scer", "pdr1Δ + pmoby_scer", data_bold$Type)

legend_title <- " "
#Scer growth curves 
Figure7_v2<- data_bold%>% filter(Condition.y=="NQO_4μM") %>% 
  ggplot(aes(y = mean_od,x= hour_Rounded,ymin = mean_od - error_estandar, ymax = mean_od + error_estandar,group=Mutation)) +
  geom_ribbon(alpha = 0.5, aes(fill=Mutation)) +
  geom_line(aes(col=Mutation)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  theme_bw()+
  ylim(0.1,1.1)+
  facet_grid(Type~Specie,labeller=label_parsed)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 19),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white") )+
  scale_linetype_manual(values=c("dashed","solid"))+
  theme(legend.position="top",
        legend.direction = "horizontal")+
  scale_color_manual(legend_title,values = c("steelblue1",
                                             "gold4","mediumaquamarine",
                                             "violet","tan1"),name = " ",
                     breaks = c("G1041/2W","G280R","G280S","M308I","WT"),
                     labels = c("G1042W","G280R","G280S","M308I","WT"))+
  scale_fill_manual(legend_title,values = c("steelblue1",
                                            "gold4","mediumaquamarine",
                                            "violet","tan1"),name = " ",
                    breaks = c("G1041/2W","G280R","G280S","M308I","WT"),
                    labels = c("G1042W","G280R","G280S","M308I","WT"))+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14))

all_data<-Supp_growth_plasmids

a1<- filter(all_data, Type=="WT + pmoby")
a2<- filter (all_data, Type=="Δpdr1 + pmoby")

all_data2<-rbind(a1,a2)

all_data2<-all_data2 %>% filter(Mutation!="G282V")

all_data2_day2<-all_data2%>% filter(Day==2)

all_data2_day2<- all_data2_day2 %<>% mutate(Type = ifelse(Type=="WT + pmoby","WT + pmoby_scer",Type))

all_data2_day2<- all_data2_day2 %<>% mutate(Type = ifelse(Type=="Δpdr1 + pmoby","Δpdr1 + pmoby_scer",Type))

all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby", "pdr1Δ + pmoby", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmobyEMPTY", "pdr1Δ + pmobyEMPTY", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby_scer", "pdr1Δ + pmoby_scer", all_data2_day2$Type)

check1_BY_day2<- all_data2_day2 %>% filter(Specie=="BY4741") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 3,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

check1_Scer_day2<- all_data2_day2 %>% filter(Specie=="Scer") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

check1_Spar_day2<- all_data2_day2 %>% filter(Specie=="Spar") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("M308I","G280S","G280R","G1041/2W","Control"),
                   labels = c("M308I","G280S","G280R","G1042W","WT")) 

#Scer heatmap
combined_plot <- plot_grid(check1_BY_day2,check1_Scer_day2,check1_Spar_day2,ncol=3)

fdata<- Supp_growth_curves_plasmids2

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr1 <- dplyr::filter(fdata1, !grepl('Δpdr1', Type))

fdataapdr2<- mutate(fdataapdr1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))

all_data<-Supp_growth_plasmids2

a1<- filter(all_data, Type=="WT + pmoby_spar")
a2<- filter (all_data, Type=="Δpdr1 + pmoby_spar")

all_data2<-rbind(a1,a2)

all_data2_day2<-all_data2%>% filter(Day==2)

all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby", "pdr1Δ + pmoby", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmobyEMPTY", "pdr1Δ + pmobyEMPTY", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby_spar", "pdr1Δ + pmoby_spar", all_data2_day2$Type)

check1_BY_day2<- all_data2_day2 %>% filter(Specie=="BY4741") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 3,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                   labels = c("G281V", "G279S", "G279R", "WT"))

all_data2_day2<-all_data2%>% filter(Day==2)

all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby", "pdr1Δ + pmoby", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmobyEMPTY", "pdr1Δ + pmobyEMPTY", all_data2_day2$Type)
all_data2_day2$Type <- ifelse(all_data2_day2$Type == "Δpdr1 + pmoby_spar", "pdr1Δ + pmoby_spar", all_data2_day2$Type)

check1_Scer_day2<- all_data2_day2 %>% filter(Specie=="Scer") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                   labels = c("G281V", "G279S", "G279R", "WT"))

check1_Spar_day2<- all_data2_day2 %>% filter(Specie=="Spar") %>% ggplot(aes(x = interaction(Specie,Type), y = Mutation, fill = aucexp)) +
  geom_tile() + 
  facet_grid(.~factor(Condition.y, levels=c('NQO_4μM','NQO_8μM','NQO_10μM','Control')))+
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "brown3", mid = "burlywood1", high = "green4", midpoint = 5,limits= c(0,36),breaks= c(1,9,18,27,36))+
  theme_bw()+
  ylab("Amino acid \n change")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank())+
  scale_y_discrete(breaks = c("G282V", "G280S", "G280R", "Control"),
                   labels = c("G281V", "G279S", "G279R", "WT"))

#Spar heatmap
combined_plot2 <- plot_grid(check1_BY_day2,check1_Scer_day2,check1_Spar_day2,ncol=3)

fdata<-Supp_growth_curves_plasmids2

fdata1 <- filter(fdata,Type!="WT")
fdata1 <- filter(fdata1,Type!="Δpdr1")

fdataapdr2<- mutate(fdata1, Type2=ifelse(grepl('^Control',Mutation), 'Control', "Mutation"))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('EMPTY', Type))
fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('CRISPR', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2, !grepl('pmoby_scer', Type))

fdataapdr2 <- dplyr::filter(fdataapdr2,Day==2)

fdataapdr2$hour<-as.numeric(fdataapdr2$hour) 
fdataapdr2$hour_Rounded <- round(fdataapdr2$hour, digits = 1)

data_new3<- fdataapdr2%>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(mean_od = mean(od))

error_estandar <- fdataapdr2 %>%
  dplyr::group_by(Type,Specie,Mutation,hour_Rounded,Condition.y,Type2) %>%
  dplyr::summarise(error_estandar = sd(od) / sqrt(n()))

data_new4<-  full_join(data_new3,error_estandar,by=c("Type"="Type","Specie"="Specie","Mutation"="Mutation","hour_Rounded"="hour_Rounded","Condition.y"="Condition.y","Type2"="Type2"))

data_new4<- data_new4 %<>% mutate(Mutation = ifelse(Mutation=="Control","WT",Mutation))

data_new4 <- filter(data_new4,Specie!="White NQO")
data_new4 <- filter(data_new4,Specie!="White control")
data_bold<- data_new4
data_bold$Specie <- factor(data_bold$Specie,        
                           labels = c("BY4741",
                                      "italic(S.cerevisiae)",
                                      "italic(S.paradoxus)"))

data_bold<- filter(data_bold,hour_Rounded <41)

data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmoby", "pdr1Δ + pmoby", data_bold$Type)
data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmobyEMPTY", "pdr1Δ + pmobyEMPTY", data_bold$Type)
data_bold$Type <- ifelse(data_bold$Type == "Δpdr1 + pmoby_spar", "pdr1Δ + pmoby_spar", data_bold$Type)

#Spar growth curves 
Figure7b_v2<- data_bold %>% filter(Condition.y=="NQO_4μM") %>% 
  ggplot(aes(y = mean_od,x= hour_Rounded,
             ymin = mean_od - error_estandar,
             ymax = mean_od + error_estandar,group=Mutation)) +
  geom_ribbon(alpha = 0.5, aes(fill=Mutation)) +
  geom_line(aes(col=Mutation)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("OD (595 nm)")+xlab("Time (hours)")+
  theme(axis.text.x = element_text()) +
  ylim(0.1,1.1)+
  theme_bw()+
  facet_grid(Type~Specie,labeller=label_parsed)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 19),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white") )+
  theme(legend.position="top",
        legend.direction = "horizontal")+
  scale_color_manual(legend_title,breaks = c("G280R", "G280S","G282V", "WT"),
                     labels = c("G279R", "G279S","G281V",  "WT"),
                     values = c("gold4","mediumaquamarine","purple","tan1",
                                "gold4","mediumaquamarine","purple","tan1"),name = " ")+
  scale_fill_manual(legend_title,breaks = c("G280R", "G280S","G282V", "WT"),
                    labels = c("G279R", "G279S","G281V",  "WT"),
                    values = c("gold4","mediumaquamarine","purple","tan1",
                               "gold4","mediumaquamarine","purple","tan1"),name = " ")+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14))
#################
############Assemble and save Supplementary Figure 15############
Fig_supp<- plot_grid(Figure7_v2,Figure7b_v2,nrow=1)
Fig_suppB <- plot_grid(combined_plot,combined_plot2,nrow=2)
Fig_supp_A<- plot_grid(Fig_supp,labels="a",label_size=26)
Fig_suppB_B<- plot_grid(Fig_suppB,labels="b",label_size=26)
FigSupp15 <- plot_grid(Fig_supp_A,Fig_suppB_B,nrow=2)

ggsave (plot = FigSupp15, filename = "Supplementary_Fig15_low_quality.jpg", units = "cm", device = "jpg",width = 40, height =40, dpi = 300,bg = "white")
ggsave (plot = FigSupp15, filename = "Supplementary_Fig15.png", units = "cm", device = "png",width = 40, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp15, filename = "Supplementary_Fig15.jpg", units = "cm", device = "jpg",width = 40, height =40, dpi = 1000,bg = "white")
ggsave (plot = FigSupp15, filename = "Supplementary_Fig15.svg", units = "cm", device = "svg",width = 40, height =40, dpi = 700,bg = "white")
ggsave (plot = FigSupp15, filename = "Supplementary_Fig15.pdf", units = "cm", device = "pdf",width = 40, height =40, dpi = 700,bg = "white")
#################

############Figure Supplementary 16############
#Figure A
all_arranged2<- Supp16_boxplot

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
Figure1E_BY <- all_arranged_BY %>% filter(Condition.y=="NQO_4μM")%>% 
  ggplot(aes(x=interaction(Type2,Specie), y=rval,group=interaction(Specie,Type2))) +
  geom_point(colour="black",pch=21,height = 0, size=2, aes(fill=as.factor(Type2)))+
  geom_boxplot(outlier.shape = NA,aes(col=Mutation,fill=Mutation,alpha=0.7, group=interaction(Specie,Type2)))+
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "Ancestral_homozygous" = "gray30", "Haploid"="antiquewhite1"))+
  scale_fill_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "Ancestral_homozygous" = "gray30", "Haploid"="antiquewhite1"))+
  ylab("Growth rate \n (OD/hour)") +
  xlab("")+
  theme_bw() +
  facet_grid(.~Specie)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=12,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  ylim(0,0.7)+
  theme(axis.text.x = element_text(angle = 90))

#make pvalue
all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="Hybrid")

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="Spar")

all_arranged_BY$Type2 <- as.factor(all_arranged_BY$Type2)
amod <- aov(rval~Type2, data=all_arranged_BY)
summary(amod)
inter.test1 <- glht(amod,  mcp(Type2= "Tukey"))
summary(amod)
summary(inter.test1)
cld(inter.test1)

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

legend_title<-""
Figure1E_BY <- all_arranged_BY %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(legend_title,
                    values = c("gray30","mediumpurple3","#21908CFF"),
                    limits=c("Heterozygous",
                             "Homozygous",
                             "Ancestral_homozygous"),
                    breaks= c("Heterozygous",
                              "Homozygous",
                              "Ancestral_homozygous"),
                    labels= c("No mutation",
                              "Heterozygous",
                              "Homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=5)+
  theme(legend.direction="horizontal",
        legend.key.width = unit(8, 'mm'),
        legend.key.height = unit(8, 'mm'),
        legend.text = element_text(size = 12))

leg <- get_legend(Figure1E_BY)

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~cerevisiae) ~ "- BY4743")
  }))
}

Figure1E_BY <- all_arranged_BY %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.05"),
           family = "", fontface = 3, size=5)

all_arranged_BY2<- all_arranged2 %>% filter(Specie=="Hybrid")

outcome_labels <- c("No mutation", 
                    "Heterozygous /n Spar", 
                    "Heterozygous /n Scer",
                    "Homozygous")


Figure1E_BY2 <- all_arranged_BY2 %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3", "mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous_BY-M308I",
             "Heterozygous_Spar-M308I",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous_BY-M308I",
              "Heterozygous_Spar-M308I",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous_BY-M308I",
                            "Heterozygous_Spar-M308I",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous_BY-M308I",
                             "Heterozygous_Spar-M308I",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous \n Scer-M308I",
                             "Heterozygous \n Spar-M307I",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.0001"),
           family = "", fontface = 3, size=5) 

all_arranged_BY3<- all_arranged2 %>% filter(Specie=="Spar")

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~paradoxus) ~ " ")
  }))
}

all_arranged_BY3$Specie<-"S.paradoxus" 
Figure1E_BY3 <- all_arranged_BY3 %>% filter(Condition.y == "NQO_4μM") %>%
  ggplot(aes(x = Type2, y = rval)) +
  geom_boxplot(colour = "black", outlier.shape = NA, aes(colour = "black", fill = Type2, alpha = 0.7)) +
  geom_point(colour = "black", pch = 21, size = 2, aes(fill = Type2)) +
  scale_fill_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  scale_x_discrete(limits=c("Ancestral_homozygous",
                            "Heterozygous",
                            "Homozygous"),
                   breaks= c("Ancestral_homozygous",
                             "Heterozygous",
                             "Homozygous"),
                   labels= c("No \n mutation",
                             "Heterozygous",
                             "Homozygous"))+
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylim(0,0.65)+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  annotate("text",
           y = c(0.64),
           x = c(1),
           label = c("p < 0.001"),
           family = "", fontface = 3, size=5) 

#Arrange a
Figure1E_BY4<-plot_grid(Figure1E_BY,Figure1E_BY2,Figure1E_BY3,nrow=1,rel_widths = c(1.2,1.3,1))
FigureSuppAleg<-plot_grid(leg, Figure1E_BY4,nrow=2,rel_heights = c(1,7))

#Figure B
all_arranged2 <- Supp16_curves
all_arranged2<- all_arranged2%>% filter(Condition.y!="Control")
all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")

Figure1E_BY <- all_arranged_BY %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),    
        strip.text.y = element_text(size = 12),      
        strip.background = element_rect(fill = "white"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

all_arranged_BY<- all_arranged2 %>% filter(Specie=="BY")
all_arranged_BY$Specie<-"S.cerevisiae - BY4743" 

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~cerevisiae) ~ "- BY4743")
  }))
}

Figure1E_BY <- all_arranged_BY %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  facet_grid(. ~ Specie, labeller = labels_custom)+
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

all_arranged_BY2<- all_arranged2 %>% filter(Specie=="Hybrid")

##add data for legend purposes
data_leg <- data.frame(
  Well = c("...1"),
  Day = c(1),
  Plate = c("Plate1"),
  Condition.x = c("Condition1"),
  time = c(0),
  od = c(NA),
  Specie = c("Hybrid"),
  Type = c("Type1"),
  Line = c("Line1"),
  Replicate = c(1),
  Mutation = c("Heterozygous"),
  Condition.y = c("Condition1"),
  Assay = c("Assay1"),
  hour = c(0),
  Type2 = c("Heterozygous")
)

all_arranged_BY2B <-all_arranged_BY2[, -c(1, 2)]

all_arranged_BY2B<-rbind(data_leg,all_arranged_BY2B)

legend_title <- ""
Figure1E_BY2 <- all_arranged_BY2B %>% dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2)) +
  scale_color_manual(legend_title,
                     values = c("mediumpurple3",
                                "darkblue",
                                "magenta3",
                                "#21908CFF",
                                "gray30"),
                     limits=c("Heterozygous",
                              "Heterozygous_BY-M308I",
                              "Heterozygous_Spar-M308I",
                              "Homozygous",
                              "Ancestral_homozygous"),
                     breaks= c("Heterozygous",
                               "Heterozygous_BY-M308I",
                               "Heterozygous_Spar-M308I",
                               "Homozygous",
                               "Ancestral_homozygous"),
                     labels= c("Heterozygous",
                               "Heterozygous Scer-M308I",
                               "Heterozygous Spar-M307I",
                               "Homozygous",
                               "No mutation")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.direction =  "horizontal")+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.y = element_text(size=14, face = "bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14),
        legend.key.width = unit(8, 'mm'),
        legend.key.height = unit(8, 'mm'),
        legend.text = element_text(size = 12))

leg <- get_legend(Figure1E_BY2)

Figure1E_BY2 <- all_arranged_BY2 %>% dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2, alpha=0.4)) +
  scale_color_manual(
    values = c("darkblue", "magenta3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous_BY-M308I",
             "Heterozygous_Spar-M308I",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous_BY-M308I",
              "Heterozygous_Spar-M308I",
              "Homozygous",
              "Ancestral_homozygous")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  facet_grid(. ~ Specie)+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size=20, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

outcome_labels <- c("No mutation", 
                    "Heterozygous /n Spar", 
                    "Heterozygous /n Scer",
                    "Homozygous")

all_arranged_BY3<- all_arranged2 %>% filter(Specie=="Spar")

all_arranged_BY3$Specie<-"S.paradoxus" 

labels_custom <- function(variable, value) {
  return(lapply(value, function(x) {
    bquote(italic(S.~paradoxus) ~ " ")
  }))
}

Figure1E_BY3 <- all_arranged_BY3 %>%  dplyr::filter(hour<21)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Type2,Type,Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Type2, alpha=0.4)) +
  facet_grid(. ~ Specie, labeller = labels_custom)+
  scale_color_manual(
    values = c("mediumpurple3",
               "#21908CFF","gray30"),
    limits=c("Heterozygous",
             "Homozygous",
             "Ancestral_homozygous"),
    breaks= c("Heterozygous",
              "Homozygous",
              "Ancestral_homozygous")) +
  ylab("Growth rate (OD/hour)") +theme_bw() +
  theme(legend.position = "none")+
  theme(axis.title = element_text(size=20, face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 19),
        axis.text.y =element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(size=14, face = "bold"),
        strip.text = element_text(size = 14))+
  theme(strip.text = element_text(face = "italic"))+
  ylab("OD (595 nm)")+xlab("Time (hours)")

#arrange b
Figure1E_BY4<-plot_grid(Figure1E_BY,Figure1E_BY2,Figure1E_BY3,nrow=1,rel_widths = c(1.2,1.3,1))
FigureSuppBleg<-plot_grid(leg, Figure1E_BY4,nrow=2,rel_heights = c(1,7))

#################
############Assemble and save Figure Supplementary 16############
###Save panels 
FigA_label<- plot_grid(FigureSuppAleg,labels="a",label_size=32)
FigB_label<- plot_grid(FigureSuppBleg,labels="b",label_size=32)

FigureSuppleg<-plot_grid(FigA_label,FigB_label,nrow=2)

Fig_Supp16<-FigureSuppleg

ggsave (plot = Fig_Supp16, filename = "Supplementary_Fig16_low_quality.jpg", units = "cm", device = "jpg",width = 45, height =25, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp16, filename = "Supplementary_Fig16.png", units = "cm", device = "png",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp16, filename = "Supplementary_Fig16.jpg", units = "cm", device = "jpg",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp16, filename = "Supplementary_Fig16.svg", units = "cm", device = "svg",width = 45, height =25, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp16, filename = "Supplementary_Fig16.pdf", units = "cm", device = "pdf",width = 45, height =25, dpi = 1000, bg = "white")
#################

############Figure Supplementary 17############
#Hybrid 13
Hybrid <- dplyr::filter(Supp_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==13)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Supp17_LOH_line13_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-Supp17_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==13)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure1<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 13", fontface='bold')
Figure1<-plot_grid(title,Figure1, ncol=1, rel_heights=c(0.1, 1))

#Hybrid 28
Hybrid <- dplyr::filter(Supp_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==28)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Supp17_LOH_line28_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day

loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-Supp17_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==28)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure2<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 28", fontface='bold')
Figure2<-plot_grid(title,Figure2, ncol=1, rel_heights=c(0.1, 1))

#Hybrid 30
Hybrid <- dplyr::filter(Supp_expevol, strain=="3Hybrid")
Hybrid_30 <- dplyr::filter(Hybrid, rep==30)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table <- Supp17_LOH_line30_final

loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- loh_table
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") 

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

loh_table<-loh_table %<>% mutate(LOH_other_positions= ifelse(LOH_1=="Heterozygous","No",
                                                             ifelse(LOH_1=="Homozygous","Yes",NA)))

data<-Supp17_LOH_parents_table_proportion_hybrid %>% dplyr::filter(Replicate==30)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure3<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("Hybrid Line 30", fontface='bold')
Figure3<-plot_grid(title,Figure3, ncol=1, rel_heights=c(0.1, 1))

#Scer 25
Hybrid <- dplyr::filter(Supp_expevol, strain=="1Scer")
Hybrid_30 <- dplyr::filter(Hybrid, rep==25)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")

Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer25_1<- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer25_1<-scer25_1 + scale_color_grey() + theme_classic() 

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Scer")
loh_tableB <- dplyr::filter(loh_tableB, Line==25)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer25_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)") + facet_grid(.~Genotype)

scer25_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("#21908CFF", "gray30")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==25)

scer25_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer25_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer25_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer25_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer25_4,labels="D",label_size=10)

Figure4<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. cerevisiae Line 25", fontface='bold')
Figure4<-plot_grid(title,Figure4, ncol=1, rel_heights=c(0.1, 1))

#Spar 22
Hybrid <- dplyr::filter(Supp_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==22)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
spar22_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

spar22_1<-spar22_1 + scale_color_grey() + theme_classic() 

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==22)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
spar22_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")

spar22_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==22)

spar22_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(spar22_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar22_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar22_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar22_4,labels="D",label_size=10)

Figure5<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 22", fontface='bold')
Figure5<-plot_grid(title,Figure5, ncol=1, rel_heights=c(0.1, 1))

#Scer 2
Hybrid <- dplyr::filter(Supp_expevol, strain=="1Scer")
Hybrid_30 <- dplyr::filter(Hybrid, rep==2)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
scer2_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

scer2_1<-scer2_1+ scale_color_grey() + theme_classic() 

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Scer")
loh_tableB <- dplyr::filter(loh_tableB, Line==2)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
scer2_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")


scer2_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==2)

scer2_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(scer2_1,labels="A",label_size=10)
Fig1B2<- plot_grid(scer2_2,labels="B",label_size=10)
Fig1C2<- plot_grid(scer2_3,labels="C",label_size=10)
Fig1D2<- plot_grid(scer2_4,labels="D",label_size=10)

Figure6<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. cerevisiae Line 2", fontface='bold')
Figure6<-plot_grid(title,Figure6, ncol=1, rel_heights=c(0.1, 1))

#Spar 24
Hybrid <- dplyr::filter(Supp_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==24)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
spar24_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

spar24_1<-spar24_1 + scale_color_grey() + theme_classic() #+ ggtitle ("Line 30") + theme(plot.title = element_text(hjust = 0.5, size=20))

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==24)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
spar24_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")

spar24_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","red")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==24)

spar24_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(spar24_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar24_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar24_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar24_4,labels="D",label_size=10)

Figure7<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 24", fontface='bold')
Figure7<-plot_grid(title,Figure7, ncol=1, rel_heights=c(0.1, 1))

#Spar 27
Hybrid <- dplyr::filter(Supp_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==27)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
spar27_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

spar27_1<-spar27_1 + scale_color_grey() + theme_classic()

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==27)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
spar27_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")

spar27_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==27)

spar27_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(spar27_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar27_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar27_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar27_4,labels="D",label_size=10)

Figure8<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 27", fontface='bold')
Figure8<-plot_grid(title,Figure8, ncol=1, rel_heights=c(0.1, 1))

#Spar 29
Hybrid <- dplyr::filter(Supp_expevol, strain=="2Spar")
Hybrid_30 <- dplyr::filter(Hybrid, rep==29)
Hybrid_30_NQO <- dplyr::filter(Hybrid_30, NQO=="Yes")
Hybrid_30_NQO_2<-Hybrid_30_NQO %<>% mutate(Generations = ifelse(day==3,15,
                                                                ifelse(day==4,20,
                                                                       ifelse(day==5,25,
                                                                              ifelse(day==6,30,
                                                                                     ifelse(day==7,35,
                                                                                            ifelse(day==8,40,
                                                                                                   ifelse(day==9,45,
                                                                                                          ifelse(day==10,50,
                                                                                                                 ifelse(day==12,55,
                                                                                                                        ifelse(day==13,65,
                                                                                                                               ifelse(day==14,70,
                                                                                                                                      ifelse(day==15,75,
                                                                                                                                             ifelse(day==16,80,
                                                                                                                                                    ifelse(day==17,85,
                                                                                                                                                           ifelse(day==18,90,
                                                                                                                                                                  ifelse(day==19,95,
                                                                                                                                                                         ifelse(day==20,100,
                                                                                                                                                                                ifelse(day==21,105,NA)))))))))))))))))))
spar29_1 <- ggplot(data=Hybrid_30_NQO, aes(x=as.numeric(day), y=rval)) +
  geom_line()+
  geom_point() + scale_x_continuous(breaks = as.numeric(Hybrid_30_NQO$day))+xlab("Time (Cycle)")+
  ylab("Growth rate (OD/hour)")

spar29_1<-spar29_1 + scale_color_grey() + theme_classic()

loh_table<-Supp17_LOH_table_parents
loh_table$day <-loh_table$Time
loh_table$Cycle <-loh_table$day
loh_tableB <- dplyr::filter(loh_table, Genotype=="Spar")
loh_tableB <- dplyr::filter(loh_tableB, Line==29)
loh_tableB<-loh_tableB %<>% mutate(Generations = ifelse(day==1,5,
                                                        ifelse(day==3,15,
                                                               ifelse(day==4,20,
                                                                      ifelse(day==5,25,
                                                                             ifelse(day==6,30,
                                                                                    ifelse(day==7,35,
                                                                                           ifelse(day==8,40,
                                                                                                  ifelse(day==9,45,
                                                                                                         ifelse(day==10,50,
                                                                                                                ifelse(day==12,55,
                                                                                                                       ifelse(day==13,65,
                                                                                                                              ifelse(day==14,70,
                                                                                                                                     ifelse(day==15,75,
                                                                                                                                            ifelse(day==16,80,
                                                                                                                                                   ifelse(day==17,85,
                                                                                                                                                          ifelse(day==18,90,
                                                                                                                                                                 ifelse(day==19,95,
                                                                                                                                                                        ifelse(day==20,100,
                                                                                                                                                                               ifelse(day==21,105,NA))))))))))))))))))))
spar29_2<- ggplot(loh_tableB, aes(fill=Codon, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() +
  ylab("Relative Frequency")+xlab("Time (Cycle)")

spar29_3<- ggplot(loh_tableB, aes(fill=Mutation, x=Cycle)) + 
  geom_bar(position="fill", stat="count") + scale_x_continuous(breaks = as.numeric(loh_table$Cycle)) +
  theme_classic() + 
  ylab("Relative Frequency")+xlab("Time (Cycle)")+
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

data<-dplyr::filter(Supp17_LOH_parents_table_proportion, Replicate==29)

spar29_4<- ggplot(data, aes(x=Time, y=percentage, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  scale_x_continuous(breaks = as.numeric(loh_table$Time))+theme_bw()+
  ylab("Relative Frequency")+xlab("Time (Cycles)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30")) 

Fig1A2<- plot_grid(spar29_1,labels="A",label_size=10)
Fig1B2<- plot_grid(spar29_2,labels="B",label_size=10)
Fig1C2<- plot_grid(spar29_3,labels="C",label_size=10)
Fig1D2<- plot_grid(spar29_4,labels="D",label_size=10)

Figure9<- plot_grid(Fig1A2,Fig1C2,Fig1B2,Fig1D2)
title <- ggdraw() + draw_label("S. paradoxus Line 29", fontface='bold')
Figure9<-plot_grid(title,Figure9, ncol=1, rel_heights=c(0.1, 1))
#################
############Assemble and save raw Figure Supplementary 17############
###Save panels 
Fig_Supp_17<- plot_grid(Figure1,Figure2,Figure3,Figure6,Figure4,Figure5,Figure7,Figure8,Figure9,nrow=3)

ggsave (plot = Fig_Supp_17, filename = "Supplementary_Fig17_raw_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =61, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp_17, filename = "Supplementary_Fig17_raw.png", units = "cm", device = "png",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17, filename = "Supplementary_Fig17_raw.jpg", units = "cm", device = "jpg",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17, filename = "Supplementary_Fig17_raw.svg", units = "cm", device = "svg",width = 70, height =61, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17, filename = "Supplementary_Fig17_raw.pdf", units = "cm", device = "pdf",width = 70, height =61, dpi = 1000, bg = "white")
#################
############Assemble and save Figure Supplementary 17############
#We added some aesthetic modifications so we import it
gpp <- rasterGrob(Supp17_Sanger)

Fig_Supp_17b<-plot_grid(gpp)

ggsave (plot = Fig_Supp_17b, filename = "Supplementary_Fig17_low_quality.jpg", units = "cm", device = "jpg",width =72, height =52, dpi = 300, bg = "white")
ggsave (plot = Fig_Supp_17b, filename = "Supplementary_Fig17.png", units = "cm", device = "png",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17b, filename = "Supplementary_Fig17.jpg", units = "cm", device = "jpg",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17b, filename = "Supplementary_Fig17.svg", units = "cm", device = "svg",width =72, height =52, dpi = 1000, bg = "white")
ggsave (plot = Fig_Supp_17b, filename = "Supplementary_Fig17.pdf", units = "cm", device = "pdf",width =72, height =52, dpi = 1000, bg = "white")
#################

############Figure Supplementary 18############
loh_table_day2<-Supp_loh_table_day2
loh_table_day2<-dplyr::filter(loh_table_day2, Mutation!="White")
loh_table_day2<-dplyr::filter(loh_table_day2, Mutation!="White NQO")

order_manual <- c("No mutation", "Heterozygous", "Homozygous")

fdata_LOH<-Supp_fdata
fdata_LOH$day <- fdata_LOH$Cycles
fdata_LOH<-fdata_LOH%<>% mutate(Generations = ifelse(day==1,5,
                                                     ifelse(day==3,15,
                                                            ifelse(day==4,20,
                                                                   ifelse(day==5,25,
                                                                          ifelse(day==6,30,
                                                                                 ifelse(day==7,35,
                                                                                        ifelse(day==8,40,
                                                                                               ifelse(day==9,45,
                                                                                                      ifelse(day==10,50,
                                                                                                             ifelse(day==12,55,
                                                                                                                    ifelse(day==13,65,
                                                                                                                           ifelse(day==14,70,
                                                                                                                                  ifelse(day==15,75,
                                                                                                                                         ifelse(day==16,80,
                                                                                                                                                ifelse(day==17,85,
                                                                                                                                                       ifelse(day==18,90,
                                                                                                                                                              ifelse(day==19,95,
                                                                                                                                                                     ifelse(day==20,100,
                                                                                                                                                                            ifelse(day==21,105,NA))))))))))))))))))))
fdata_LOH <- dplyr::filter(fdata_LOH,Day==2)
fdata_LOH <- dplyr::filter(fdata_LOH,Mutation!="White")
fdata_LOH <- na.omit(fdata_LOH)
fdata_LOH_NQO <- dplyr::filter(fdata_LOH,Condition.y=="NQO_4μM")

#Data everyone until 23h
Figure_Extended_A<- fdata_LOH_NQO %>% dplyr::filter(hour<23.5)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Replicate,Mutation,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  facet_grid(Line~Generations)+
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=8,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12), 
        strip.background = element_rect(fill = "white"))

Figure_Extended_B<- fdata_LOH_NQO %>% dplyr::filter(hour<23.5)%>%
  ggplot( aes(y=od,x=hour,group=interaction(Replicate,Mutation,Generations,Assay))) +
  scale_shape_identity() +  
  geom_line(aes(col=Mutation, alpha=0.4)) +
  ylab("OD \n (595 nm)")+xlab("Time (hours)")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=7)) +
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(strip.text.x = element_text(size = 12),    
        strip.text.y = element_text(size = 12),      
        strip.background = element_rect(fill = "white"))

loh_table_day2b<-dplyr::filter(loh_table_day2,Codon!="White")
Figure_Extended_C <- loh_table_day2b %>% dplyr::filter(Condition.y=="NQO_4μM")%>% 
  ggplot(aes(x=as.numeric(Generations), y=rval)) +
  geom_boxplot(outlier.shape = NA,aes(col=Mutation,fill=Mutation,alpha=0.7, group=interaction(Generations,Mutation)))+
  scale_x_continuous(breaks = as.numeric(Supp_loh_table_day2$Generations))+xlab("Time (Generations)")+
  scale_color_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                     labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  scale_fill_manual(values = c("Heterozygous" = "mediumpurple3", "Homozygous" = "#21908CFF", "No mutation" = "gray30"),
                    labels = c("Heterozygous" = "Heterozygous PDR1 mutation","No mutation" = "No PDR1 mutation", "Homozygous" = "Homozygous PDR1 mutation")) +
  ylab("Growth rate \n (OD/hour)") +theme_bw() +
  facet_grid(.~Line)+
  theme(axis.title = element_text(size=20, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  theme(legend.position="none")

#Disclaimer: data and graph are invented for legend purposes
FigureE_B<-Supp_loh_table_day2%>% dplyr::filter(Mutation!="White")%>%
  ggplot(aes(x=as.numeric(Cycles), y=rval, fill=Mutation)) + 
  geom_area(alpha=0.6 , size=1, colour="black",pos = "stack") +
  ylab("Relative Frequency")+xlab("Time (Generations)") +
  scale_fill_manual(values = c("mediumpurple3", "#21908CFF","gray30","gray30")) +
  theme(legend.position = "top")+
  theme(axis.title = element_text(size=18, face = "bold"),
        strip.text = element_text(color = "black", face = "bold",size = 19))+
  theme(axis.title = element_text(size=20, face = "bold"),
        axis.text.x = element_text(size=16,face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black", face = "bold"))+
  theme(strip.text.x = element_text(size = 19))+
  theme(strip.text = element_text(size = 14))+
  theme(legend.text = element_text(size = 14),  
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 14, face = "bold"))

leg <- get_legend(FigureE_B)

#################
############Assemble and save Figure Supplementary 18############
FigA_label<- plot_grid(Figure_Extended_A,labels="a",label_size=32)
FigB_label<- plot_grid(Figure_Extended_B,labels="b",label_size=32)
FigC_label<- plot_grid(Figure_Extended_C,labels="c",label_size=32)
Fig_Extended_Hyb_up <-plot_grid(FigA_label,FigC_label,nrow=2)
Fig_Extended_Hyb <-plot_grid(Fig_Extended_Hyb_up,FigB_label,nrow=1,rel_widths = c(1,0.5))  
Figure18leg<-plot_grid(leg, Fig_Extended_Hyb,nrow=2,rel_heights = c(0.6,7))

ggsave (plot = Figure18leg, filename = "Supplementary_Fig18_low_quality.jpg", units = "cm", device = "jpg",width = 43, height =25, dpi = 300, bg = "white")
ggsave (plot = Figure18leg, filename = "Supplementary_Fig18.png", units = "cm", device = "png",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure18leg, filename = "Supplementary_Fig18.jpg", units = "cm", device = "jpg",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure18leg, filename = "Supplementary_Fig18.svg", units = "cm", device = "svg",width = 43, height =25, dpi = 1000, bg = "white")
ggsave (plot = Figure18leg, filename = "Supplementary_Fig18.pdf", units = "cm", device = "pdf",width = 43, height =25, dpi = 1000, bg = "white")
#################

############Figure Supplementary 19############
Supp_raw_cov$Replicate<-as.numeric(Supp_raw_cov$Replicate)
all_mean_genome_chrom<-Supp_raw_cov%<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon","Spar chrI",
                                                               ifelse(chrom=="utg1271_pilon","Spar chrII",
                                                                      ifelse(chrom=="utg584_pilon","Spar chrIII",
                                                                             ifelse(chrom=="utg988_pilon","Spar chrIV",
                                                                                    ifelse(chrom=="utg69_pilon","Spar chrV",
                                                                                           ifelse(chrom=="utg639_pilon","Spar chrVI",
                                                                                                  ifelse(chrom=="utg199_pilon","Spar chrVII",
                                                                                                         ifelse(chrom=="utg176_pilon","Spar chrVIII",
                                                                                                                ifelse(chrom=="utg1121_pilon","Spar chrIX",
                                                                                                                       ifelse(chrom=="utg245_pilon","Spar chrX",
                                                                                                                              ifelse(chrom=="utg298_pilon","Spar chrXI",
                                                                                                                                     ifelse(chrom=="utg110_pilon","Spar chrXII 1/2",
                                                                                                                                            ifelse(chrom=="utg11_pilon","Spar chrXII 2/2",
                                                                                                                                                   ifelse(chrom=="utg210_pilon","Spar chrXIII",
                                                                                                                                                          ifelse(chrom=="utg48_pilon","Spar chrXIV",
                                                                                                                                                                 ifelse(chrom=="utg155_pilon","Spar chrXV 1/2",
                                                                                                                                                                        ifelse(chrom=="utg122_pilon","Spar chrXV 2/2",
                                                                                                                                                                               ifelse(chrom=="utg675_pilon","Spar chrXVI",
                                                                                                                                                                                      ifelse(chrom=="chrI","Scer chrI",
                                                                                                                                                                                             ifelse(chrom=="chrII","Scer chrII",
                                                                                                                                                                                                    ifelse(chrom=="chrIII","Scer chrIII",
                                                                                                                                                                                                           ifelse(chrom=="chrIV","Scer chrIV",
                                                                                                                                                                                                                  ifelse(chrom=="chrV","Scer chrV",
                                                                                                                                                                                                                         ifelse(chrom=="chrVI","Scer chrVI",
                                                                                                                                                                                                                                ifelse(chrom=="chrVII","Scer chrVII",
                                                                                                                                                                                                                                       ifelse(chrom=="chrVIII","Scer chrVIII",
                                                                                                                                                                                                                                              ifelse(chrom=="chrIX","Scer chrIX",
                                                                                                                                                                                                                                                     ifelse(chrom=="chrX","Scer chrX",
                                                                                                                                                                                                                                                            ifelse(chrom=="chrXI","Scer chrXI",
                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXII","Scer chrXII",
                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXIII","Scer chrXIII",
                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXIV","Scer chrXIV",
                                                                                                                                                                                                                                                                                        ifelse(chrom=="chrXV","Scer chrXV",
                                                                                                                                                                                                                                                                                               ifelse(chrom=="chrXVI","Scer chrXVI",NA)))))))))))))))))))))))))))))))))))


all_pools_samples<- all_mean_genome_chrom 
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))

#Hybrid - Ancestor
Anc_hyb<-all_pools_samples %>%filter(Type=="Ancestor")
Anc_hyb<-Anc_hyb %>% filter(Strain=="3Hybrid")
Anc_hyb23 <- Anc_hyb %>% filter(Replicate == 23) 
Anc_hyb24 <- Anc_hyb %>% filter(Replicate == 24) 
Anc_hyb27 <- Anc_hyb %>% filter(Replicate == 27)

hybanc <-rbind(Anc_hyb23,Anc_hyb24,Anc_hyb27)

hybanc_read <- hybanc  %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Hybrid - Evolved in Control
Con_hyb<-all_pools_samples %>% filter(Type=="Evolved_control")
Con_hyb<-Con_hyb %>% filter(Strain=="3Hybrid")
Con_hyb20 <- Con_hyb %>% filter(Replicate == 20)
Con_hyb23 <- Con_hyb %>% filter(Replicate == 23) 
Con_hyb26 <- Con_hyb %>% filter(Replicate == 26) 
Con_hyb27 <- Con_hyb %>% filter(Replicate == 27)
Con_hyb28 <- Con_hyb %>% filter(Replicate == 28)

hybcon <-rbind(Con_hyb20,Con_hyb23,Con_hyb26,Con_hyb27,Con_hyb28)

hybcon_read <- hybcon %>%
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Hybrid - Evolved in UV mimetic
NQ_hyb<-all_pools_samples%>% filter(Type=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(Strain=="3Hybrid")
NQ_hyb1 <- NQ_hyb %>% filter(Replicate == 1) 
NQ_hyb26 <- NQ_hyb %>% filter(Replicate == 26) 
NQ_hyb27 <- NQ_hyb %>% filter(Replicate == 27) 
NQ_hyb28 <- NQ_hyb %>% filter(Replicate == 28)

hybnq <-rbind(NQ_hyb1,NQ_hyb26,NQ_hyb27,NQ_hyb28)

hybnq_read <-hybnq %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(col="#FF9999",size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(Replicate~.)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

#Ploidy graphs
ploidy <- Supp_ploidy

x<-ploidy%>%filter(!(plate=="A"&evolved=="4n"))
x<-x%>%filter(!(plate=="C"&evolved=="4n"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"&other_info=="WT ploidy2"))
x<-x%>%filter(!(plate=="B"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="A"&evolved=="2n"&replicate=="Control"))
x<-x%>%filter(!(plate=="C"&evolved=="2n"&specie=="Spar haploid nat alfa"))
x<-x%>%filter(!(plate=="A"&evolved=="3n"))
x<-x%>%filter(!(plate=="C"&evolved=="3n"))
x<-x%>%filter(!(plate=="A"&evolved=="n"))
x<-x%>%filter(!(plate=="B"&evolved=="n"))
x<-x%>%filter(!(plate=="C"&evolved=="n"&specie=="Spar haploid hyg alfa"))

fdata <- x
newdat53n <- filter(newdat5, LMH>184)

triploides <- newdat53n %>% group_by(evolved,same_Replicate,specie) %>% dplyr:::summarise(n()) %>% ungroup()

colnames(newdat)[5]  <- "counts"
control3n<-filter(x, evolved=="3n")
control2n<-filter(x, evolved=="2n")
controles <- filter(x, evolved=="2n" | evolved=="3n")

fdata$log_FL1A = log1p(fdata$FL1_A)
control2n$log_FL1A = log1p(control2n$FL1_A)
control3n$log_FL1A = log1p(control3n$FL1_A)

legend_title <- "Sample"

Anc_scer<-fdata %>% filter(evolved=="Ancestor")
Anc_scer<-Anc_scer %>% filter(specie=="1Scer")
Anc_scer3 <- Anc_scer %>% filter(replicate == 3) 
Anc_scer5 <- Anc_scer %>% filter(replicate == 5) 
Anc_scer11 <- Anc_scer %>% filter(replicate == 11) 
Anc_scer12 <- Anc_scer %>% filter(replicate == 12) 
Anc_scer19 <- Anc_scer %>% filter(replicate == 19)


Con_scer<-fdata %>% filter(evolved=="Evolved_control")
Con_scer<-Con_scer %>% filter(specie=="1Scer")
Con_scer1 <- Con_scer %>% filter(replicate == 1) 
Con_scer3 <- Con_scer %>% filter(replicate == 3)
Con_scer5 <- Con_scer %>% filter(replicate == 5)
Con_scer11 <- Con_scer %>% filter(replicate == 11)
Con_scer18 <- Con_scer %>% filter(replicate == 18)


NQ_scer<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_scer<-NQ_scer %>% filter(specie=="1Scer")
NQ_scer3 <- NQ_scer %>% filter(replicate == 3) 
NQ_scer5 <- NQ_scer %>% filter(replicate == 5) 
NQ_scer11 <- NQ_scer %>% filter(replicate == 11)
NQ_scer13 <- NQ_scer %>% filter(replicate == 13) 
NQ_scer18 <- NQ_scer %>% filter(replicate == 18)
NQ_scer21 <- NQ_scer %>% filter(replicate == 21)

Anc_spar<-fdata %>% filter(evolved=="Ancestor")
Anc_spar<-Anc_spar %>% filter(specie=="2Spar")
Anc_spar1 <- Anc_spar %>% filter(replicate == 1) 
Anc_spar2 <- Anc_spar %>% filter(replicate == 2) 
Anc_spar4 <- Anc_spar %>% filter(replicate == 4)
Anc_spar23 <- Anc_spar %>% filter(replicate == 23)

Con_spar<-fdata %>% filter(evolved=="Evolved_control")
Con_spar<-Con_spar %>% filter(specie=="2Spar")
Con_spar1 <- Con_spar %>% filter(replicate == 1)
Con_spar2 <- Con_spar %>% filter(replicate == 2)
Con_spar4 <- Con_spar %>% filter(replicate == 4)
Con_spar23 <- Con_spar %>% filter(replicate == 23)
Con_spar26 <- Con_spar %>% filter(replicate == 26)

NQ_spar<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_spar<-NQ_spar %>% filter(specie=="2Spar")
NQ_spar2 <- NQ_spar %>% filter(replicate == 2)
NQ_spar5 <- NQ_spar %>% filter(replicate == 5)
NQ_spar16 <- NQ_spar %>% filter(replicate == 16)
NQ_spar23 <- NQ_spar %>% filter(replicate == 23)
NQ_spar24 <- NQ_spar %>% filter(replicate == 24)
NQ_spar28 <- NQ_spar %>% filter(replicate == 28)

Anc_hyb<-fdata %>% filter(evolved=="Ancestor")
Anc_hyb<-Anc_hyb %>% filter(specie=="3Hybrid")
Anc_hyb1 <- Anc_hyb %>% filter(replicate == 19)
Anc_hyb23 <- Anc_hyb %>% filter(replicate == 23) 
Anc_hyb24 <- Anc_hyb %>% filter(replicate == 24) 
Anc_hyb27 <- Anc_hyb %>% filter(replicate == 27)

Con_hyb<-fdata %>% filter(evolved=="Evolved_control")
Con_hyb<-Con_hyb %>% filter(specie=="3Hybrid")
Con_hyb20 <- Con_hyb %>% filter(replicate == 20)
Con_hyb23 <- Con_hyb %>% filter(replicate == 23) 
Con_hyb26 <- Con_hyb %>% filter(replicate == 26) 
Con_hyb27 <- Con_hyb %>% filter(replicate == 27)
Con_hyb28 <- Con_hyb %>% filter(replicate == 28)

NQ_hyb<-fdata %>% filter(evolved=="Evolved_NQO")
NQ_hyb<-NQ_hyb %>% filter(specie=="3Hybrid")
NQ_hyb1 <- NQ_hyb %>% filter(replicate == 1) 
NQ_hyb26 <- NQ_hyb %>% filter(replicate == 26) 
NQ_hyb27 <- NQ_hyb %>% filter(replicate == 27) 
NQ_hyb28 <- NQ_hyb %>% filter(replicate == 28)

tri1graph <-rbind(Anc_scer3,control3n,control2n)
tri2graph <-rbind(Anc_scer11,control3n,control2n)
tri3graph <-rbind(Anc_scer19,control3n,control2n)
triagraph <-rbind(Anc_scer5,control3n,control2n)
tribgraph <-rbind(Anc_scer12,control3n,control2n)

tri4graph <-rbind(Con_scer1,control3n,control2n)
tri5graph <-rbind(Con_scer11,control3n,control2n)
tricgraph <-rbind(Con_scer3,control3n,control2n)
tridgraph <-rbind(Con_scer5,control3n,control2n)
triegraph <-rbind(Con_scer18,control3n,control2n)

tri6graph <-rbind(NQ_scer5,control3n,control2n)
tri7graph <-rbind(NQ_scer11,control3n,control2n)
tri8graph <-rbind(NQ_scer13,control3n,control2n)
tri9graph <-rbind(NQ_scer21,control3n,control2n)
trifgraph <-rbind(NQ_scer3,control3n,control2n)
triggraph <-rbind(NQ_scer18,control3n,control2n)

tri10graph <-rbind(Anc_spar1,control3n,control2n)
tri11graph <-rbind(Anc_spar2,control3n,control2n)
tri12graph <-rbind(Anc_spar4,control3n,control2n)
trihgraph <-rbind(Anc_spar23,control3n,control2n)

triigraph <-rbind(Con_spar1,control3n,control2n)
trijgraph <-rbind(Con_spar2,control3n,control2n)
trikgraph <-rbind(Con_spar4,control3n,control2n)
trilgraph <-rbind(Con_spar23,control3n,control2n)
trimgraph <-rbind(Con_spar26,control3n,control2n)

tringraph <-rbind(NQ_spar23,control3n,control2n)
triograph <-rbind(NQ_spar24,control3n,control2n)
tripgraph <-rbind(NQ_spar2,control3n,control2n)
triqgraph <-rbind(NQ_spar5,control3n,control2n)
trisgraph <-rbind(NQ_spar16,control3n,control2n)
triugraph <-rbind(NQ_spar28,control3n,control2n)

tri15graph <-rbind(Anc_hyb19,control3n,control2n)
tri16graph <-rbind(Anc_hyb23,control3n,control2n)
tri17graph <-rbind(Anc_hyb24,control3n,control2n)
tri18graph <-rbind(Anc_hyb27,control3n,control2n)

tri19graph <-rbind(Con_hyb20,control3n,control2n)
tri20graph <-rbind(Con_hyb23,control3n,control2n)
tri21graph <-rbind(Con_hyb26,control3n,control2n)
tri22graph <-rbind(Con_hyb27,control3n,control2n)
tri23graph <-rbind(Con_hyb28,control3n,control2n)

tri24graph <-rbind(NQ_hyb10,control3n,control2n)
tri25graph <-rbind(NQ_hyb25,control3n,control2n)
tri26agraph <-rbind(NQ_hyb1,control3n,control2n)
tri26graph <-rbind(NQ_hyb26,control3n,control2n)
tri27graph <-rbind(NQ_hyb27,control3n,control2n)
tri28graph <-rbind(NQ_hyb28,control3n,control2n)

#Scer - Ancestor
tri1<- tri1graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+
  ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 3"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri2<- tri2graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+theme(strip.text.y = element_blank())

tri3<- tri3graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 19"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tria<- triagraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 5"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

trib<- tribgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 12"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


#Scer - Evolved in Control
tri4<- tri4graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "green4"),
                    labels = (c("2n", "Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri5<- tri5graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tric<- tricgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 3"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trid<- tridgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 5"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trie<- triegraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 18"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Scer - Evolved in UV mimetic
tri6<- tri6graph  %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 5"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri7<- tri7graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 11"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri8<- tri8graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 13"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri9<- tri9graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 21"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

trif<- trifgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 3"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

trig<- triggraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey","green4"),
                    labels = (c("2n", "Line 18"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Spar - Ancestor
tri10<- tri10graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n","Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri11<- tri11graph %>%filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 2"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri12<- tri12graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 4"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trih<- trihgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Spar - Evolved in control
trii<- triigraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trij<- trijgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 2"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trik<- trikgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 4"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tril<- trilgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trim<- trimgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 26"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Spar - Evolved in UV mimetic
trin<- tringraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

trio<- triograph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 24"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

trip<- tripgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 2"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

triq<- triqgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 5"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tris<- trisgraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 16"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

triu<- triugraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "dodgerblue1"),
                    labels = (c("2n", "Line 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Ancestor
tri16<- tri16graph  %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri17<- tri17graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 24"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri18<- tri18graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Evolved in Control
tri19<- tri19graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 20"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri20<- tri20graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 23"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri21<- tri21graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 26"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri22<- tri22graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tri23<- tri23graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#Hybrid - Evolved in UV mimetic
tri26a<- tri26agraph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 1"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


tri26<- tri26graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 26"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri27<- tri27graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 27"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

tri28<- tri28graph %>% filter(evolved!="3n")%>%
  ggplot(aes(log(FL1_A)), fill=evolved, group=interaction(same_Replicate,evolved)) +
  xlab("DNA content \n (Arb. units fluorescence)")+ylab("Cell count \n (density)")+
  scale_fill_manual(legend_title,values=c("grey", "#FF9999"),
                    labels = (c("2n", "Line 28"))) + xlim(12,16) +
  geom_density(alpha = 0.7, aes(group=interaction(same_Replicate,evolved),fill=evolved)) +
  theme_bw() +
  theme(panel.spacing = unit(0.72, "cm"))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=14, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#################
############Assemble and save Supplementary Figure 19############
tri123<- plot_grid(tri1,tria,tri2,trib,tri3, nrow=5)
ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces cerevisiae"), "Ancestor")), 
             size = 15)
scer_ancestor<-tri123
scer_ancestor_title <-plot_grid(title,scer_ancestor,nrow=2,rel_heights=c(0.2,2))

tri45<- plot_grid(tri4,tric,trid,tri5,trie, nrow=5)
title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces cerevisiae"), "Evolved in control")), 
             size = 15)
scer_con<-tri45
scer_con_title<-plot_grid(title,scer_con,nrow=2,rel_heights=c(0.2,2))

tri6789<- plot_grid(trif,tri6,tri7,tri8,trig,tri9, nrow=6)
title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces cerevisiae"), "Evolved in UV mimetic")), 
             size = 15)
scer_nq<-tri6789
scer_nq_title<-plot_grid(title,scer_nq,nrow=2,rel_heights=c(0.2, 2))

tri101112<- plot_grid(tri10,tri11,tri12,trih, nrow=4)
title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces paradoxus"), "Ancestor")), 
             size = 15)
spar_ancestor<-tri101112
spar_ancestor_title<-plot_grid(title,spar_ancestor,nrow=2,rel_heights=c(0.2, 2))

triijklm<- plot_grid(trii,trij,trik,tril,trim, nrow=5)
title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces paradoxus"), "Evolved in control")), 
             size = 15)
spar_control<-triijklm
spar_control_title<-plot_grid(title,spar_control,nrow=2,rel_heights=c(0.2, 2))

tri1314<- plot_grid(trip,triq,tris,trin,trio,triu, nrow=6)
title <- ggdraw() + 
  draw_label(expression(atop(bolditalic("Saccharomyces paradoxus"), "Evolved in UV mimetic")), 
             size = 15)
spar_nq<-tri1314
spar_nq_title<-plot_grid(title,spar_nq,nrow=2,rel_heights=c(0.2, 2))

tri161718<- plot_grid(tri16,tri17,tri18, nrow=3)
title <- ggdraw() + 
  draw_label(expression(atop(bold("Hybrid"), "Ancestor")), 
             fontface = 'bold', size = 15)
hyb_ancestor<-plot_grid(hybanc_read,tri161718,nrow=1)
hyb_ancestor_title<-plot_grid(title,hyb_ancestor,nrow=2,rel_heights=c(0.2, 2))

tri1920212223<- plot_grid(tri19,tri20,tri21,tri22,tri23, nrow=5)
title <- ggdraw() + 
  draw_label(expression(atop(bold("Hybrid"), "Evolved in control")), 
             fontface = 'bold', size = 15)
hyb_con<-plot_grid(hybcon_read,tri1920212223,nrow=1)
hyb_con_title<-plot_grid(title,hyb_con,nrow=2,rel_heights=c(0.08, 2))

tri262728<- plot_grid(tri26a,tri26,tri27,tri28, nrow=4)
title <- ggdraw() + 
  draw_label(expression(atop(bold("Hybrid"), "Evolved in UV mimetic")), 
             fontface = 'bold', size = 15)
hyb_nq<-plot_grid(hybnq_read,tri262728,nrow=1)
hyb_nq_title<-plot_grid(title,hyb_nq,nrow=2,rel_heights=c(0.2, 2))

A<-plot_grid(scer_ancestor_title,scer_con_title,scer_nq_title,ncol=3)
z<- plot_grid(spar_ancestor_title,spar_control_title,ncol=2)
B<-plot_grid(z,spar_nq_title,ncol=2,rel_widths = c(2,1))

AB<-plot_grid(A,B,nrow=1)

C<-plot_grid(hyb_ancestor_title,hyb_nq_title,ncol=1,rel_heights = c(1,1.5))
CD<-plot_grid(C,hyb_con_title,ncol=2)
BC <-plot_grid(AB,CD,ncol=1,rel_heights=c(1.1,1.4))

#make graph legend
set.seed(123)  
data <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  category = sample(c("Category 1", "Category 2", "Category 3"), 100, replace = TRUE)
)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

plot_legend <- ggplot(data, aes(x = x, y = y, fill= category)) +
  geom_tile(aes(fill=category))+
  theme_prism() + 
  theme_bw(base_size=24) +
  scale_fill_manual("",values=c("green4", "dodgerblue1", "#FF9999"), 
                    labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme_minimal() +
  theme(legend.position = "top") + 
  guides(color = guide_legend(nrow = 1)) +
  theme(
    legend.position = "top",       
    legend.text = element_text(size = 40),  
    legend.title = element_text(size = 40),  
    legend.key.size = unit(3, "cm") 
  )

leg <- get_legend(plot_legend)

BC1<-plot_grid(leg,BC,nrow=2,rel_heights = c(0.4,7))

ggsave (plot = BC1, filename = "Supplementary_Fig19_low_quality.jpg", units = "cm", device = "jpg",width =70, height =96, dpi = 300,bg = "white")
ggsave (plot = BC1, filename = "Supplementary_Fig19.png", units = "cm", device = "png",width =70, height =96, dpi = 800,bg = "white")
ggsave (plot = BC1, filename = "Supplementary_Fig19.jpg", units = "cm", device = "jpg",width =70, height =96, dpi = 800,bg = "white")
ggsave (plot = BC1, filename = "Supplementary_Fig19.svg", units = "cm", device = "svg",width =70, height =96, dpi = 500,bg = "white")
ggsave (plot = BC1, filename = "Supplementary_Fig19.pdf", units = "cm", device = "pdf",width =70, height =96, dpi = 500,bg = "white")
#################

############Figure Supplementary 20,21 and 22############
############Figure read coverage######## 
all_pools_samples <- Supp_raw_cov
all_pools_samples$log2 <- log2(all_pools_samples$mean_window)

all_pools_samples$gainloss <-ifelse(all_pools_samples$log2>7,"Gain",NA)
all_pools_samples$gainloss <-ifelse(all_pools_samples$log2<2,"Loss",all_pools_samples$gainloss)

all_pools_samples<-all_pools_samples %<>% mutate(chrom2 = ifelse(chrom=="utg351_pilon",17,
                                                                 ifelse(chrom=="utg1271_pilon",18,
                                                                        ifelse(chrom=="utg584_pilon",19,
                                                                               ifelse(chrom=="utg988_pilon",20,
                                                                                      ifelse(chrom=="utg69_pilon",21,
                                                                                             ifelse(chrom=="utg639_pilon",22,
                                                                                                    ifelse(chrom=="utg199_pilon",23,
                                                                                                           ifelse(chrom=="utg176_pilon",24,
                                                                                                                  ifelse(chrom=="utg1121_pilon",25,
                                                                                                                         ifelse(chrom=="utg245_pilon",26,
                                                                                                                                ifelse(chrom=="utg298_pilon",27,
                                                                                                                                       ifelse(chrom=="utg110_pilon",28,
                                                                                                                                              ifelse(chrom=="utg11_pilon",29,
                                                                                                                                                     ifelse(chrom=="utg210_pilon",30,
                                                                                                                                                            ifelse(chrom=="utg48_pilon",31,
                                                                                                                                                                   ifelse(chrom=="utg155_pilon",32,
                                                                                                                                                                          ifelse(chrom=="utg122_pilon",33,
                                                                                                                                                                                 ifelse(chrom=="utg675_pilon",34,
                                                                                                                                                                                        ifelse(chrom=="chrI",1,
                                                                                                                                                                                               ifelse(chrom=="chrII",2,
                                                                                                                                                                                                      ifelse(chrom=="chrIII",3,
                                                                                                                                                                                                             ifelse(chrom=="chrIV",4,
                                                                                                                                                                                                                    ifelse(chrom=="chrV",5,
                                                                                                                                                                                                                           ifelse(chrom=="chrVI",6,
                                                                                                                                                                                                                                  ifelse(chrom=="chrVII",7,
                                                                                                                                                                                                                                         ifelse(chrom=="chrVIII",8,
                                                                                                                                                                                                                                                ifelse(chrom=="chrIX",9,
                                                                                                                                                                                                                                                       ifelse(chrom=="chrX",10,
                                                                                                                                                                                                                                                              ifelse(chrom=="chrXI",11,
                                                                                                                                                                                                                                                                     ifelse(chrom=="chrXII",12,
                                                                                                                                                                                                                                                                            ifelse(chrom=="chrXIII",13,
                                                                                                                                                                                                                                                                                   ifelse(chrom=="chrXIV",14,
                                                                                                                                                                                                                                                                                          ifelse(chrom=="chrXV",15,
                                                                                                                                                                                                                                                                                                 ifelse(chrom=="chrXVI",16,NA)))))))))))))))))))))))))))))))))))










all_pools_samples1 <- filter(all_pools_samples, Unique_value=="Yes")
all_pools_samples1$Replicate <- as.numeric(all_pools_samples1$Replicate) 

#Scer - Ancestor
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

#Spar - Ancestor
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Ancestor")))))#+theme(strip.text.y = element_blank())

#Hybrid - Ancestor
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Ancestor") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Ancestor")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M1<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)

#Scer - Evolved in Control
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in control")))))#+theme(strip.text.y = element_blank())

#Spar - Evolved in Control
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Evolved in control")))))#+theme(strip.text.y = element_blank())

#Hybrid - Evolved in Control
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Evolved_control") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Evolved in control")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M2<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)

#Scer - Evolved in UV mimetic
exampleScer <- filter(all_pools_samples1, Strain=="1Scer")
ancestor <- exampleScer %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces cerevisiae"),bold(" - Evolved in UV mimetic")))))#+theme(strip.text.y = element_blank())

#Spar - Evolved in UV mimetic
exampleSpar <- filter(all_pools_samples1, Strain=="2Spar")
ancestor2 <- exampleSpar %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Saccharomyces paradoxus")+ theme(plot.title = element_text(face = "italic"))+
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold.italic"))+
  ggtitle(expression(atop(paste(bolditalic("Saccharomyces paradoxus"),bold(" - Evolved in UV mimetic")))))#+theme(strip.text.y = element_blank())

#Hybrid - Evolved in UV mimetic
exampleHybrid <- filter(all_pools_samples1, Strain=="3Hybrid")
ancestor3 <- exampleHybrid %>% filter(Type=="Evolved_NQO") %>% 
  ggplot(aes(x=interaction(start,as.numeric(chrom2)), y=log2, col = chrom)) +
  xlab("Chromosome") +
  geom_point(aes(col = chrom),size=0.00002) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey", "grey",
                               "grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey","grey","grey","grey",
                               "grey","grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Hybrid - Evolved in UV mimetic")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")+
  ylab("Mean read depth (log2)")+
  ylim(0,10)+
  scale_x_discrete(breaks=c("50000.1","50000.2","50000.3","50000.4","50000.5","50000.6","50000.7","50000.8",
                            "50000.9","50000.10","50000.11","50000.12","50000.13","50000.14","50000.15","50000.16",
                            "50000.17","50000.18","50000.19","50000.20",
                            "50000.21","50000.22","50000.23","50000.24",
                            "50000.25","50000.26","50000.27","50000.28",
                            "50000.29","50000.30","50000.31","50000.32",
                            "50000.33","50000.34"),
                   labels=c("Scer chrI","Scer chrII","Scer chrIII","Scer chrIV","Scer chrV",
                            "Scer chrVI","Scer chrVII","Scer chrVIII","Scer chrIX","Scer chrX",
                            "Scer chrXI","Scer chrXII","Scer chrXIII","Scer chrXIV","Scer chrXV","Scer chrXVI",
                            "Spar chrI","Spar chrII","Spar chrIII","Spar chrIV","Spar chrV",
                            "Spar chrVI","Spar chrVII","Spar chrVIII","Spar chrIX","Spar chrX",
                            "Spar chrXI","Spar chrXII 1/2","Spar chrXII 2/2","Spar chrXIII","Spar chrXIV","Spar chrXV 1/2","Spar chrXV 2/2","Spar chrXVI")) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 20),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15),
        axis.title=element_text(size=30,"bold"))+ facet_wrap(.~ Replicate, ncol=5)+
  theme(plot.title = element_text(hjust = 0.5, size=40, face = "bold"))

M3<- plot_grid(ancestor,ancestor2,ancestor3,nrow=3)
#################
############Assemble and save Supplementary Figure 20############
ggsave (plot = M1, filename = "Supplementary_Fig20_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig20.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig20.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig20.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M1, filename = "Supplementary_Fig20.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################
############Assemble and save Supplementary Figure 21############
ggsave (plot = M2, filename = "Supplementary_Fig21_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig21.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig21.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig21.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M2, filename = "Supplementary_Fig21.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################
############Assemble and save Supplementary Figure 22############
ggsave (plot = M3, filename = "Supplementary_Fig22_low_quality.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 300,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig22.png", units = "cm", device = "png",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig22.jpg", units = "cm", device = "jpg",width = 70, height =90, dpi = 800,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig22.svg", units = "cm", device = "svg",width = 70, height =90, dpi = 500,bg = "white")
ggsave (plot = M3, filename = "Supplementary_Fig22.pdf", units = "cm", device = "pdf",width = 70, height =90, dpi = 500,bg = "white")
#################
