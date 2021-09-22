######## Contribution of opsins and chromophores to cone pigment variation across populations of Lake Victoria cichlids 

# Clean Environment
rm(list = ls())

# load packages
library(tidyr)
library(tidyverse)
library(MASS)
library(lme4)
library(ggplot2)
library(stats)
library(effects)
library(splines)
library(xlsx)
library(car)
library(lattice)
library(grid)
library(gridExtra)
library(pbkrtest)
library(multcomp)

## load data
str(wild)
## define factors and numerical variables
wild <-  mutate(wild,
                island = factor(island),
                sample.id = factor(sample.id),
                phenotype=factor(phenotype))

## separate dataset for outlier checks 

#Makobe
wild.Mn <- wild[which(wild$island=="Makobe" & wild$phenotype=="N"),]
wild.Mp <- wild[which(wild$island=="Makobe" & wild$phenotype=="P"),]

#Anchor
wild.An <- wild[which(wild$island=="Anchor" & wild$phenotype=="N"),]
wild.Ap <- wild[which(wild$island=="Anchor" & wild$phenotype=="P"),]

#Kissenda
wild.Kn <- wild[which(wild$island=="Kissenda" & wild$phenotype=="N"),]
wild.Kp <- wild[which(wild$island=="Kissenda" & wild$phenotype=="P"),]

#Python
wild.Pn <- wild[which(wild$island=="Python" & wild$phenotype=="N"),]
wild.Pp <- wild[which(wild$island=="Python" & wild$phenotype=="P"),]

#Luanso
wild.Ln <- wild[which(wild$island=="Luanso" & wild$phenotype=="N"),]
wild.Lp <- wild[which(wild$island=="Luanso" & wild$phenotype=="P"),]

## outlier checks
source("http://goo.gl/UUyEzD")

## Makobe
outlierKD(wild.Mp,cyp) 
outlierKD(wild.Mp, sws2b)
outlierKD(wild.Mp, sws2a)
outlierKD(wild.Mp, rh2) 
outlierKD(wild.Mp, lws)

outlierKD(wild.Mn,cyp) 
outlierKD(wild.Mn, sws2b)
outlierKD(wild.Mn, sws2a)
outlierKD(wild.Mn, rh2)
outlierKD(wild.Mn, lws)

## Anchor
outlierKD(wild.Ap, cyp)
outlierKD(wild.Ap, sws2b)
outlierKD(wild.Ap, sws2a)
outlierKD(wild.Ap, rh2)
outlierKD(wild.Ap, lws)

outlierKD(wild.An, cyp)
outlierKD(wild.An, sws2b)
outlierKD(wild.An, sws2a)
outlierKD(wild.An, rh2)
outlierKD(wild.An, lws)

## Kissenda
outlierKD(wild.Kp, cyp) 
outlierKD(wild.Kp, sws2b)
outlierKD(wild.Kp, sws2a)
outlierKD(wild.Kp, rh2)
outlierKD(wild.Kp, lws)

outlierKD(wild.Kn, cyp) 
outlierKD(wild.Kn, sws2b)
outlierKD(wild.Kn, sws2a)
outlierKD(wild.Kn, rh2)
outlierKD(wild.Kn, lws)

## Python
outlierKD(wild.Pp, cyp) 
outlierKD(wild.Pp, sws2b) 
outlierKD(wild.Pp, sws2a)
outlierKD(wild.Pp, rh2)
outlierKD(wild.Pp, lws)

outlierKD(wild.Pn, cyp) 
outlierKD(wild.Pn, sws2b)
outlierKD(wild.Pn, sws2a)
outlierKD(wild.Pn, rh2)
outlierKD(wild.Pn, lws)

## Luanso
outlierKD(wild.Lp, cyp) 
outlierKD(wild.Lp, sws2b) 
outlierKD(wild.Lp, sws2a)
outlierKD(wild.Lp, rh2)
outlierKD(wild.Lp, lws)

outlierKD(wild.Ln, cyp) 
outlierKD(wild.Ln, sws2b)
outlierKD(wild.Ln, sws2a)
outlierKD(wild.Ln, rh2)
outlierKD(wild.Ln, lws)

## recombine all data into dataset excluding outliers (1,5 *IQR)
d <- rbind(wild.Mn, wild.Mp,
           wild.An, wild.Ap,
           wild.Kn, wild.Kp,
           wild.Pn, wild.Pp,
           wild.Ln,wild.Lp)

## recombine data excluding Lusanso island & excluding outliers (1.5 * IQR) 
d_n <- rbind(wild.Mn, wild.Mp,
             wild.An, wild.Ap,
             wild.Kn, wild.Kp,
             wild.Pn, wild.Pp)

### calculate mean & sd
aggregate(d[, 16], list(d$phenotype:d$island), mean,na.rm = TRUE)
aggregate(d[, 16], list(d$phenotype:d$island), sd,na.rm = TRUE)

###### new working dataset ###### 
######################### Figures #########################

##  Figure 2:  Cyp27c1 expression in sympatric Pundamilia phenotypes at 5 locations

## separate datasets by phenotype
P <- d[which(d$phenotype=="P"),]
N <- d[which(d$phenotype=="N"),]

## Blue phenotype 
ggplot(P, aes(spectra, cyp, fill=island)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_fill_manual(values=c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue","cornflowerblue"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.025))+
  theme(legend.position="none",legend.background = element_rect(colour = 'white'),
        legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"),legend.text=element_text(size=8),legend.key.size = unit(0.5, "cm"))+
  theme(axis.text.x = element_text(size=15))+ 
  theme(axis.text.y = element_text(size=15))+
  labs(x="", y=" ") +
  scale_x_reverse()

## Red phenotype 
ggplot(N, aes(spectra, cyp, fill=island,alpha=factor(island))) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_alpha_manual(values = c(1, 1, 0.2,1,1))+
  scale_fill_manual(values=c("lightsalmon2", "lightsalmon2", "lightsalmon2", "lightsalmon2","lightsalmon2"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.025))+
  theme(legend.position="none",legend.background = element_rect(colour = 'white'),
        legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"),legend.text=element_text(size=8),legend.key.size = unit(0.5, "cm"))+
  theme(axis.text.x = element_text(size=15))+ 
  theme(axis.text.y = element_text(size=15))+
  labs(x="", y=" ") +
  scale_x_reverse()

## Figure 3: No covariation between cyp27c1 expression and population-level orange ratio (from Seehausen et al. 2008)
## load data with mean & standard error (previously calculated; combined in excel as new dataset)
fig_1<- read_excel("wild_data_final.xlsx",2)
str(fig_1)

## exclude Luanso island
fig_1.1 <- fig_1[-which(fig_1$island=="Luanso"),] 
fig_1.1 <-  mutate(fig_1.1,
                 island = factor(island),
                 phenotype=factor(phenotype))

ggplot(fig_1.1, aes(OR_pop, cyp_mean, colour=phenotype))+
  geom_point(aes(shape=island), size=6)+
  geom_pointrange(aes(ymin=cyp_mean-cyp_se, ymax=cyp_mean+cyp_se))+
  geom_pointrange(aes(xmin=OR_pop-OR_pop_se, xmax=OR_pop+OR_pop_se))+
  scale_colour_manual(values=c("blue3","red3"))+
  scale_shape_manual(values=c(8, 15, 17, 19))+
  ylim(0,0.015)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  labs(x=" ", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position = "none")

## Figure S3: Relationship between cyp27c1 and orange ratio based on the capture depth of the individuals used in this study
ggplot(fig_1.1, aes(OR_ind, cyp_mean, colour=phenotype))+
  geom_point(aes(shape=island), size=6)+
  geom_pointrange(aes(ymin=cyp_mean-cyp_se, ymax=cyp_mean+cyp_se))+
  geom_pointrange(aes(xmin=OR_ind-OR_ind_se, xmax=OR_ind+OR_ind_se))+
  scale_colour_manual(values=c("blue3","red3"))+
  scale_shape_manual(values=c(8, 15, 17, 19))+
  ylim(0,0.015)+
  theme_linedraw()+ ### creates boxes
  theme(panel.grid = element_blank())+
  labs(x=" ", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position = "none")

## Figure 4a: Cyp27c1 expression in sympatric blue and red Pundamilia phenotypes (excluding Luanso island)

ggplot(d_n, aes(island, cyp, fill=phenotype)) +
  facet_grid(~island,scales = "free")+
  theme(strip.text.x = element_text(size=10,face="bold"),
        strip.background = element_rect(colour="white",fill="white"))+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_fill_manual(values=c("cornflowerblue", "lightsalmon2"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.025))+
  theme(legend.position="none",legend.background = element_rect(colour = 'white'),
        legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"),legend.text=element_text(size=8),legend.key.size = unit(0.5, "cm"))+
  labs(y="")+
  labs(x=" ")+
  theme(axis.text.y = element_text(size=15))+
  theme(axis.text.x = element_text(size=15))

## Figure 4b  Difference in opsin gene expression patterns between phenotypes (excluding Luanso island) 
ggplot(d_n, aes(x = opsin, y = opsin_expression, fill=phenotype)) +
  facet_grid(~island) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_fill_manual(values=c("cornflowerblue", "lightsalmon2"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1.1))+
  theme( legend.position="none",strip.background = element_blank(),
         strip.text.x = element_blank())+
  theme(axis.text.x = element_text(size=15))+ 
  theme(axis.text.y = element_text(size=15))+
  labs(x="", y=" ") 


## Figure 5:  Relationship between cyp27c1 and opsin gene expression (including population means & individual data points)

## load data 
individual <-  mutate(individual,
                  island = factor(island),
                  phenotype = factor(phenotype))

population <-  mutate(population,
               island = factor(island),
               phenotype = factor(phenotype))

## define factors and numerical variables

## LWS
ggplot(population, aes(lws , cyp, colour=phenotype, shape=island, group=island)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_pointrange(aes(ymin=cyp-cyp_se, ymax=cyp+cyp_se), width=.1) +
  geom_pointrange(aes(xmin=lws -lws_se, xmax=lws +lws_se)) +
  geom_line(aes(linetype=island, colour=island),size=0.7)+
  geom_point(aes(shape=island), size=4)+
  geom_point(data = individual, alpha=0.4, size=3)+
  scale_shape_manual(values=c(17, 8, 16, 15))+
  scale_colour_manual(values=c("black" ,"black","black" ,"red3","blue3" ,"black"))+
  scale_linetype_manual(values=c("dotted", "solid","dotted", "solid"))+
  labs(x="", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  ylim(0,0.015)+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text( size=10, face="bold"),
        axis.title.y = element_text( size=10, face="bold"))

## Rh2
ggplot(population, aes(rh2 , cyp, colour=phenotype, shape=island, group=island)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_pointrange(aes(ymin=cyp-cyp_se, ymax=cyp+cyp_se), width=.1) +
  geom_pointrange(aes(xmin=rh2 -rh2.se, xmax=rh2 +rh2.se)) +
  geom_line(aes(linetype=island, colour=island),size=0.7)+
  geom_point(aes(shape=island), size=4)+
  geom_point(data = individual, alpha=0.4, size=3)+
  scale_shape_manual(values=c(17, 8, 16, 15))+
  scale_colour_manual(values=c("black" ,"black","black" ,"red3","blue3" ,"black"))+
  scale_linetype_manual(values=c("dotted", "solid","dotted", "solid"))+
  labs(x="", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  ylim(0,0.015)+
  theme(legend.position = "none")

## SWS2a
ggplot(population, aes(sws2a , cyp, colour=phenotype, shape=island, group=island)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_pointrange(aes(ymin=cyp-cyp_se, ymax=cyp+cyp_se), width=.1) +
  geom_pointrange(aes(xmin=sws2a -sws2a.se, xmax=sws2a +sws2a.se)) +
  geom_line(aes(linetype=island, colour=island),size=0.7)+
  geom_point(aes(shape=island), size=4)+
  geom_point(data = individual, alpha=0.4, size=3)+
  scale_shape_manual(values=c(17, 8, 16, 15))+
  scale_colour_manual(values=c("black" ,"black","black" ,"red3","blue3" ,"black"))+
  scale_linetype_manual(values=c("dotted", "solid","dotted", "solid"))+
  labs(x="", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  ylim(0,0.015)+
  theme(legend.position = "none")

## SWS2b
ggplot(population, aes(sws2b , cyp, colour=phenotype, shape=island, group=island)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_pointrange(aes(ymin=cyp-cyp_se, ymax=cyp+cyp_se), width=.1) +
  geom_pointrange(aes(xmin=sws2b -sws2b.se, xmax=sws2b +sws2b.se)) +
  geom_line(aes(linetype=island, colour=island),size=0.7)+
  geom_point(aes(shape=island), size=4)+
  geom_point(data = individual, alpha=0.4, size=3)+
  scale_shape_manual(values=c(17, 8, 16, 15))+
  scale_colour_manual(values=c("black" ,"black","black" ,"red3","blue3" ,"black"))+
  scale_linetype_manual(values=c("dotted", "solid","dotted", "solid"))+
  labs(x="", y="")+
  theme(axis.text.x = element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+
  ylim(0,0.015)+
  theme(legend.position = "none")

## Figure S4: Relationship between cyp27c1 expression and opsin gene expression at individual level

##LWS
ggplot(d_n, aes(lws , cyp, group=species)) +
  facet_grid(~island)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_smooth(method=lm, se=T, color="grey56", size=0.5)+
  theme(panel.grid = element_blank())+ 
  geom_point(aes(colour=species), size=1.5)+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("red3","blue3"))+
  theme(axis.text.x = element_text(size=11))+
  theme(axis.text.y = element_text(size=11))+
  labs(x="", y="")+
  xlim(0.4,0.9)

##Rh2
ggplot(d_n, aes(rh2 , cyp, group=species)) +
  facet_grid(~island)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_smooth(method=lm, se=T, color="grey56", size=0.5)+
  theme(panel.grid = element_blank())+ 
  geom_point(aes(colour=species), size=1.5)+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("red3","blue3"))+
  theme(axis.text.x = element_text(size=11))+
  theme(axis.text.y = element_text(size=11))+
  labs(x="", y="")+
  xlim(0,0.4)

##SWS2b
ggplot(d_n, aes(sws2b , cyp, group=species)) +
  facet_grid(~island)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_smooth(method=lm, se=T, color="grey56", size=0.5)+
  theme(panel.grid = element_blank())+ 
  geom_point(aes(colour=species), size=1.5)+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("red3","blue3"))+
  theme(axis.text.x = element_text(size=11))+
  theme(axis.text.y = element_text(size=11))+
  labs(x="", y="")+
  xlim(0, 0.015)

##SWS2a
ggplot(d_n, aes(sws2a , cyp, group=species)) +
  facet_grid(~island)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_smooth(method=lm, se=T, color="grey56", size=0.5)+
  theme(panel.grid = element_blank())+ 
  geom_point(aes(colour=species), size=1.5)+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("red3","blue3"))+
  theme(axis.text.x = element_text(size=11))+
  theme(axis.text.y = element_text(size=11))+
  labs(x="", y="")+
  xlim(0,0.02)


ggplot(d_n, aes(lws , rh2, group=species)) +
  facet_grid(~island)+
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+ 
  geom_smooth(method=lm, se=T, color="grey56", size=0.5)+
  theme(panel.grid = element_blank())+ 
  geom_point(aes(colour=species), size=1.5)+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("red3","blue3"))+
  theme(axis.text.x = element_text(size=11))+
  theme(axis.text.y = element_text(size=11))+
  labs(x="", y="")+
  xlim(0.4, 0.9)

## Figure 6: Variation in A2 proportion does not maximize visual performance for Pundamilia phenotypes at three locations (based on population capture depth; from Seehausen et al. 2008))

## load data on quantum catch estimates (visual modeling work done in excel)
## define factors and numerical variables
## define factors and numerical variables
qc <-  mutate(qc,
              island = factor(island),
              phenotype=factor(phenotype),
              qc=factor(qc),
              genotype=factor(genotype))

limits <- aes(ymax = mean_pop + se_pop, ymin=mean_pop - se_pop)
ggplot(qc, aes(x = phenotype, y = mean_pop, fill=phenotype:qc)) +
  facet_grid(~island)+
  theme_light()+
  theme(panel.grid = element_blank())+
  theme(strip.text.x = element_text(size=10,face="bold"),strip.background = element_rect(colour="white",fill="white"))+
  theme(axis.text.x = element_text(size=9))+
  theme(axis.text.y = element_text(size=9))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 100),breaks=c(0,20,40,60,80,100))+ 
  scale_x_discrete(expand = c(0.15,0.15))+
  geom_bar(stat = "identity",position=position_dodge(width = 0.64),width=0.59)+
  geom_errorbar(limits2, position=position_dodge(width=0.65), width=0.15)+
  scale_fill_manual(values=c("lightskyblue3", "dodgerblue", "blue","red3","indianred2","indianred"))+
  theme(legend.background = element_rect(colour = 'white'),legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"))+
  labs(y="")+
  theme(plot.title = element_text(size=10, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5,vjust=-0.1))+theme(axis.text.x = element_text(size=13,face="bold"))+ 
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position='',legend.background = element_rect(colour = 'white'),
        legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"))

### Figure S5:Quantum catch estimates based on individual capture depth (individuals captured in this study)
limits <- aes(ymax = mean + se, ymin=mean - se)
ggplot(qc, aes(x = phenotype, y = mean, fill=phenotype:qc)) +
  facet_grid(~island)+
  theme_light()+
  theme(panel.grid = element_blank())+
  theme(strip.text.x = element_text(size=10,face="bold"),strip.background = element_rect(colour="white",fill="white"))+
  theme(axis.text.x = element_text(size=9))+
  theme(axis.text.y = element_text(size=9))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 100),breaks=c(0,20,40,60,80,100))+ 
  scale_x_discrete(expand = c(0.15,0.15))+
  geom_bar(stat = "identity",position=position_dodge(width = 0.64),width=0.59)+
  geom_errorbar(limits, position=position_dodge(width=0.65), width=0.15)+
  scale_fill_manual(values=c("lightskyblue3", "dodgerblue", "blue","red3","indianred2","indianred"))+
  theme(legend.background = element_rect(colour = 'white'),legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"))+
  labs(y="")+
  theme(plot.title = element_text(size=10, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5,vjust=-0.1))+theme(axis.text.x = element_text(size=13,face="bold"))+ 
  theme(axis.text.y = element_text(size=15))+
  theme(legend.position='',legend.background = element_rect(colour = 'white'),
        legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm"),
        panel.spacing = unit(0.5,"lines"))

### Figure 7: Cyp27c1 expression in wild-caught and laboratory-reared fish - from Python island

## load data
labwild <- mutate(labwild,
                  treatment=factor(treatment),
                  origin=factor(origin),
                  family=factor(family),
                  phenotype=factor(phenotype))
#### separate data
d.wildN <- labwild[which(labwild$origin=="wild"& labwild$phenotype=="N"),]
d.wildP <- labwild[which(labwild$origin=="wild"& labwild$phenotype=="P"),]
d.labP.s <- labwild[which(labwild$origin=="lab"& labwild$phenotype=="P" & labwild$treatment=="shallow"),]
d.labP.d <- labwild[which(labwild$origin=="lab"& labwild$phenotype=="P"& labwild$treatment=="deep"),]
d.labN.s <- labwild[which(labwild$origin=="lab"& labwild$phenotype=="N"& labwild$treatment=="shallow"),]
d.labN.d <- labwild[which(labwild$origin=="lab"& labwild$phenotype=="N"& labwild$treatment=="deep"),]

## outlier checks
source("http://goo.gl/UUyEzD")

outlierKD(d.wildN,cyp) 
outlierKD(d.wildP,cyp)

outlierKD(d.labP.s,cyp) 
outlierKD(d.labP.d,cyp) 
outlierKD(d.labN.s,cyp) 
outlierKD(d.labN.d,cyp)


## Figure 7a
## recombine data 
labvswild <- rbind(d.wildN, d.wildP,d.labP.d, d.labP.s, d.labN.s,d.labN.d)
lab <- rbind(d.labP.d, d.labP.s, d.labN.s,d.labN.d)

### compare cyp27c1 expression from wild-caught fish to lab-reared fish
ggplot(labvswild, aes(x = origin, y = cyp, fill=phenotype)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_fill_manual(values=c("cornflowerblue", "lightsalmon2"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.025))+
  theme( legend.position="none",strip.background = element_blank(),
         strip.text.x = element_blank())+
  theme(axis.text.y = element_text(size=15))+ 
  theme(axis.text.x= element_blank())+
  labs(x="", y=" ") 

### compare lab data (broad vs. red-shifted light)
ggplot(lab, aes(x = treatment, y = cyp, fill=phenotype)) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank())+
  geom_boxplot ( )+
  scale_fill_manual(values=c("lightskyblue", "gold3"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.025))+
  theme( legend.position="none",strip.background = element_blank(),
         strip.text.x = element_blank())+
  theme(axis.text.y = element_text(size=15))+ 
  theme(axis.text.x= element_blank())+
  labs(x="", y=" ") 

######################### Statistical analysis #########################

## 1) Does cyp27c1 expression decrease with water transparency?
### Blue phenotype 
blue_water <- glm(cyp~spectra,data=P)
plot(blue_water,which=1)
plot(allEffects(blue_water))
Anova(blue_water,type=3,singular.ok = TRUE) ## not significant

### Red phenotype 
red_water <- glm(cyp~spectra,data=N)
plot(red_water,which=1)
plot(allEffects(red_water))
Anova(red_water,type=3,singular.ok = TRUE) ##  significant

## 2) Difference between islands and/or phenoytpes? & 2) Can cyp27c1 expression be predicted by the local light environment (i.e. orange ratio)?
##  data distribution
hist(log(d$cyp))# looks good if log transformed
qqnorm(log(d$cyp)) 
qqline(log(d$cyp))

# log transform data
lcyp <- log(d$cyp)

## individual level orange ratio
cyp1.1 <- lm(lcyp~species*island+OR,data=d) 
drop1(cyp1.1,test="Chisq") 

cyp2.1 <- update(cyp1.1,.~.-OR)
drop1(cyp2.1,test="Chisq") 
Anova(cyp2.1,type=3,singular.ok = TRUE) 

cyp2.2 <- update(cyp2.1,.~.-species:island)
drop1(cyp2.2,test="Chisq") 

cyp2.3 <- update(cyp2.2,.~.-species)
drop1(cyp2.3,test="Chisq")

cyp2.4 <- update(cyp2.3,.~.-island)
drop1(cyp2.4,test="Chisq") 

cyp2.5 <- update(cyp2.4,.~.+OR)
drop1(cyp2.5,test="Chisq") 
Anova(cyp2.5) 

cyp2.6 <- update(cyp2.4,.~.+species)
drop1(cyp2.6,test="Chisq") 
Anova(cyp2.6) 

cyp2.7 <- update(cyp2.4,.~.+island)
drop1(cyp2.7,test="Chisq") 
Anova(cyp2.7) 

#### population level orange ratio 
cyp3 <- lm(lcypl~species*island+OR_pop,data=d)
plot(cyp3,which=1)
hist(resid(cyp3))

drop1(cyp3,test="Chisq") 
cyp3.1 <- update(cyp3,.~.-OR_pop)
drop1(cyp3.1,test="Chisq") 
cyp3.2 <- update(cyp3.1,.~.-species:island)
drop1(cyp3.2,test="Chisq")

cyp3.3 <- update(cyp3.2,.~.-species)
drop1(cyp3.3,test="Chisq")

cyp3.4 <- update(cyp3.3,.~.-island)
drop1(cyp3.4,test="Chisq")

cyp3.5 <- update(cyp3.4,.~.+OR_pop)
drop1(cyp3.5,test="Chisq")
Anova(cyp3.5)

## posthoc test using multcomp - package - Within island between phenotype difference in cyp27c1 expression? 
library(multcomp)
d$inter <- interaction (d$island,d$phenotype)
m1 <- lm(log(cyp)~inter,data=d)
summary(m1)
posthoc1<-glht(m1, linfct=mcp(inter="Tukey"))
summary(posthoc1) 

## set contrasts
contrMat(table(d$inter))
contrast.matrix <- rbind(  "Anchor:N vs. Anchor:P" =       c(0, 0, 0, 0, 0, 1, 0, 0,0,0),
                           "Kissenda:N vs. Kissenda:P" =   c(0, -1, 0, 0, 0, 0, 1, 0,0,0),
                           "Makobe:N vs. Makobe:P" =       c(0, 0, 0, -1, 0, 0, 0, 0,1,0),
                           "Python:N vs. Python:P" =       c(0, 0, 0, 0, -1, 0, 0, 0,0,1))
summary(glht(m1, contrast.matrix, type="Tukey",test = adjusted("bonferroni"))) 

##### 4) Relationship between cyp27c1 and opsin gene expression at individual level (excluding Luanso island)
#seperate dataset
#Makobe
d.MN <- d[which(d$island=="Makobe" & d$phenotype=="N"),]
d.MP <- d[which(d$island=="Makobe" & d$phenotype=="P"),]

#Anchor
d.An <- d[which(d$island=="Anchor" & d$phenotype=="N"),]
d.Ap <-d[which(d$island=="Anchor" &d$phenotype=="P"),]

#Kissenda
d.Kn <- d[which(d$island=="Kissenda" & d$phenotype=="N"),]
d.Kp <- d[which(d$island=="Kissenda" & d$phenotype=="P"),]

#Python
d.Pn <-d[which(d$island=="Python" & d$phenotype=="N"),]
d.Pp <- d[which(d$island=="Python" &d$phenotype=="P"),]

##LWS
summary(lm(cyp~lws,data=d.MN))
summary(lm(cyp~lws,data=d.An)) 
summary(lm(cyp~lws,data=d.Kn))
summary(lm(cyp~lws,data=d.Pn)) 

summary(lm(cyp~lws,data=d.MP)) 
summary(lm(cyp~lws,data=d.Ap)) 
summary(lm(cyp~lws,data=d.Kp)) 
summary(lm(cyp~lws,data=d.Pp)) 

##RH2
summary(lm(cyp~rh2,data=d.MN)) 
summary(lm(cyp~rh2,data=d.An)) 
summary(lm(cyp~rh2,data=d.Kn)) 
summary(lm(cyp~rh2,data=d.Pn)) 

summary(lm(cyp~rh2,data=d.MP)) 
summary(lm(cyp~rh2,data=d.Ap)) 
summary(lm(cyp~rh2,data=d.Kp)) 
summary(lm(cyp~rh2,data=d.Pp)) 

##SWS2a
summary(lm(cyp~sws2a,data=d.MN)) 
summary(lm(cyp~sws2a,data=d.An)) 
summary(lm(cyp~sws2a,data=d.Kn)) 
summary(lm(cyp~sws2a,data=d.Pn)) 

summary(lm(cyp~sws2a,data=d.MP)) 
summary(lm(cyp~sws2a,data=d.Ap)) 
summary(lm(cyp~sws2a,data=d.Kp)) 
summary(lm(cyp~sws2a,data=d.Pp)) 

##SWS2b
summary(lm(cyp~sws2b,data=d.MN)) 
summary(lm(cyp~sws2b,data=d.An)) 
summary(lm(cyp~sws2b,data=d.Kn)) 
summary(lm(cyp~sws2b,data=d.Pn)) 
cor.test(d.Pn$sws2b,d.Pn$cyp)

summary(lm(cyp~sws2b,data=d.MP)) 
summary(lm(cyp~sws2b,data=d.Ap)) 
summary(lm(cyp~sws2b,data=d.Kp)) 
summary(lm(cyp~sws2b,data=d.Pp)) 
cor.test(d.MP$sws2a,d.MP$cyp)

## 5) Do the observed patterns in cyp27c1 and opsin expression optimize visual performance (i.e. quantum catch estimates)? 
##### Quantum catch -  statistical analysis

#seperate data

qc_statMn <- qc_stat[which(qc_stat$island=="Makobe"& qc_stat$species=="NN"),]
qc_statMp <- qc_stat[which(qc_stat$island=="Makobe"& qc_stat$species=="PP"),]

qc_statPn <- qc_stat[which(qc_stat$island=="Python"& qc_stat$species=="NN"),]
qc_statPp <- qc_stat[which(qc_stat$island=="Python"& qc_stat$species=="PP"),]

qc_statKn <- qc_stat[which(qc_stat$island=="Kissenda"& qc_stat$species=="NN"),]
qc_statKp <- qc_stat[which(qc_stat$island=="Kissenda"& qc_stat$species=="PP"),]

#### Quantum catch estimates based on population depth
r<- lm(output_popdepth~qc*phenotype*island,data=d)
Anova(r)
drop1(r,test="Chisq") 

r.1 <- update(r,.~.-qc:phenotype:island)
drop1(r.1,test="Chisq") 

r.2 <- update(r.1,.~.-qc:island)
drop1(r.2,test="Chisq") 

r.3 <- update(r.2,.~.-qc:phenotype)
drop1(r.3,test="Chisq") 

r.4 <- update(r.3,.~.-phenotype:island)
drop1(r.4,test="Chisq") 

r.5 <- update(r.4,.~.-phenotype)
drop1(r.5,test="Chisq") 

r.6 <- update(r.5,.~.-island)
drop1(r.6,test="Chisq") 

r.7 <- update(r.6,.~.-qc)
drop1(r.7,test="Chisq") 

r.8 <- update(r.7,.~.+qc)
drop1(r.8,test="Chisq")
Anova(r.8)

r.9 <- update(r.8,.~.+phenotype)
drop1(r.9,test="Chisq")

## Quantum catch estimates based based on individual depth 
i<- lm(logit(output_inddepth)~qc*phenotype*island,data=d)
Anova(i)
drop1(r,test="Chisq") 

i.1 <- update(r,.~.-qc:phenotype:island)
drop1(i.1,test="Chisq") 

i.2 <- update(i.1,.~.-qc:island)
drop1(i.2,test="Chisq") 

i.3 <- update(i.2,.~.-qc:phenotype)
drop1(i.3,test="Chisq") 

i.4 <- update(i.3,.~.-phenotype:island)
drop1(i.4,test="Chisq") 

i.5 <- update(i.4,.~.-phenotype)
drop1(i.5,test="Chisq") 

i.6 <- update(i.5,.~.-island)
drop1(i.6,test="Chisq") 

i.7 <- update(i.6,.~.-qc)
drop1(i.7,test="Chisq") 

i.8 <- update(i.7,.~.+qc)
drop1(i.8,test="Chisq")
Anova(i.8)

i.9 <- update(i.8,.~.+phenotype)
drop1(i.9,test="Chisq")

## 6) Difference between wild-caught and lab-reared fish (Python island)
## only lab data from shallow P and deep N
### check data distribution 
hist(log(labvswild$cyp))
qqnorm(log(labvswild$cyp)) 
qqline(log(labvswild$cyp))

labvswild_1 <- glm(cyp~phenotype*treatment+origin*phenotype,data=labvswild)
Anova(test1)

labvswild_1.1 <- update(labvswild_1,.~.-phenotype:treatment)
drop1(labvswild_1.1,test="Chisq")

labvswild_1.2 <- update(labvswild_1.1,.~.-treatment)
drop1(labvswild_1.2,test="Chisq")
Anova(labvswild_1.2)

labvswild_1.3 <- update(labvswild_1.2,.~.-phenotype:origin)
drop1(labvswild_1.3,test="Chisq")

labvswild_1.4 <- update(labvswild_1.3,.~.-origin)
Anova(labvswild_1.4)

labvswild_1.5 <- update(labvswild_1.4,.~.-phenotype)
Anova(labvswild_1.5)

#### posthoc test
labvswild$inter_lab <- interaction (labvswild$origin,labvswild$phenotype)
labwild.6 <- lm(log(cyp)~inter_lab,data=labvswild)
summary(labwild.6)
library(multcomp)

## set contrasts
contrast.matrix <- rbind(  "wild:NN vs. lab:NN" =       c(0, 1, 0, 0),
                           "wild:PP vs. lab:PP" =       c(0, 0, -1, 1))
summary(glht(labwild.6, contrast.matrix, type="Tukey",test = adjusted("bonferroni"))) 

####  Does rearing light influence cyp27c1 expression in laboratory-reared fish (analysis only including laboratory-reared fish)?
lab.1 <- lmer(log(cyp)~phenotype*treatment+(1|mother_ID)+(1|father_ID),data=lab)
Anova(lab.1)

lab.2 <- update(lab.1,.~.-phenotype:treatment)
drop1(lab.2,test="Chisq")

lab.3 <- update(lab.2,.~.-phenotype)
drop1(lab.3,test="Chisq")

lab.4 <- update(lab.3,.~.-treatment)
drop1(lab.4,test="Chisq")

lab.5 <- update(lab.4,.~.+treatment)
Anova(lab.5)

lab.6 <- update(lab.5,.~.+phenotype)
Anova(lab.6)



