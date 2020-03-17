library(searchable)
library(reshape)
library(rcarbon)
library(ggpubr)
library(zoo)
library(outliers)
library(broom)
library(robustbase)
library(moments)
library(matrixStats)
library(readr)
library(psych)
library(gmodels)
library(plyr)
library(tidyverse)

####Note: C:/Users/Darcy/Dropbox/R/ to ~/Dropbox/R/

setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2")
RawData <- read_csv("Bird2020_Holocene.csv")
AgID <- read_csv("AgID.csv")
Directory <- read_csv("Directory_Sbox.csv")
Summary <- read_csv("Summary.csv")

Directory$ag <- as.factor(Directory$ag)

cbPalette <- c("#D55E00" , "#009E73", "#56B4E9")
cbbPalette <- c("#56B4E9", "#D55E00")
agnames <- list(
  "0" = "Hunter-Gatherers",
  "1" = "Agriculturalists"
)
variable_labeller <- function(variable,value){
  return(agnames[value])
}

#### Just Booms/Busts ----
library(plyr) ##This is necessary to make the code, as it currently is, work.
setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/3_FirstDiff/")

Dif50long <- read_csv ("Dif50long.csv")
Dif100long <- read_csv ("Dif100long.csv")
Dif200long <- read_csv ("Dif200long.csv")


#### Boom/Bust Density plots ====
Dif50long_pos <- subset(Dif50long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif50long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif50long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif50long_pos, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif50long_pos, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette)+
  labs(x="Positive 50 Year First Difference Values (Booms)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox50_posamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
Dif100long_pos <- subset(Dif100long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif100long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif100long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif100long_pos, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif100long_pos, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Positive 100 Year First Difference Values (Booms)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox100_posamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last
Dif200long_pos <- subset(Dif200long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif200long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif200long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif200long_pos, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif200long_pos, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Positive 200 Year First Difference Values (Booms)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox200_posamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd


##### 7.22 Negative now.
Dif50long_neg <- subset(Dif50long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif50long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif50long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif50long_neg, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif50long_neg, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette)+
  labs(x="Negative 50 Year First Difference Values (Busts)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox50_negamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
Dif100long_neg <- subset(Dif100long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif100long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif100long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif100long_neg, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif100long_neg, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Negative 100 Year First Difference Values (Busts)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox100_negamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last
Dif200long_neg <- subset(Dif200long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-plyr::ddply(Dif200long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif200long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif200long_neg, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif200long_neg, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Negative 200 Year First Difference Values (Busts)", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox200_negamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd


## Boom/Bust Skewness ====

library(tidyverse)
##Pos/Neg 50
ampHGpos <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value > 0)

ampHGneg <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value < 0)

ampAGpos <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value > 0)

ampAGneg <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value < 0)

skewness(ampHGpos$value)
skewness(ampAGpos$value)

skewness(abs(ampHGneg$value))
skewness(abs(ampAGneg$value))


##Pos/Neg 100
ampHGpos <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value > 0)

ampHGneg <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value < 0)

ampAGpos <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value > 0)

ampAGneg <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value < 0)

skewness(ampHGpos$value)
skewness(ampAGpos$value)

skewness(abs(ampHGneg$value))
skewness(abs(ampAGneg$value))


##Pos/Neg 200
ampHGpos <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value > 0)

ampHGneg <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "0") %>%
  filter(value < 0)

ampAGpos <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value > 0)

ampAGneg <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif=value) %>%
  filter(ag== "1") %>%
  filter(value < 0)

skewness(ampHGpos$value)
skewness(ampAGpos$value)

skewness(abs(ampHGneg$value))
skewness(abs(ampAGneg$value))



#### Climate median scatterplots ----

setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/9_Climate")
q <- ggplot(data=Directory, aes(y=(minvamp200), x=(precstab200)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "200 Year Precipitation Stability", y="200 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, 0.14))
#facet_wrap(~(ew))
q

jpeg("200prec_minvamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(minvamp100), x=(precstab100)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "100 Year Precipitation Stability", y="100 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("100prec_minvamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(minvamp50), x=(precstab50)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "50 Year Precipitation Stability", y="50 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("50prec_minvamp.jpeg", width=500, height=300)
q
dev.off()


###Temperature now
q <- ggplot(data=Directory, aes(y=(minvamp50), x=(log(tempstab50))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 50 Year Temperature Stability", y="50 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("50temp_minvamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(minvamp100), x=(log(tempstab100))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 100 Year Temperature Stability", y="100 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("100temp_minvamp.jpeg", width=500, height=300)
q
dev.off()



q <- ggplot(data=Directory, aes(y=(minvamp200), x=(log(tempstab200))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 200 Year Temperature Stability", y="200 Year Median SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("200temp_minvamp.jpeg", width=500, height=300)
q
dev.off()




### Mean Boxplots ----

dir.create("~/Dropbox/R/thesis_trial2/5_boxplots/")
setwd("~/Dropbox/R/thesis_trial2/5_boxplots/")

p <- ggplot(Directory, aes(x=(ag), abs(invamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("50 Year Mean SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("invamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()

p <- ggplot(Directory, aes(x=(ag), abs(invamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Mean SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("invamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(invamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Mean SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("invamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()


###Positive
p <- ggplot(Directory, aes(x=(ag), abs(posamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("50 Year Mean SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("posamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()



p <- ggplot(Directory, aes(x=(ag), abs(posamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Mean SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("posamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(posamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Mean SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("posamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()

####Negative 

p <- ggplot(Directory, aes(x=(ag), abs(negamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("50 Year Mean SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("negamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(negamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Mean SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("negamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()

p <- ggplot(Directory, aes(x=(ag), abs(negamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Mean SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
p+stat_compare_means(label.x=0.7, label.y=.0025)
p

jpeg("negamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()

### Median Boxplots -----
setwd("~/Dropbox/R/thesis_trial2/5_boxplots/")

p <- ggplot(Directory, aes(x=(ag), abs(minvamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  ylab("50 Year Median SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
      #axis.title.y=element_blank(),
      text = element_text(size=15))
 #axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("minvamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()

p <- ggplot(Directory, aes(x=(ag), abs(minvamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Median SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("minvamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(minvamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Median SPD Stability")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("minvamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()


###Positive
p <- ggplot(Directory, aes(x=(ag), abs(mposamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("50 Year Median SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mposamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()



p <- ggplot(Directory, aes(x=(ag), abs(mposamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Median SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mposamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(mposamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Median SPD Increase")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        #axis.title.y=element_blank(),
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mposamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()

####Negative 

p <- ggplot(Directory, aes(x=(ag), abs(mnegamp50)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("50 Year Median SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mmnegamp50_boxplot.jpeg", width=205, height=300)
p
dev.off()


p <- ggplot(Directory, aes(x=(ag), abs(mnegamp100)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("100 Year Median SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mnegamp100_boxplot.jpeg", width=205, height=300)
p
dev.off()

p <- ggplot(Directory, aes(x=(ag), abs(mnegamp200)))+
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill=cbbPalette, colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  ylab("200 Year Median SPD Decrease")+
  scale_x_discrete(labels=c("Forager", "Mixed Ag"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        text = element_text(size=15))
#  axis.title.y = element_text(size=14))
#p+stat_compare_Medians(label.x=0.7, label.y=.0025)
p

jpeg("mnegamp200_boxplot.jpeg", width=205, height=300)
p
dev.off()

##### Remove 6 borderlands ----
setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/3_FirstDiff/")

Dif50long <- read_csv ("Dif50long.csv")
Dif50longb <- Dif50long %>%
  rename(Sbox = variable) %>%
  select(Sbox, dates, value, ag) %>%
  filter(!Sbox %in% c("10", "25", "29", "30", "31", "37"))

mu<-plyr::ddply(Dif50longb, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif50longb, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif50longb, "ag", summarise, grp.sd=sd(abs(value)))

p <-ggplot(Dif50longb, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = agnames, values = cbbPalette)+
  labs(x="Absolute Value of 50 Year First Difference Values", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~.)+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

mu
me
sd
table(Dif50longb$ag)
Dif50longb_pos <- subset(Dif50longb, value >0) 
table(Dif50longb_pos$ag)

library(tidyverse)

#Skewness
ampHG <- Dif50longb %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "0") 

ampAG <- Dif50longb %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "1") 

skewness(abs(ampHG$value))
skewness(abs(ampAG$value))





##### Scatter Log N vs. invamp -----
dir.create("~/Dropbox/R/thesis_trial2/7_lognCheck/")
setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/7_lognCheck/")

q <- ggplot(data=Directory, aes(x=(minvamp200), y=log(n)))+
  theme_bw() +
  #aes(shape=factor(ag), colour=factor(ag))+
  #scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  #scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "200 Year Median SPD Stability", y="Logged Number of 14C dates")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, 8))
#facet_wrap(~(ew))
q

stab<-lm((invamp200)~log(n), data=Directory)
summary(stab)

jpeg("200minvamp_n.jpeg", width=500, height=300)
q
dev.off()


##### Histogram: N 14C DATES VS TIME ------
dir.create("~/Dropbox/R/thesis_trial2/8_N14CHist/")
setwd("~/Dropbox/R/thesis_trial2/8_N14CHist/")

by_sbox2 <- RawData %>%
  left_join(AgID, by="Sbox") 

by_sbox2 <- by_sbox2[!is.na(by_sbox2$ag), ]

p <-ggplot(by_sbox2, aes(x=(abs(date)), fill=factor(ag))) +
  geom_histogram(binwidth = 200, boundary=8000)+
  theme_bw() +
  scale_fill_manual(name="Number of 14C dates before present", labels = agnames, values = cbbPalette)+
  labs(x="14C BP", y = "Count")+
  theme(legend.position="top")+
  scale_x_reverse()+
  facet_wrap(~Sbox)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p

#jpeg("HistN14CSbox.jpeg", width=912, height=800)
#p
#dev.off()

pg <- ggplot_build(p)
pg2 <- as.data.frame(pg$data)

N14CSbox <- pg2 %>%
  left_join(y=pg$layout$layout, by= "PANEL") %>%
  select(Sbox, count, xmin)%>%
  #mutate(fill=replace(fill, fill=="#56B4E9", "0"),  
   #      fill=replace(fill, fill=="#D55E00", "1"))%>%
  #rename(ag=fill)%>%
  spread(key = Sbox, value = count)

write.csv(N14CSbox, "HistN14CSbox.csv", row.names = FALSE)

###Divided by ag
p <-ggplot(by_sbox2, aes(x=(abs(date)), fill=factor(ag))) +
  geom_histogram(binwidth = 200, boundary=8000,)+
  theme_bw() +
  scale_fill_manual(name="Number of 14C dates before present", labels = agnames, values = cbbPalette)+
  labs(x="14C BP", y = "Count")+
  theme(legend.position="top")+
  scale_x_reverse()+
  facet_grid(factor(ag)~.)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

#jpeg("HistN14CAg.jpeg", width=912, height=390)
#p
#dev.off()

pg <- ggplot_build(p)
pg2 <- as.data.frame(pg$data)

N14CAg <- pg2 %>%
  left_join(y=pg$layout$layout, by= "PANEL") %>%
  select("factor(ag)", count, xmin)%>%
  rename(ag="factor(ag)") %>%
  spread(key = ag, value = count) %>%
  rename(AG="1", HG="0")

write.csv(N14CAg, "HistN14CAg.csv", row.names = FALSE)

##### Histogram: N bins vs time -----

#Sometimes r treats values within a dataframe in a way you cannot use. These lines ensure our calibration will work.
RawData$date <- as.numeric(RawData$date)
RawData$sd <- as.numeric(RawData$sd)
RawData$labnumber <- as.character(RawData$labnumber)

Directory_small<- aggregate(data.frame(count = RawData$Sbox), list(value = RawData$Sbox), length) 
names(Directory_small) <- c("Sbox", "n")
Directory_small<-Directory_small[!(Directory_small$n<200),] #The following code works better with a simple Directory

##Turn each sampling unit into a data.frame and make a list of them.
SboxList <- list()
for(i in 1:length(unique(Directory_small$Sbox))){
  nam <- make.names(paste("Sbox", Directory_small[i,"Sbox"]))
  assign(nam, RawData[RawData$Sbox == Directory_small[i,"Sbox"],]) #This line makes a dataframe for each sampling unit
  SboxList[i] <- lapply(make.names(paste("Sbox",Directory_small[i,"Sbox"])), get)  #this makes a list of the sampling units.
}
remove(nam)
remove(i)

##Calibration function for each hemisphere according to the different calibration curves
north <- function(Sbox){
  Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
  print(length(unique(Sboxbins))) ###Manually add results to Add to summary csv.
  
}

#Run to print the number of bins for each box
for(i in 1:length(unique(Directory_small$Sbox))){
  north(data.frame(SboxList[i])) 
}



##### False Positive Model ------
dir.create("~/Dropbox/R/thesis_trial2/6_model/")
setwd("~/Dropbox/R/thesis_trial2/6_model/")

model_base <- Directory %>%
  select(Sbox, ag, invamp50, invamp100, invamp200)

df2<- replicate(1000, sample((0:1),39,replace=T))

df<- cbind(model_base, df2)

results50 <- lapply(df[6:1005],function(x) wilcox.test(invamp50~x, data=df, alternative = "two.sided"))
results50_1 <- do.call(cbind,lapply(results50,function(v){v$p.value}))
results50_2 <- as.data.frame(rbind(results50_1,df2))
results50_1 <- as.list(do.call(cbind,lapply(results50,function(v){v$p.value})))

results100 <- lapply(df[6:1005],function(x) wilcox.test(invamp100~x, data=df, alternative = "two.sided"))
results100_1 <- do.call(cbind,lapply(results100,function(v){v$p.value}))
results100_2 <- as.data.frame(rbind(results100_1,df2))
results100_1<- as.list(do.call(cbind,lapply(results100,function(v){v$p.value})))

results200 <- lapply(df[6:1005],function(x) wilcox.test(invamp200~x, data=df, alternative = "two.sided"))
results200_1 <- do.call(cbind,lapply(results200,function(v){v$p.value}))
results200_2 <- as.data.frame(rbind(results200_1,df2))
results200_1 <- as.list(do.call(cbind,lapply(results200,function(v){v$p.value})))

result_50 <- data.frame(matrix(nrow = 3, ncol = 2))
result_100 <- data.frame(matrix(nrow = 3, ncol = 2))
result_200 <- data.frame(matrix(nrow = 3, ncol = 2))

for (i in 1:1000){
  if(results200_1[i] <= 0.05){
    staty <- Yule(table(df2[,i],df[,2]))
    result_200[i, 1] <- i
    result_200[i, 2] <- staty
  }
}

result_200 <- na.omit(result_200)


for (i in 1:1000){
  if(results100_1[i] <= 0.05){
    staty <- Yule(table(df2[,i],df[,2]))
    result_100[i, 1] <- i
    result_100[i, 2] <- staty
  }
}

result_100 <- na.omit(result_100)


for (i in 1:1000){
  if(results50_1[i] <= 0.05){
    staty <- Yule(table(df2[,i],df[,2]))
    result_50[i, 1] <- i
    result_50[i, 2] <- staty
  }
}

result_50 <- na.omit(result_50)

write.csv(result_50, "result50.csv", row.names = FALSE)
write.csv(result_100, "result100.csv",  row.names = FALSE)
write.csv(result_200, "result200.csv",  row.names = FALSE)

##### Taylor's Law with Non-normalized Data ----
setwd("~/Dropbox/R/thesis_trial2")
AgID <- read_csv("AgID.csv")

setwd("~/Dropbox/R/thesis_trial2/nonnormalized/2_Sbox_Bins/")
Sum50 <- read_csv("SboxSum50.csv")
Sum100 <- read_csv("SboxSum100.csv")
Sum200 <- read_csv("SboxSum200.csv")

cbbPalette <- c("#56B4E9", "#D55E00")

### 50 year
norm_sums50 <- as.data.frame(cbind(colnames(Sum50),
                                   as.numeric(colMeans(abs(Sum50), na.rm = TRUE)),
                                   as.numeric(colMedians(as.matrix(abs(Sum50)), na.rm = TRUE)),
                                   as.numeric(colSds(as.matrix(Sum50), na.rm = TRUE))))
colnames(norm_sums50) <- c("Sbox", "SPDmean", "SPDmedian", "SPDsd")
norm_sums50$SPDmean <-as.numeric(as.character(norm_sums50$SPDmean))
norm_sums50$SPDmedian <- as.numeric(as.character(norm_sums50$SPDmedian))
norm_sums50$SPDsd <- as.numeric(as.character(norm_sums50$SPDsd))
norm_sums50$SPDvar <- (norm_sums50$SPDsd)^2
norm_sums50 <- merge(norm_sums50, AgID, by.x = "Sbox", by.y = "Sbox") 


q <- ggplot(data=norm_sums50, aes(x=log(SPDmean), y=log(SPDvar)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

q <- ggplot(data=norm_sums50, aes(x=log(SPDmedian), y=log(SPDvar)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

m <- lm(log(SPDvar)~log(SPDmean)*ag, data=norm_sums50)
summary(m)
norm_sums50$meanResiduals <- residuals(m)

md <- lm(log(SPDvar)~log(SPDmedian)*ag, data=norm_sums50)
summary(md)
plot(m)
norm_sums50$medianResiduals <- residuals(md)

ggplot(data=norm_sums50, aes(x=(meanResiduals)))+
  geom_histogram()
ggplot(data=norm_sums50, aes(x=(medianResiduals)))+
  geom_histogram()

shapiro.test(norm_sums50$meanResiduals)
shapiro.test(norm_sums50$medianResiduals)

norm_sums50$zmeanResiduals <- scores((norm_sums50$meanResiduals))
norm_sums50$zmedianResiduals <- scores((norm_sums50$medianResiduals))

norm_sums50$zmean2 <- as.factor(ifelse(abs(norm_sums50$zmeanResiduals) >2, 1,0))
norm_sums50$zmedian2 <- as.factor(ifelse(abs(norm_sums50$zmedianResiduals) >2, 1,0))
norm_sums50_zmean2 <- subset(norm_sums50, zmean2 == 1)
norm_sums50_zmedian2 <- subset(norm_sums50, zmedian2 == 1)

q <- ggplot(data=norm_sums50_zmean2, aes(y=log(SPDmean), x=log(SPDsd)^2))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

q <- ggplot(data=norm_sums50_zmedian2, aes(y=log(SPDmedian), x=log(SPDsd)^2))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

####200 year
norm_sums200 <- as.data.frame(cbind(colnames(Sum200),
                                    as.numeric(colMeans(abs(Sum200), na.rm = TRUE)),
                                    as.numeric(colMedians(as.matrix(abs(Sum200)), na.rm = TRUE)),
                                    as.numeric(colSds(as.matrix(Sum200), na.rm = TRUE))))
colnames(norm_sums200) <- c("Sbox", "SPDmean", "SPDmedian", "SPDsd")
norm_sums200$SPDmean <-as.numeric(as.character(norm_sums200$SPDmean))
norm_sums200$SPDmedian <- as.numeric(as.character(norm_sums200$SPDmedian))
norm_sums200$SPDsd <- as.numeric(as.character(norm_sums200$SPDsd))
norm_sums200 <- merge(norm_sums200, AgID, by.x = "Sbox", by.y = "Sbox") 

q <- ggplot(data=norm_sums200, aes(x=log(SPDmean), y=log(SPDsd)^2))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

q <- ggplot(data=norm_sums200, aes(x=log(SPDmedian), y=log(SPDsd)^2))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  #labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, .9))
#facet_wrap(~(ew))
q

m <- lm(log(SPDsd)~log(SPDmean)*ag, data=norm_sums200)
summary(m)
norm_sums200$meanResiduals <- residuals(m)

md <- lm(log(SPDsd)~log(SPDmedian)*ag, data=norm_sums200)
summary(md)
#plot(m)
norm_sums200$medianResiduals <- residuals(md)

ggplot(data=norm_sums200, aes(x=(meanResiduals)))+
  geom_histogram()
ggplot(data=norm_sums200, aes(x=(medianResiduals)))+
  geom_histogram()

shapiro.test(norm_sums200$meanResiduals)
shapiro.test(norm_sums200$medianResiduals)





##### Alternative Growth Rate -----

####50 year bins

Rawdata1 <- read.csv("Darcy50Long.csv", header=TRUE)

Rawdata<- subset(Rawdata1, NegID=="1")

library(plyr)
mu<-plyr::ddply(Rawdata, "ID", summarise, grp.mean=mean((grow)))
me<-plyr::ddply(Rawdata, "ID", summarise, grp.median=median((grow)))
sd<-plyr::ddply(Rawdata, "ID", summarise, grp.sd=sd((grow)))
mu
me

p <-ggplot(Rawdata, aes(x=(((grow))), fill=factor(ID))) +
  #scale_x_log10()+
  #geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_histogram(binwidth=0.01)+
  geom_rug(aes(x =((grow)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  labs(x="50 year growth rate", y = "Frequency of observation")+
  theme(legend.position="top")+
  facet_grid(factor(ID)~.)
p


p <-ggplot(Rawdata1, aes(x=(((grow))), fill=factor(ID))) +
  #geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_histogram(binwidth=0.01)+
  geom_rug(aes(x =((grow)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  labs(x="50 year growth rate", y = "Frequency of observation")+
  theme(legend.position="top")+
  facet_grid(factor(ID)~.)
p

Growth<- subset(Rawdata1, sbox=="X8")
Growth2<- subset(Rawdata1, sbox=="X16")

p <- ggplot(Growth, aes((spd), (grow)))
p +theme_bw() +
  geom_point(size=2)+ 
  theme(axis.text = element_text(size = rel(1.8), colour = "black"))+
  labs(x="SPD", y = "Growth rate")+
  geom_smooth(se=FALSE)
# facet_grid(factor(ID)~.)

p <- ggplot(Growth2, aes((spd), (grow)))
p +theme_bw() +
  geom_point(size=2)+ 
  theme(axis.text = element_text(size = rel(1.8), colour = "black"))+
  labs(x="SPD", y = "Growth rate")+
  geom_smooth(se=FALSE)
# facet_grid(factor(ID)~.)




######200 year bins

Rawdata1 <- read.csv("Darcy200Long.csv", header=TRUE)

Rawdata<- subset(Rawdata1, NegID=="1")

library(plyr)
mu<-plyr::ddply(Rawdata, "ID", summarise, grp.mean=mean((grow)))
me<-plyr::ddply(Rawdata, "ID", summarise, grp.median=median((grow)))
sd<-plyr::ddply(Rawdata, "ID", summarise, grp.sd=sd((grow)))
mu
me

p <-ggplot(Rawdata, aes(x=(((grow))), fill=factor(ID))) +
  #scale_x_log10()+
  #geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_histogram(binwidth=0.04)+
  geom_rug(aes(x =((grow)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  labs(x="200 year growth rate", y = "Frequency of observation")+
  theme(legend.position="top")+
  facet_grid(factor(ID)~.)
p


p <-ggplot(Rawdata1, aes(x=(((grow))), fill=factor(ID))) +
  #geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_histogram(binwidth=0.02)+
  geom_rug(aes(x =((grow)) , y = 0), position = position_jitter(height = 0))+
  #geom_vline(data=mu, aes(xintercept=grp.mean),
  # linetype="dashed", show.legend = FALSE)+
  #geom_vline(data=me, aes(xintercept=grp.median),
  # linetype="solid", show.legend = FALSE)+
  theme_bw() +
  labs(x="50 year growth rate", y = "Frequency of observation")+
  theme(legend.position="top")+
  facet_grid(factor(ID)~.)
p


Growth<- subset(Rawdata1, sbox=="X8")
Growth2<- subset(Rawdata1, sbox=="X16")

p <- ggplot(Growth, aes((spd), (grow)))
p +theme_bw() +
  geom_point(size=2)+ 
  theme(axis.text = element_text(size = rel(1.8), colour = "black"))+
  labs(x="SPD", y = "Growth rate")+
  geom_smooth(se=FALSE)
# facet_grid(factor(ID)~.)

p <- ggplot(Growth2, aes((spd), (grow)))
p +theme_bw() +
  geom_point(size=2)+ 
  theme(axis.text = element_text(size = rel(1.8), colour = "black"))+
  labs(x="SPD", y = "Growth rate")+
  geom_smooth(se=FALSE)
# facet_grid(factor(ID)~.)
