setwd() ##set your working directory to a folder with just “radiocarbon_data.csv” in it.
RawData <- read.csv(“radiocarbon_data.csv")

library(plyr)
library(searchable)
library(reshape)
library(rcarbon)
library(ggplot2)
library(rcarbon)
library(ggpubr)
library(zoo)
library(robustbase)
library(moments)
#library(tidyverse)
#remove.packages(tidyverse)

#### 2. Make Directory for all sample units with more than 200 lab numbers -------
Directory <- count(RawData, vars = "Sbox")
names(Directory) <- c("Sbox", "n")
Directory<-Directory[!(Directory$n<200),]
write.csv(Directory, "Directory_Sbox.csv", row.names = FALSE)

#### 3. Calibrate each sampling unit recursively --------

#Sometimes r treats values within a dataframe in a way you cannot use. The lines ensure our calibration will work.
RawData$date <- as.numeric(RawData$date)
RawData$sd <- as.numeric(RawData$sd)
RawData$labnumber <- as.character(RawData$labnumber)

##Turn each sampling unit into a data.frame and make a list of them. 
SboxList <- list()
for(i in 1:length(unique(Directory$Sbox))){
  nam <- make.names(paste("Sbox", Directory[i,"Sbox"]))
  assign(nam, RawData[RawData$Sbox == Directory[i,"Sbox"],]) #This line makes a dataframe for each sampling unit
  SboxList[i] <- lapply(make.names(paste("Sbox",Directory[i,"Sbox"])), get)  #this makes a list of the sampling units.
}
remove(nam)
remove(i)

cptcal <- calibrate(x = Sbox$date,  errors = Sbox$sd) #This calibrates the dates using the default intcal13
Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
Sboxspd <- spd(x=cptcal, timeRange=c(6000,300), spdnormalised = TRUE) #This produces normalized SPD values
write.csv(Sboxspd,file = paste("Sbox TX.csv")) #This writes the SPD values to the working directory, allowing you to view them outside of R and pull them back in later.



##Calibration function for each hemisphere according to the different calibration curves
north <- function(Sbox){
  cptcal <- calibrate(x = Sbox$date,  errors = Sbox$sd) #This calibrates the dates using the default intcal13
  Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
  Sboxspd <- spd(x=cptcal, timeRange=c(6000,300), spdnormalised = TRUE) #This produces normalized SPD values
  write.csv(Sboxspd,file = paste("Sbox", Directory[i,"Sbox"], ".csv")) #This writes the SPD values to the working directory, allowing you to view them outside of R and pull them back in later.
}

south <- function(Sbox){
  cptcal <- calibrate(x = Sbox$date,  errors = Sbox$sd,calCurves = 'shcal13') #This calibrates the dates using the default shcal13
  Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
  Sboxspd <- spd(x=cptcal, timeRange=c(6000,300), spdnormalised = TRUE) #This produces normalized SPD values.  
  write.csv(Norm,file = paste("Sbox", Directory[i,"Sbox"], ".csv")) #This writes the SPD values to the working directory, allowing you to view them outside of R and pull them back in later.
}

##Create a directory to store the SPD results
dir.create("~/R/THESIS_2/1_Sbox_SPD")
setwd("~/R/THESIS_2/1_Sbox_SPD")

#Run the code to calibrate recursively, This may take awhile.
for(i in 1:length(unique(Directory$Sbox))){
  if(nrow(data.frame(SboxList[i])) >= 200){
     if(RawData$Lat >= 0){
      north(data.frame(SboxList[i]))
    }
    else{
      (south(data.frame(SboxList[i])))
    }
  }
}

#Clean the environment.
rm(list=ls())


#### 4. Recursively bin all SPDs -----
setwd("~/R/THESIS_2/1_Sbox_SPD")

temp <- list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)

files <- list.files(path="./")
dates <- read.table(files[1], sep=",", header=TRUE)[,11]     # gene names
df    <- do.call(cbind,lapply(files,function(fn)read.table(fn, header=TRUE, sep=",")[,12]))
df2 <- gsub('Sbox', '', files)
df2 <- sub(' ', '', df2)
df2 <- sub(' .csv', '', df2)
colnames(df) <- df2
df3 <- cbind(dates,df)

dir.create("~/R/THESIS_2/2_Sbox_Bins/")
setwd("~/R/THESIS_2/2_Sbox_Bins/")

###Sum the spd data at different bin widths
out10 <- rollapply(df3,50,(sum),by=50,by.column=TRUE,align='right')
out10 <- as.data.frame(out10)
out10$dates <- ((out10$dates / 50) -25.5)
write.table(out10, file = "SboxSum50.csv", sep = ",", row.names = FALSE)

out20 <- rollapply(df3,100,(sum),by=100,by.column=TRUE,align='right')
out20 <-as.data.frame(out20)
out20$dates<-((out20$dates/100) - 50.5)
write.table(out20, file = "SboxSum100.csv", sep = ",", row.names = FALSE)

out50 <- rollapply(df3,200,(sum),by=200,by.column=TRUE,align='right')
out50<-as.data.frame(out50)
out50$dates<-((out50$dates/200)-100.5)
write.table(out50, file = "SboxSum200.csv", sep = ",", row.names = FALSE)


rm(list=ls())


#### 5. Calculate First Difference Values ----
setwd("~/R/THESIS_2/2_Sbox_Bins/")
Sum50 <- read.csv("SboxSum50.csv")
Sum100 <- read.csv("SboxSum100.csv")
Sum200 <- read.csv("SboxSum200.csv")

setwd("~/R/THESIS_2/")
Directory <- read.csv("Directory_Sbox.csv")

dir.create("~/R/THESIS_2/3_FirstDiff/")
setwd("~/R/THESIS_2/3_FirstDiff/")

###Calculate First Difference Values for each time series. Let's start with the 50 year time scale
dates<- Sum50$dates[-c(114)] ##extract the dates for the time series
Sum50 <- Sum50[-c(1)] ###Remove the dates
Sum50_2 <- Sum50[-c(1),] ###Remove the first row (year 5950 BP)
rownames(Sum50_2) <- 1:113 ### Renumber row names so we can properly subtract
Sum50 <- Sum50[-c(114),] ### Remove the last row (year 300BP)
SboxDif50 <- Sum50_2 - Sum50 ### younger SPD values (starting with 5900BP) - older SPD values, so positive numbers will demonstrate SPD increase, and negative will be decrease.
SboxDif50<- cbind(dates,SboxDif50) ###Recombine dates with the SPD difference values
SboxDif50[SboxDif50 == 0] <- NA ### All first difference values of "0" will be replaced with NA
remove(Sum50) #Clean up environment
remove(Sum50_2) #Clean up environment
remove(dates)
write.csv(SboxDif50, "SboxDif50.csv", row.names = FALSE) ##Write to a csv.

###Let's move on to the 100 year time scale. The only change will be how long the dataframe will be, which will affect the math.
dates<- Sum100$dates[-c(57)] ##extract the dates for the time series
Sum100 <- Sum100[-c(1)] ###Remove the dates
Sum100_2 <- Sum100[-c(1),] ###Remove the first row (year 5900 BP)
rownames(Sum100_2) <- 1:56 ### Renumber row names so we can properly subtract
Sum100 <- Sum100[-c(57),] ### Remove the last row (year 300BP)
SboxDif100 <- Sum100_2 - Sum100 ### younger SPD values (starting with 5900BP) - older SPD values, so positive numbers will demonstrate SPD increase, and negative will be decrease.
SboxDif100<- cbind(dates,SboxDif100) ###Recombine dates with the SPD difference values
SboxDif100[SboxDif100 == 0] <- NA ### All first difference values of "0" will be replaced with NA
remove(Sum100) #Clean up environment
remove(Sum100_2) #Clean up environment
remove(dates) #Clean up environment
write.csv(SboxDif100, "SboxDif100.csv", row.names = FALSE) ##Write to a csv.

###Finally, let's do the 200 year time scale.
dates<- Sum200$dates[-c(28)] ##extract the dates for the time series
Sum200 <- Sum200[-c(1)] ###Remove the dates
Sum200_2 <- Sum200[-c(1),] ###Remove the first row (bin 6000-5800 BP)
rownames(Sum200_2) <- 1:27 ### Renumber row names so we can properly subtract
Sum200 <- Sum200[-c(28),] ### Remove the last row (year bin 600-400BP)
SboxDif200 <- Sum200_2 - Sum200 ### younger SPD values (starting with 5600BP) - older SPD values, so positive numbers will demonstrate SPD increase, and negative will be decrease.
SboxDif200<- cbind(dates,SboxDif200) ###Recombine dates with the SPD difference values
SboxDif200[SboxDif200 == 0] <- NA ### All first difference values of "0" will be replaced with NA
remove(Sum200) #Clean up environment
remove(Sum200_2) #Clean up environment
remove(dates) #Clean up environment
write.csv(SboxDif200, "SboxDif200.csv", row.names = FALSE) ##Write to a csv.

#### 6. Calculate mean and median first difference values, then add to directory ----

##Let's start with 50 year time scale
NegDif50 <- as.data.frame(apply(SboxDif50,  MARGIN = c(1,2), function(x) {ifelse(x < 0, NA, x)})) ###Create a dataframe of all positive first difference values
PosDif50 <- as.data.frame(apply(SboxDif50,  MARGIN = c(1,2), function(x) {ifelse(x > 0, NA, x)})) ###Create a dataframe of all negative first difference values.

###This code calculates average mean amplitude, median amplitude, mean positive first difference values, median positive first difference values, mean negative first difference values, and median negative first difference values.
amp50<- as.data.frame(cbind(colnames(SboxDif50), 
                             colMeans(abs(SboxDif50), na.rm = TRUE), 
                             colMedians(as.matrix(abs(SboxDif50)), na.rm = TRUE), 
                             colMeans(PosDif50, na.rm = TRUE), 
                             colMedians(as.matrix(PosDif50), na.rm = TRUE), 
                             colMeans(NegDif50, na.rm = TRUE), 
                             colMedians(as.matrix(NegDif50), na.rm = TRUE))) 
colnames(amp50) <- c("Sbox", "amp50", "mamp50", "posamp50", "mposamp50", "negamp50", "mnegamp50") ###Rename column headings to clarify

amp50$amp50 <- as.numeric(as.character(amp50$amp50)) ##Convert the column from factors to non-integer numbers
amp50$mamp50 <- as.numeric(as.character(amp50$mamp50)) 
amp50$invamp50 <- sapply(amp50$amp50, FUN=function(x) 1/x) ##Take the inverse of the mean amplitude values to represent mean stability
amp50$minvamp50 <- sapply(amp50$mamp50, FUN=function(x) 1/x) ##Take the inverse of the median amplitude values
amp50 <- amp50[-c(1),] #Get rid of the dates row
amp50$Sbox <- sub('X', '', amp50$Sbox) 

Directory <- merge(Directory, amp50, by.x = "Sbox", by.y = "Sbox") ##Merge all of our stability measurements to the directory.

###And then move to the 100 year time scale
NegDif100 <- as.data.frame(apply(SboxDif100,  MARGIN = c(1,2), function(x) {ifelse(x < 0, NA, x)})) ###Create a dataframe of all positive first difference values
PosDif100 <- as.data.frame(apply(SboxDif100,  MARGIN = c(1,2), function(x) {ifelse(x > 0, NA, x)})) ###Create a dataframe of all negative first difference values.

###This code calculates average mean amplitude, median amplitude, mean positive first difference values, median positive first difference values, mean negative first difference values, and median negative first difference values.
amp100<- as.data.frame(cbind(colnames(SboxDif100), 
                             colMeans(abs(SboxDif100), na.rm = TRUE), 
                             colMedians(as.matrix(abs(SboxDif100)), na.rm = TRUE), 
                             colMeans(PosDif100, na.rm = TRUE), 
                             colMedians(as.matrix(PosDif100), na.rm = TRUE), 
                             colMeans(NegDif100, na.rm = TRUE), 
                             colMedians(as.matrix(NegDif100), na.rm = TRUE))) 
colnames(amp100) <- c("Sbox", "amp100", "mamp100", "posamp100", "mposamp100", "negamp100", "mnegamp100") ###Rename column headings to clarify

amp100$amp100 <- as.numeric(as.character(amp100$amp100)) ##Convert the column from factors to non-integer numbers
amp100$mamp100 <- as.numeric(as.character(amp100$mamp100)) 
amp100$invamp100 <- sapply(amp100$amp100, FUN=function(x) 1/x) ##Take the inverse of the mean amplitude values to represent mean stability
amp100$minvamp100 <- sapply(amp100$mamp100, FUN=function(x) 1/x) ##Take the inverse of the median amplitude values
amp100 <- amp100[-c(1),] #Get rid of the dates row
amp100$Sbox <- sub('X', '', amp100$Sbox) 

Directory <- merge(Directory, amp100, by.x = "Sbox", by.y = "Sbox") ##Merge all of our stability measurements to the directory.


###Finally 200 year time scale
NegDif200 <- as.data.frame(apply(SboxDif200,  MARGIN = c(1,2), function(x) {ifelse(x < 0, NA, x)})) ###Create a dataframe of all positive first difference values
PosDif200 <- as.data.frame(apply(SboxDif200,  MARGIN = c(1,2), function(x) {ifelse(x > 0, NA, x)})) ###Create a dataframe of all negative first difference values.

###This code calculates average mean amplitude, median amplitude, mean positive first difference values, median positive first difference values, mean negative first difference values, and median negative first difference values.
amp200<- as.data.frame(cbind(colnames(SboxDif200), 
                             colMeans(abs(SboxDif200), na.rm = TRUE), 
                             colMedians(as.matrix(abs(SboxDif200)), na.rm = TRUE), 
                             colMeans(PosDif200, na.rm = TRUE), 
                             colMedians(as.matrix(PosDif200), na.rm = TRUE), 
                             colMeans(NegDif200, na.rm = TRUE), 
                             colMedians(as.matrix(NegDif200), na.rm = TRUE))) 
colnames(amp200) <- c("Sbox", "amp200", "mamp200", "posamp200", "mposamp200", "negamp200", "mnegamp200") ###Rename column headings to clarify

amp200$amp200 <- as.numeric(as.character(amp200$amp200)) ##Convert the column from factors to non-integer numbers
amp200$mamp200 <- as.numeric(as.character(amp200$mamp200)) 
amp200$invamp200 <- sapply(amp200$amp200, FUN=function(x) 1/x) ##Take the inverse of the mean amplitude values to represent mean stability
amp200$minvamp200 <- sapply(amp200$mamp200, FUN=function(x) 1/x) ##Take the inverse of the median amplitude values
amp200 <- amp200[-c(1),] #Get rid of the dates row
amp200$Sbox <- sub('X', '', amp200$Sbox) 

Directory <- merge(Directory, amp200, by.x = "Sbox", by.y = "Sbox") ##Merge all of our stability measurements to the directory.


####Finally, let's export that directory to reference is later.
setwd("~/R/THESIS_2/")
write.csv(Directory, "Directory_Sbox.csv", row.names = FALSE)

rm(list=ls())

#### 7. Let's begin the analysis.  ----
##First, you will need to open your Directory CSV and add agriculture values (0,1). Then load it below.

setwd("~/R/THESIS_2/")
Directory <- read.csv("Directory_Sbox.csv")

cbbPalette <- c("#56B4E9", "#D55E00")
labels <- c("0" = "Hunter-Gatherers", "1" = "Agriculturalists")

###  7.1 Look at all SPDs at once. 50-year below, change 50 to 100 and 200 to see those instead. ----
setwd("~/R/THESIS_2/2_Sbox_Bins/")
Sum50 <- read.csv("SboxSum50.csv")
Sum100 <- read.csv("SboxSum100.csv")
Sum200 <- read.csv("SboxSum200.csv")

Sum50long <- melt.data.frame(Sum50, id=c("dates"))
Sum100long <- melt.data.frame(Sum100, id=c("dates"))
Sum200long <- melt.data.frame(Sum200, id=c("dates"))

Sum50long$variable <- gsub('X', 'Sbox ', Sum50long$variable)
Sum100long$variable <- gsub('X', 'Sbox ', Sum100long$variable)
Sum200long$variable <- gsub('X', 'Sbox ', Sum200long$variable)

dir.create("~/R/THESIS_2/4_SPDs/")
setwd("~/R/THESIS_2/4_SPDs/")

#50 year first

p.list = lapply(sort(unique(Sum50long$variable)), function(i) {
    ggplot(Sum50long[Sum50long$variable==i,], aes((dates), (value))) +
    geom_line(show.legend=FALSE) +
    theme_bw() +
    theme(axis.text = element_text(angle=45, size=12, colour = "black"), axis.title=element_text(size=18))+
    labs(x = "Cal years BP", y="Summed probability")+
    ggtitle(paste(i, "SPD"))+
    #geom_point(colour= ifelse(value < value, "red", "blue"))+
    #geom_hline(yintercept = mean(value))+
    scale_x_reverse(breaks = seq(500,6000,500))
    #geom_vline(aes(xintercept=650, colour= "red"), linetype="solid", show.legend = FALSE)
})

p.list ##View your SPDs

pdf("Sum50_200.pdf", width = 6, height = 3) ##Export them as a pdf.
p.list 
dev.off()

#Now let's do 100 year time scale.
p.list = lapply(sort(unique(Sum100long$variable)), function(i) {
  ggplot(Sum100long[Sum100long$variable==i,], aes((dates), (value))) +
    geom_line(show.legend=FALSE) +
    theme_bw() +
    theme(axis.text = element_text(angle=45, size=12, colour = "black"), axis.title=element_text(size=18))+
    labs(x = "Cal years BP", y="Summed probability")+
    ggtitle(paste(i, "SPD"))+
    #geom_point(colour= ifelse(value < value, "red", "blue"))+
    #geom_hline(yintercept = mean(value))+
    scale_x_reverse(breaks = seq(500,6000,500))
})
p.list

pdf("SboxSum100.pdf", width = 6, height = 3) ##Export them as a pdf.
p.list 
dev.off()


#Finally, 200 year time scale.
p.list = lapply(sort(unique(Sum200long$variable)), function(i) {
  ggplot(Sum200long[Sum200long$variable==i,], aes((dates), (value))) +
    geom_line(show.legend=FALSE) +
    theme_bw() +
    theme(axis.text = element_text(angle=45, size=12, colour = "black"), axis.title=element_text(size=18))+
    labs(x = "Cal years BP", y="Summed probability")+
    ggtitle(paste(i, "SPD"))+
    #geom_point(colour= ifelse(value < value, "red", "blue"))+
    #geom_hline(yintercept = mean(value))+
    scale_x_reverse(breaks = seq(500,6000,500))
})
p.list

pdf("SboxSum200.pdf", width = 6, height = 3) ##Export them as a pdf.
p.list 
dev.off()


###  7.2 Make histograms ----
setwd("~/R/THESIS_2/3_FirstDiff/")
SboxDif50 <- read.csv ("SboxDif50.csv")
SboxDif100 <- read.csv ("SboxDif100.csv")
SboxDif200 <- read.csv ("SboxDif200.csv")

colnames(SboxDif50) <- gsub(x=colnames(SboxDif50), pattern= "X", replacement = "")
colnames(SboxDif100) <- gsub(x=colnames(SboxDif100), pattern= "X", replacement = "")
colnames(SboxDif200) <- gsub(x=colnames(SboxDif200), pattern= "X", replacement = "")

Dif50long <- melt.data.frame(SboxDif50, id.vars = c("dates"))
Dif100long <- melt.data.frame(SboxDif100, id.vars = c("dates"))
Dif200long <- melt.data.frame(SboxDif200, id.vars = c("dates"))

Dif50long <- merge.data.frame(x=na.omit(Dif50long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")
Dif100long <- merge.data.frame(x=na.omit(Dif100long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")
Dif200long <- merge.data.frame(x=na.omit(Dif200long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")

write.csv(Dif50long, "Dif50long.csv", row.names = FALSE)
write.csv(Dif100long, "Dif100long.csv", row.names = FALSE)
write.csv(Dif200long, "Dif200long.csv", row.names = FALSE)

Dif50long <- read.csv ("Dif50long.csv")
Dif100long <- read.csv ("Dif100long.csv")
Dif200long <- read.csv ("Dif200long.csv")

##Table Organization out of the way, let's view the histograms! As always, 50 year scale first.

mu<-ddply(Dif50long, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif50long, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif50long, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif50long, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette)+
  labs(x="Absolute Value of 50 Year First Difference Values", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~.)
p

jpeg("Sbox50_invamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
mu<-ddply(Dif100long, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif100long, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif100long, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif100long, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Absolute Value of 100 Year First Difference Values", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~.)
p

jpeg("Sbox100_invamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last

mu<-ddply(Dif200long, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif200long, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif200long, "ag", summarise, grp.sd=sd(abs(value)))

p<-ggplot(Dif200long, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = c("Hunter-Gatherers", "Agriculturalists"), values = cbbPalette, guide=FALSE)+
  labs(x="Absolute Value of 200 Year First Difference Values", y = "Density")+
  theme(legend.position="top")+
  facet_grid(factor(ag)~.)
p

jpeg("Sbox200_invamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd


##### 7.21 Let's do positive population changes next
Dif50long_pos <- subset(Dif50long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif50long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif50long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif50long_pos, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox50_posamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
Dif100long_pos <- subset(Dif100long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif100long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif100long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif100long_pos, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox100_posamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last
Dif200long_pos <- subset(Dif200long, value >0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif200long_pos, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif200long_pos, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif200long_pos, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox200_posamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd


##### 7.22 Negative now.
Dif50long_neg <- subset(Dif50long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif50long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif50long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif50long_neg, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox50_negamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
Dif100long_neg <- subset(Dif100long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif100long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif100long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif100long_neg, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox100_negamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last
Dif200long_neg <- subset(Dif200long, value <0) ###Uncomment to view only population increases. Change > to < to view only decreases.

mu<-ddply(Dif200long_neg, "ag", summarise, grp.mean=mean(abs(value)))
me<-ddply(Dif200long_neg, "ag", summarise, grp.median=median(abs(value)))
sd<-ddply(Dif200long_neg, "ag", summarise, grp.sd=sd(abs(value)))

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
  facet_grid(factor(ag)~.)
p

jpeg("Sbox200_negamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd








###  7.3 Skewness WORK IN PROGRESS. This code currently requires the creation of a new cvs (AmpT200.csv) that has one column with only Agriculturalist first difference values (ampAG) and the other column with only hunter gatherer first difference values (ampHG). Whichever one has fewer values will have an error and the excess values in the other column must be deleted, so this code needs to be improved. ———
ampT <- read.csv("ampT200.csv")

ampHGpos <- subset(ampT, ampHG >0)
ampAGpos <- subset(ampT, ampAG >0)

ampHGneg <- subset(ampT, ampHG <0)
ampAGneg <- subset(ampT, ampAG <0)

skewness(abs(ampT$ampHG))
skewness(abs(ampT$ampAG))

skewness(ampHGpos$ampHG)
skewness(ampAGpos$ampAG)

skewness(abs(ampHGneg$ampHG))
skewness(abs(ampAGneg$ampAG))




###  7.4 Subsistence Strategy Boxplots. Change the second variable, title, and label to match desired Y axis ----

dir.create("~/R/THESIS_2/5_boxplots/")
setwd("~/R/THESIS_2/5_boxplots/")

options(scipen=999)

p <- ggplot(Directory, aes(factor(ag), abs(posamp50)))+
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(fill = "steelblue2", colour = "black")+
  #scale_y_continuous(trans="log", labels = comma_format(digits = 1))+
  #coord_flip() +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Agriculture Factor", y="50 Year Mean Population Increase")+
  scale_x_discrete(labels=c("Hunter-Gatherers", "Agriculturalists"))+
  #facet_wrap(~(ew))+
  theme_bw()+
  stat_compare_means(label.x=0.7, label.y=.0025)
p 

jpeg("posamp50_boxplot.jpeg", width=305, height=440)
p
dev.off()

###  7.5 Scatterplots with Climate. Change the second variable, title, and label to match desired Y axis ----

setwd("~/R/THESIS_2/")
Directory <- read.csv("Directory_Sbox.csv")

dir.create("~/R/THESIS_2/9_Climate")
setwd("~/R/THESIS_2/9_Climate")

q <- ggplot(data=Directory, aes(y=(invamp200), x=(precstab200)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) + 
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "200 Year Precipitation Stability", y="Population Stability")+
  geom_smooth()+
  theme(legend.position = c(0.14, 0.12))
#facet_wrap(~(ew))
q

jpeg("200prec_invamp.jpeg", width=656, height=440)
q
dev.off()

##### Scatter Log N vs. invamp -----
dir.create("~/R/THESIS_2/7_lognCheck/")
setwd("~/R/THESIS_2/7_lognCheck/")

q <- ggplot(data=Directory, aes(x=(invamp200), y=log(n)))+
  theme_bw() +
  #aes(shape=factor(ag), colour=factor(ag))+
  #scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  #scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) + 
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "200 Year Mean Population Stability", y="Logged Number of 14C dates")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, 8))
#facet_wrap(~(ew))
q

stab<-lm((minvamp200)~log(n), data=Directory)
summary(stab)

jpeg("200invamp_n.jpeg", width=656, height=440)
q
dev.off()


##### Model ------
dir.create("~/R/THESIS_2/6_SI/")
setwd("~/R/THESIS_2/6_model/")

library(psych)
library(gmodels)

df <- read.csv("model_base.csv")

df2<- replicate(1000, sample((0:1),40,replace=T))

df<- cbind(df, df2)

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





