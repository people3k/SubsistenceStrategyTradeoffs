library(readr)
library(plyr)
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
library(tidyverse)

setwd("~/Dropbox/R/thesis_trial2")
RawData <- read_csv("Bird2020_Holocene.csv")
AgID <- read_csv("AgID.csv")

#### 2. Make Directory for all sample units with more than 200 lab numbers -------
#setwd("~/Dropbox/R/Thesis_trial2/")

Directory <- aggregate(data.frame(count = RawData$Sbox), list(value = RawData$Sbox), length)
names(Directory) <- c("Sbox", "n")
Directory<-Directory[!(Directory$n<200),]
Directory <- merge(Directory, AgID, by.x = "Sbox", by.y = "Sbox") 
write.csv(Directory, "Directory_Sbox.csv", row.names = FALSE)

#### 3. Calibrate each sampling unit recursively --------

#Sometimes r treats values within a dataframe in a way you cannot use. These lines ensure our calibration will work.
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


##Calibration function for each hemisphere according to the different calibration curves
north <- function(Sbox){
  cptcal <- calibrate(x = Sbox$date,  errors = Sbox$sd) #This calibrates the dates using the default intcal13
  Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
  Sboxspd <- spd(x=cptcal, bins=Sboxbins, timeRange=c(6000,300), spdnormalised = FALSE) #This produces normalized SPD values
  write.csv(Sboxspd$grid,file = paste("Sbox", Directory[i,"Sbox"], ".csv")) #This writes the SPD values to the working directory, allowing you to view them outside of R and pull them back in later.
}

south <- function(Sbox){
  cptcal <- calibrate(x = Sbox$date,  errors = Sbox$sd,calCurves = 'shcal13') #This calibrates the dates using the default shcal13
  Sboxbins <- binPrep(sites = Sbox$SiteID, ages = Sbox$date, h = 100) #This bins the values.
  Sboxspd <- spd(x=cptcal,bins=Sboxbins,  timeRange=c(6000,300), spdnormalised = TRUE) #This produces normalized SPD values.
  write.csv(Sboxspd$grid, file = paste("Sbox", Directory[i,"Sbox"], ".csv")) #This writes the SPD values to the working directory, allowing you to view them outside of R and pull them back in later.
}

##Create a directory to store the SPD results
dir.create("~/Dropbox/R/thesis_trial2/1_Sbox_SPD")
setwd("~/Dropbox/R/thesis_trial2/1_Sbox_SPD")

#Run the code to calibrate recursively, This may take awhile.
for(i in 1:length(unique(Directory$Sbox))){
  #if(nrow(data.frame(SboxList[i])) >= 200){
   #  if(RawData$Lat >= 0){
      north(data.frame(SboxList[i]))
    #}
    #else{
     # (south(data.frame(SboxList[i])))
  #  }
#  }
}

#Clean the environment.
rm(list=ls())


#### 4. Recursively bin all SPDs -----
setwd("~/Dropbox/R/thesis_trial2/1_Sbox_SPD")

temp <- list.files(pattern="*.csv") ##Select all files in the working directory
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), #Remove the .csv at the end of files
         read_csv), envir = .GlobalEnv)

files <- list.files(path="./")
dates <- read.table(files[1], sep=",", header=TRUE)[,2]     # Extract the dates column
df    <- do.call(cbind,lapply(files,function(fn)read.table(fn, header=TRUE, sep=",")[,3])) #Extract SPD value
df2 <- gsub('Sbox', '', files) #rename
df2 <- sub(' ', '', df2) #rename
df2 <- sub(' .csv', '', df2) #rename
colnames(df) <- df2 #headers
df3 <- cbind(dates,df) #merge dates with SPD values to create one table

dir.create("~/Dropbox/R/thesis_trial2/2_Sbox_Bins/")
setwd("~/Dropbox/R/thesis_trial2/2_Sbox_Bins/")

write.csv(df3, file="SboxSPD.csv", sep = ",", row.names=FALSE)###Store the new df (bin = 1)

###Sum the spd data at different bin widths
out10 <- rollapply(df3,50,(sum),by=50,by.column=TRUE,align='right') #create new dataset that sums all SPD values every 50 years
out10 <- as.data.frame(out10)
out10$dates <- ((out10$dates / 50) -25.5) #Fix the dates column
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
setwd("~/Dropbox/R/thesis_trial2/2_Sbox_Bins/")
Sum50 <- read_csv("SboxSum50.csv")
Sum100 <- read_csv("SboxSum100.csv")
Sum200 <- read_csv("SboxSum200.csv")

setwd("~/Dropbox/R/thesis_trial2/")
Directory <- read_csv("Directory_Sbox.csv")

dir.create("~/Dropbox/R/thesis_trial2/3_FirstDiff/")
setwd("~/Dropbox/R/thesis_trial2/3_FirstDiff/")

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
setwd("~/Dropbox/R/thesis_trial2/")
write.csv(Directory, "Directory_Sbox.csv", row.names = FALSE)

rm(list=ls())

#### 7. Extract and calculate Climate Stability values ----

## The climate stability code currently takes a lot of time and only works on a windows computer. 
## Because of this, we've made it a separate file on github. 
## You only need it for section 8.5, so 8-8.4 should run fine without it.

#### 8. Let's begin the analysis.  ----
##First, you will need to open your Directory CSV and add agriculture values (0,1). Then load it below.

setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/")
Directory <- read_csv("Directory_Sbox.csv")

Directory$ag <-as.factor(Directory$ag)

cbbPalette <- c("#56B4E9", "#D55E00") #Ensures everything is always color-coded the same

agnames <- list( #Necessary for facet labels (note that these are depreciated and may fail soon.)
  "0" = "Hunter-Gatherers",
  "1" = "Agriculturalists"
)
variable_labeller <- function(variable,value){
  return(agnames[value])
}

###  8.1 Look at all SPDs at once. ----
setwd("~/Dropbox/R/thesis_trial2/2_Sbox_Bins/")
Sum50 <- read_csv("SboxSum50.csv")
Sum100 <- read_csv("SboxSum100.csv")
Sum200 <- read_csv("SboxSum200.csv")

Sum50 <- as.data.frame(Sum50)
Sum100 <- as.data.frame(Sum100)
Sum200 <- as.data.frame(Sum200)

Sum50long <- melt.data.frame(Sum50, id=c("dates"))
Sum100long <- melt.data.frame(Sum100, id=c("dates"))
Sum200long <- melt.data.frame(Sum200, id=c("dates"))

Sum50long$variable <- gsub('X', 'Sbox ', Sum50long$variable)
Sum100long$variable <- gsub('X', 'Sbox ', Sum100long$variable)
Sum200long$variable <- gsub('X', 'Sbox ', Sum200long$variable)

dir.create("~/Dropbox/R/thesis_trial2/4_SPDs/")
setwd("~/Dropbox/R/thesis_trial2/4_SPDs/")

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


###  8.2 Make histograms ----
setwd("~/Dropbox/R/thesis_trial2/3_FirstDiff/")
SboxDif50 <- read_csv ("SboxDif50.csv")
SboxDif100 <- read_csv ("SboxDif100.csv")
SboxDif200 <- read_csv ("SboxDif200.csv")

colnames(SboxDif50) <- gsub(x=colnames(SboxDif50), pattern= "X", replacement = "")
colnames(SboxDif100) <- gsub(x=colnames(SboxDif100), pattern= "X", replacement = "")
colnames(SboxDif200) <- gsub(x=colnames(SboxDif200), pattern= "X", replacement = "")

SboxDif50 <- as.data.frame(SboxDif50)
SboxDif100 <- as.data.frame(SboxDif100)
SboxDif200 <- as.data.frame(SboxDif200)

Dif50long <- melt.data.frame(SboxDif50, id.vars = c("dates"))
Dif100long <- melt.data.frame(SboxDif100, id.vars = c("dates"))
Dif200long <- melt.data.frame(SboxDif200, id.vars = c("dates"))

Dif50long <- merge.data.frame(x=na.omit(Dif50long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")
Dif100long <- merge.data.frame(x=na.omit(Dif100long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")
Dif200long <- merge.data.frame(x=na.omit(Dif200long), y=Directory[,c("Sbox", "ag")], by.x= "variable", by.y = "Sbox")

write.csv(Dif50long, "Dif50long.csv", row.names = FALSE)
write.csv(Dif100long, "Dif100long.csv", row.names = FALSE)
write.csv(Dif200long, "Dif200long.csv", row.names = FALSE)

Dif50long <- read_csv ("Dif50long.csv")
Dif100long <- read_csv ("Dif100long.csv")
Dif200long <- read_csv ("Dif200long.csv")



##Table Organization out of the way, let's view the histograms! As always, 50 year scale first.

mu<-plyr::ddply(Dif50long, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif50long, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif50long, "ag", summarise, grp.sd=sd(abs(value)))


p <-ggplot(Dif50long, aes(x=(abs(value)), fill=factor(ag))) +
  geom_density(adjust=1, alpha=1,  position = "identity")+
  geom_rug(aes(x =(abs(value)) , y = 0), position = position_jitter(height = 0))+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed", show.legend = FALSE)+
  geom_vline(data=me, aes(xintercept=grp.median),
             linetype="solid", show.legend = FALSE)+
  theme_bw() +
  scale_fill_manual(name="First Difference", labels = agnames, values = cbbPalette)+
  labs(x="Absolute Value of 50 Year First Difference Values", y = "Density")+
  theme(legend.position="top",
              text = element_text(size=20))+
  #facet_grid(factor(ag)~.)+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox50_invamp.jpeg", width=912, height=390)
p
dev.off()

mu
me
sd

###100 year next
mu<-plyr::ddply(Dif100long, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif100long, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif100long, "ag", summarise, grp.sd=sd(abs(value)))

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
  theme(text = element_text(size=20))+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox100_invamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###200 year last

mu<-plyr::ddply(Dif200long, "ag", summarise, grp.mean=mean(abs(value)))
me<-plyr::ddply(Dif200long, "ag", summarise, grp.median=median(abs(value)))
sd<-plyr::ddply(Dif200long, "ag", summarise, grp.sd=sd(abs(value)))

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
  theme(text = element_text(size=20))+
  facet_grid(factor(ag)~., labeller=variable_labeller)
p

jpeg("Sbox200_invamp.jpeg", width=912, height=345)
p
dev.off()

mu
me
sd

###  8.3 Skewness ----

setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/3_FirstDiff/")
library(tidyverse) #Writes over 8.2's library(plyr)

Dif50long <- read_csv("Dif50long.csv")
Dif100long <- read_csv("Dif100long.csv")
Dif200long <- read_csv("Dif200long.csv")

#Raw 50
ampHG <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "0") 

ampAG <- Dif50long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "1") 

skewness(abs(ampHG$value))
skewness(abs(ampAG$value))


#Raw 100
ampHG <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "0") 

ampAG <- Dif100long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "1") 

skewness(abs(ampHG$value))
skewness(abs(ampAG$value))

#Raw 200
ampHG <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "0") 

ampAG <- Dif200long %>%
  select(ag, value)%>%
  #rename(dif = value) %>%
  filter(ag== "1") 

skewness(abs(ampHG$value))
skewness(abs(ampAG$value))





###  8.4 Scatterplots with Climate ----

setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/")
Directory <- read_csv("Directory_Sbox.csv")

dir.create("C:/Users/Darcy/Dropbox/R/thesis_trial2/9_Climate")
setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2/9_Climate")

q <- ggplot(data=Directory, aes(y=(invamp200), x=(precstab200)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "200 Year Precipitation Stability", y="200 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = c(0.14, 0.14))
#facet_wrap(~(ew))
q

jpeg("200prec_invamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(invamp100), x=(precstab100)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "100 Year Precipitation Stability", y="100 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("100prec_invamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(invamp50), x=(precstab50)))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "50 Year Precipitation Stability", y="50 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("50prec_invamp.jpeg", width=500, height=300)
q
dev.off()


###Temperature now
q <- ggplot(data=Directory, aes(y=(invamp50), x=(log(tempstab50))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 50 Year Temperature Stability", y="50 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("50temp_invamp.jpeg", width=500, height=300)
q
dev.off()


q <- ggplot(data=Directory, aes(y=(invamp100), x=(log(tempstab100))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 100 Year Temperature Stability", y="100 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
  #facet_wrap(~(ew))
q

jpeg("100temp_invamp.jpeg", width=500, height=300)
q
dev.off()



q <- ggplot(data=Directory, aes(y=(invamp200), x=(log(tempstab200))))+
  theme_bw() +
  aes(shape=factor(ag), colour=factor(ag))+
  scale_colour_manual("", labels=c("Hunter-Gatherer","Agriculturalists" ), values= cbbPalette)+
  scale_shape_manual("", labels=c( "Hunter-Gatherer", "Agriculturalists"), values= c(16,17))+
  geom_point(size=2.5) +
  theme(axis.text = element_text(size = rel(1.9), colour = "black"), axis.title=element_text(size=16))+
  labs(x = "Logged 200 Year Temperature Stability", y="200 Year Mean SPD Stability")+
  geom_smooth(method="lm")+
  theme(legend.position = "none")
#facet_wrap(~(ew))
q

jpeg("200temp_invamp.jpeg", width=500, height=300)
q
dev.off()
