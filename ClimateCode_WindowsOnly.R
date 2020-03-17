library(raster)
library(rgdal)

#This line of code prompts us to select a Directory; this reduces the need for making
#changes.
climDirectory <- vector(length = 2L)
climDirectory[1] <-
  choose.dir("", "Please select a folder with raw TEMPERATURE rasters")
climDirectory[2] <-
  choose.dir("", "Please select a folder with raw PRECIPITATION rasters.")
climType <- c("Temp", "Precip")

#Select the path of your "Boxes" shapefile
Boxes <- shapefile(choose.files("","Please select your \"Boxes\" shapefile"))

#This line of code prompts us to select the file to be updated ("Directory" csv)
DirectoryFile <- choose.files(caption = "Please select your \"Directory\" file.")
Directory <- read.table(DirectoryFile, header = T, sep = ",")

#This code was designed to work with different shapefiles that had different
#names for ID columns; we assign a new column to standardize these names
Boxes$ID <- (1:length(Boxes[,1]))


for(iteration in 1:length(climDirectory)){
  
  climateFiles1 <- climDirectory[iteration]
  
  ifelse(iteration == 1, adjustment <- 274.15, adjustment <- 0)

  #We build a list of all files in the Directory with climate data
  climateFiles <-
    list.files(path = climateFiles1,
               "grid_data_.*",
               no.. = T)
  #and a list of the years represented; here, we rely on the files using the
  #PaleoView naming scheme, and that no extra files will be present
  sortedList <- gsub("grid_data_", "", climateFiles)
  sortedList <- gsub("BP.asc", "", sortedList)
  sortedList <- sort(as.numeric(sortedList), decreasing = F)
  
  
  #####invCV for 50 Year Timeseries####
  #Our climate data are ten year averages, so to get fifty years, we use 5; a
  #sequence of intervals of five will help. 
  fiftyInterval <-
    seq(1, (length(sortedList)-1), by = 5)
  
  for (i in 1:(length(fiftyInterval))) {
    
    #Load the rasters, aggregating each into 5 degree grid cells
    cat(paste("Loading Rasters for 50 YR ",sortedList[fiftyInterval[i]]," BP ", climType[iteration],"\n",sep = "" ))
    r1 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[fiftyInterval[i]],
          "BP.asc",
          sep = ""
          )) +
          adjustment, 
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r2 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[fiftyInterval[i]] + 10,
          "BP.asc",
          sep = ""
          )) + 
          adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
        
    r3 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[fiftyInterval[i]] + 20,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r4 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[fiftyInterval[i]] + 30,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r5 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[fiftyInterval[i]] + 40,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    
    #Ensure that the rasters are in WGS84
    proj4string(r1) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r2) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r3) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r4) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r5) <- CRS("+proj=longlat +datum=WGS84")
    
    #invCV is mean over standard deviation, so we need Mean and SD.
    cat("calculating mean...\n")
    meanRaster <- mean(r1, r2, r3, r4, r5)
    cat("calculating standard deviation...\n")
    sdRaster <-
      sqrt((((r1 - meanRaster) ^ 2) + ((r2 - meanRaster) ^ 2) + ((r3 - meanRaster) ^
                                                                   2) + ((r4 - meanRaster) ^ 2) + ((r5 - meanRaster) ^ 2)) / 5)
    
    #Calculate the Coefficient of Variation.
    cat("calculating stability...\n")
    invCVRaster <- meanRaster / sdRaster
    
    #Extract bounded rasters to boxes
    cat("extracting...\n")
    BoxesExtracted <- extract(invCVRaster, Boxes, fun = mean)
    
    #If our storage dataframe does not exist, create one.
    if(!exists("BoxesExtracted1")){BoxesExtracted1 <- as.data.frame(Boxes$ID)}
    
    #add a new column to the dataframe for each year extracted
    BoxesExtracted1 <- cbind(BoxesExtracted1, BoxesExtracted)
    
    #This line tracks column names
    if(!exists("columnNames")){columnNames <- "ID"}
    columnNames <- rbind(columnNames,as.character(paste("50YR_",sortedList[fiftyInterval[i]],"BP_",climType[iteration],"\n",sep = "")))
    colnames(BoxesExtracted1) <- columnNames
    #This line reports progress
    cat(paste("Completed Extraction: 50 year increments (",climType[iteration],"); ",sortedList[fiftyInterval[i]],"BP\n\n", sep = ""))
    
    
  }
  
  #####invCV for 100 Year Timeseries####
  hundredInterval <- seq(1, (length(sortedList)-1), by = 10)
  for (i in 1:(length(hundredInterval))) {
  
    #load and aggregate rasters into 5 degree grid
    cat(paste("Loading Rasters for 100 YR ",sortedList[hundredInterval[i]]," BP ", climType[iteration],"\n",sep = ""))
    r1 <-
      aggregate(
        raster(paste(climateFiles1, "\\grid_data_", sortedList[hundredInterval[i]], "BP.asc", sep = "")) +
        adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r2 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 10,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r3 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 20,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r4 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 30,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r5 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 40,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r6 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 50,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r7 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 60,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r8 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 70,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r9 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 80,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r10 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[hundredInterval[i]] + 90,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    
    proj4string(r1) <- CRS("+proj=longlat +datu=WGS84")
    proj4string(r2) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r3) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r4) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r5) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r6) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r7) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r8) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r9) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r10) <- CRS("+proj=longlat +datum=WGS84")
    
    
    #We calculated the mean and standard deviation long-hand in 100 and 200 year
    #data after it appeared that mean() and sd() functions had processing
    #limitations with raster-sized datasets; it is uncertain whether any such
    #limitation actually exists, but we kept the long-hand calculations to limit
    #that uncertainty.
    
    #mean = sum/n
    cat("calculating mean...\n")
    meanRaster <- r1 + r2
    meanRaster <- meanRaster + r3
    meanRaster <- meanRaster + r4
    meanRaster <- meanRaster + r5
    meanRaster <- meanRaster + r6
    meanRaster <- meanRaster + r7
    meanRaster <- meanRaster + r8
    meanRaster <- meanRaster + r9
    meanRaster <- meanRaster + r10
    meanRaster <- meanRaster / 10
    
    #sd = sqrt( sum( (x_i - mean)^2 for each x in i= 1 to n ) /n)
    cat("calculating standard deviation...\n")
    sumDev <- ((r1 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r2 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r3 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r4 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r5 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r6 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r7 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r8 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r9 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r10 - meanRaster) ^ 2)
    sdRaster <- sqrt(sumDev / 10)
    
    cat("calculating stability...\n")
    invCVRaster <- meanRaster / sdRaster
    
    #Extract bounded rasters to boxes
    cat("extracting...\n")
    BoxesExtracted <- extract(invCVRaster, Boxes, fun = mean)
    
    #If our storage dataframe does not exist, create one.
    if(!exists("BoxesExtracted1")){BoxesExtracted1 <- as.data.frame(Boxes$ID)}
    
    #add a new column to the dataframe for each year extracted
    BoxesExtracted1 <- cbind(BoxesExtracted1, BoxesExtracted)
    
    #This line tracks column names
    if(!exists("columnNames")){columnNames <- "ID"}
    columnNames <- rbind(columnNames,as.character(paste("100YR_",sortedList[hundredInterval[i]],"BP_",climType[iteration],sep = "")))
    colnames(BoxesExtracted1) <- columnNames
    #This line reports progress
    cat(paste("Completed Extraction: 100 year increments (",climType[iteration],"); ",sortedList[hundredInterval[i]],"BP\n\n", sep = ""))
    
  }
  
  
  #####invCV for 200 Year Timeseries####
  twoInterval <- seq(11, (length(sortedList)), by = 20)
  
  for (i in 1:(length(twoInterval)-1)) {
    
    #load and aggregate rasters into 5 degree grid
    cat(paste("Loading Rasters for 200 YR ",sortedList[twoInterval[i]]," BP ", climType[iteration],"\n",sep = ""))
    
    r1 <-
      aggregate(
        raster(paste(climateFiles1, "\\grid_data_", sortedList[twoInterval[i]], "BP.asc", sep = "")) +
        adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r2 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 10,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r3 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 20,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r4 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 30,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r5 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 40,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r6 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 50,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r7 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 60,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r8 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 70,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r9 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 80,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r10 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 90,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r11 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 100,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r12 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 110,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r13 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 120,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r14 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 130,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r15 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 140,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r16 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 150,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r17 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 160,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r18 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 170,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r19 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 180,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    r20 <-
      aggregate(
        raster(paste(
          climateFiles1,
          "\\grid_data_",
          sortedList[twoInterval[i]] + 190,
          "BP.asc",
          sep = ""
        )) + adjustment,
        fact = 2,
        fun = mean,
        na.rm = TRUE
      )
    
    
    proj4string(r1) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r2) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r3) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r4) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r5) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r6) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r7) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r8) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r9) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r10) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r11) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r12) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r13) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r14) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r15) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r16) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r17) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r18) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r19) <- CRS("+proj=longlat +datum=WGS84")
    proj4string(r20) <- CRS("+proj=longlat +datum=WGS84")
    
    #meanRaster <- mean(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
    cat("calculating mean...\n")
    meanRaster <- r1 + r2
    meanRaster <- meanRaster + r3
    meanRaster <- meanRaster + r4
    meanRaster <- meanRaster + r5
    meanRaster <- meanRaster + r6
    meanRaster <- meanRaster + r7
    meanRaster <- meanRaster + r8
    meanRaster <- meanRaster + r9
    meanRaster <- meanRaster + r10
    meanRaster <- meanRaster + r11
    meanRaster <- meanRaster + r12
    meanRaster <- meanRaster + r13
    meanRaster <- meanRaster + r14
    meanRaster <- meanRaster + r15
    meanRaster <- meanRaster + r16
    meanRaster <- meanRaster + r17
    meanRaster <- meanRaster + r18
    meanRaster <- meanRaster + r19
    meanRaster <- meanRaster + r20
    meanRaster <- meanRaster / 20
    #sdRaster <- sqrt((((r1-meanRaster)^2)+((r2-meanRaster)^2)+((r3-meanRaster)^2)+((r4-meanRaster)^2)+((r5-meanRaster)^2)+((r6-meanRaster)^2)+((r7-meanRaster)^2)+((r8-meanRaster)^2)+((r9-meanRaster)^2)+((r10-meanRaster)^2)+((r11-meanRaster)^2)+((r12-meanRaster)^2)+((r13-meanRaster)^2)+((r14-meanRaster)^2)+((r15-meanRaster)^2)+((r16-meanRaster)^2)+((r17-meanRaster)^2)+((r18-meanRaster)^2)+((r19-meanRaster)^2)+((r20-meanRaster)^2))/20)
    cat("calculating standard deviation...\n")
    sumDev <- ((r1 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r2 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r3 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r4 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r5 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r6 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r7 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r8 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r9 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r10 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r11 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r12 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r13 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r14 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r15 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r16 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r17 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r18 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r19 - meanRaster) ^ 2)
    sumDev <- sumDev + ((r20 - meanRaster) ^ 2)
    sdRaster <- sqrt(sumDev / 20)
    
    cat("calculating stability...\n")
    invCVRaster <- meanRaster / sdRaster
    
    #Extract bounded rasters to boxes
    cat("extracting...\n")
    BoxesExtracted <- extract(invCVRaster, Boxes, fun = mean)
    
    #If our storage dataframe does not exist, create one.
    if(!exists("BoxesExtracted1")){BoxesExtracted1 <- as.data.frame(Boxes$ID)}
    
    #add a new column to the dataframe for each year extracted
    BoxesExtracted1 <- cbind(BoxesExtracted1, BoxesExtracted)
    
    #This line tracks column names
    if(!exists("columnNames")){columnNames <- "ID"}
    columnNames <- rbind(columnNames,as.character(paste("200YR_",sortedList[twoInterval[i]],"BP_",climType[iteration],sep = "")))
    colnames(BoxesExtracted1) <- columnNames
    #This line reports progress
    cat(paste("Completed Extraction: 200 year increments (",climType[iteration],"); ",sortedList[twoInterval[i]],"BP\n\n", sep = ""))
  
  }
  

}

write.csv(BoxesExtracted1, file = "BoxesExtracted1.csv", row.names = F)


setwd("C:/Users/Darcy/Dropbox/R/thesis_trial2")
BoxesExtracted1 <- read_csv("BoxesExtracted1.csv")
Directory <- read_csv("Directory_Sbox.csv")

BoxesExtracted50Temp <- BoxesExtracted1[,grepl( "^50YR.*Temp*", names( BoxesExtracted1 ) )] #extract all columns with 50YR preceding and Temp somewhere in there
BoxesExtracted100Temp <- BoxesExtracted1[,grepl( "^100YR.*Temp*", names( BoxesExtracted1 ) )] 
BoxesExtracted200Temp <- BoxesExtracted1[,grepl( "^200YR.*Temp*", names( BoxesExtracted1 ) )] 

BoxesExtracted50Precip <- BoxesExtracted1[,grepl( "^50YR.*Precip*", names( BoxesExtracted1 ) )] 
BoxesExtracted100Precip <- BoxesExtracted1[,grepl( "^100YR.*Precip*", names( BoxesExtracted1 ) )] 
BoxesExtracted200Precip <- BoxesExtracted1[,grepl( "^200YR.*Precip*", names( BoxesExtracted1 ) )] 

ClimData <- as.data.frame(cbind(BoxesExtracted1$ID,
                                rowMeans(as.matrix(BoxesExtracted50Temp), na.rm = TRUE), #Calculate means of all year stability rates
                                rowMeans(as.matrix(BoxesExtracted100Temp), na.rm = TRUE),
                                rowMeans(as.matrix(BoxesExtracted200Temp), na.rm = TRUE),
                                rowMeans(as.matrix(BoxesExtracted50Precip), na.rm = TRUE),
                                rowMeans(as.matrix(BoxesExtracted100Precip), na.rm = TRUE),
                                rowMeans(as.matrix(BoxesExtracted200Precip), na.rm = TRUE)))

colnames(ClimData) <- c("Sbox", "tempstab50", "tempstab100", "tempstab200", "precstab50", "precstab100", "precstab200") ###Rename column headings to clarify




#Merge the new data to the Directory
cat("merging to Directory...")
Directory <- merge(Directory,ClimData, by.x = "Sbox", by.y = "Sbox")
cat("Saving final Directory file...")
#This script overwrites your original "Directory" file. You can change this by
#uncommenting the line below.
write.csv(Directory, file = "Directory_Sbox_Clim.csv", row.names = F)
