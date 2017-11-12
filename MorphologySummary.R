author <- "Haley E. Speed, Ph.D."
analysisDate <- date()
project <- "Project"
info <- data.frame(Author = author, Date = analysisDate, Project = project)
print(info)

library(nlme)
library (dplyr)
library(rcompanion)

## FUNCTION DEFINITIONS

# Get directory and file information 
getUserInput <- function (){
        file.inDir <- readline("Enter the path to the working directory: ")
        file.inName1 <- readline("Enter the filename of group 1: ")
        file.inName2 <- readline("Enter the filename of group 2: ")
        file.inName3 <- readline("Enter the filename of group 3: ")
        file.inName4 <- readline("Enter the filename of group 4: ")
        return (c(file.inName1, file.inName2, file.inName3, file.inName4, 
                  file.inDir))
}

getDataTables <- function (file.name,file.dir) {
        # Constructs path to each file
        file.path <- paste(file.dir,file.name, sep = "\\")
        # Reads in the files into data frames        
        data.group <- read.csv(file.path, header = TRUE)
        # Gets group name from csv filename
        file.groupName <- strsplit(file.name, "_summary.csv")
        # Adds a column named "group" to the existing data set"
        # If you don't transform to character here, you will get a list error 
        # When you try to write to file or do anovas
        data.group <- mutate(data.group, group = as.character(file.groupName))
        return(data.group)
}

# Function to calculate descriptive statistics
getStats <- function (data.group){
        data.n    <- length(data.group$cell)
        data.mean <- c()
        data.sd   <- c()
        data.se   <- c()
        data.row  <- colnames(data.group)
        
        # While loop to cycle through all columns containing 
        # numeric data
        i = 2
        maxCol <- length(colnames(data.group)) - 1
        while (i < maxCol){
                data.mean[i-1] <- mean(as.numeric(as.character
                                  (data.group[,i])))
                data.sd[i-1]   <- sd(as.numeric(as.character
                                  (data.group[,i])))
                data.se[i-1]   <- data.sd[i-1]/sqrt(data.n)
                i <- i + 1
        }
        
       return (as.data.frame(cbind(parameter = data.row[2:13], mean = data.mean, 
                                   stdDev = data.sd, stdErr = data.se, data.n)))
}

# Calculate one-way ANOVA results
getANOVA <- function (data.compare){
        
        # Define the Model
        lme.dendriteLength      = lme(dendriteLength ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.dendriteSurfaceArea = lme(dendriteSurfaceArea ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.dendriteVolume      = lme(dendriteVolume ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.dendriteBranches    = lme(dendriteBranches ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.dendriteOrders      = lme(dendriteOrders ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.spineArea           = lme(spineArea ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.spineVolume         = lme(spineVolume ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.totalSpines         = lme(totalSpines ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.thinSpines          = lme(thinSpines ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.stubbySpines        = lme(stubbySpines ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.mushroomSpines      = lme(mushroomSpines ~ group, 
                                  random = ~1|cell, data = data.compare)
        lme.filopodia           = lme(filopodia ~ group, 
                                  random = ~1|cell, data = data.compare)
        
        # Run the ANOVA
        anova.dendriteLength <- anova(lme.dendriteLength)
        anova.dendriteSurfaceArea <- anova(lme.dendriteSurfaceArea)
        anova.dendriteVolume <- anova(lme.dendriteVolume)
        anova.dendriteBranches <- anova(lme.dendriteBranches)
        anova.dendriteOrders <- anova(lme.dendriteOrders)
        anova.spineArea <- anova(lme.spineArea)
        anova.spineVolume <- anova(lme.spineVolume)
        anova.totalSpines <- anova(lme.totalSpines)
        anova.thinSpines <- anova(lme.thinSpines)
        anova.stubbySpines <- anova(lme.stubbySpines)
        anova.mushroomSpines <- anova(lme.mushroomSpines)
        anova.filopodia <- anova(lme.filopodia)
        return(as.data.frame(rbind(dendriteLength = anova.dendriteLength, 
                    dendriteSurfaceArea = anova.dendriteSurfaceArea, 
                    dendriteVolume = anova.dendriteVolume, 
                    dendriteBranches = anova.dendriteBranches, 
                    dendriteOrders = anova.dendriteOrders, 
                    spineArea = anova.spineArea, 
                    spineVolume = anova.spineVolume,
                    totalSpines = anova.totalSpines, 
                    thinSpines = anova.thinSpines, 
                    stubbySpines = anova.stubbySpines,
                    mushroomSpines = anova.mushroomSpines, 
                    filopodia = anova.filopodia)))
}

# Function to determine 95% confidence intervals
getConfidence <- function (data.compare){
        data.row  <- colnames(data.compare)
        conf.dendriteLength       <-groupwiseMean(dendriteLength ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.dendriteSurfaceArea  <-groupwiseMean(dendriteSurfaceArea ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.dendriteVolume       <-groupwiseMean(dendriteVolume ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.dendriteBranches     <-groupwiseMean(dendriteBranches ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.dendriteOrders        <-groupwiseMean(dendriteOrders ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.spineArea            <-groupwiseMean(spineArea ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.spineVolume          <-groupwiseMean(spineVolume ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.totalSpines          <-groupwiseMean(totalSpines ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.thinSpines           <-groupwiseMean(thinSpines ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.stubbySpines         <-groupwiseMean(stubbySpines ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.mushroomSpines       <-groupwiseMean(mushroomSpines ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        conf.filopodia            <-groupwiseMean(filopodia ~ group, 
                                    data=data.compare, conf=0.95, digits=4)
        return(as.data.frame(parameter = data.row, 
                                   rbind(dendriteLength = conf.dendriteLength, 
                                   dendriteSurfaceArea = conf.dendriteSurfaceArea, 
                                   dendriteVolume = conf.dendriteVolume, 
                                   dendriteBranches = conf.dendriteBranches, 
                                   dendriteOrders = conf.dendriteOrders, 
                                   spineArea = conf.spineArea, 
                                   spineVolume = conf.spineVolume,
                                   totalSpines = conf.totalSpines, 
                                   thinSpines = conf.thinSpines, 
                                   stubbySpines = conf.stubbySpines,
                                   mushroomSpines = conf.mushroomSpines, 
                                   filopodia = conf.filopodia)))
}

# Function to write descriptive statistics results to file
writeDescriptive <- function (file.data, file.name, file.path){
        file.name <- paste(file.name, "_descriptive.csv", sep = "")
        file.path <- paste(file.path, "analyzed", sep = "\\")
        if (!dir.exists(file.path)) {dir.create(file.path)}
        file.path <- paste(file.path, file.name, sep = "\\")
        write.csv(file.data, file.path, row.names = FALSE) 
}

# Function to write anova results to file
writeANOVA <- function (file.data, file.name, file.path){
        file.name <- paste(file.name, "_anova.csv", sep = "")
        file.path <- paste(file.path, "analyzed", sep = "\\")
        if (!dir.exists(file.path)) {dir.create(file.path)}
        file.path <- paste(file.path, file.name, sep = "\\")
        write.csv(file.data, file.path, row.names = TRUE) 
}

# Function to write 95% confidence intervals to file
writeConfidence <- function (file.data, file.name, file.path){
        file.name <- paste(file.name, "_confidence.csv", sep = "")
        file.path <- paste(file.path, "analyzed", sep = "\\")
        if (!dir.exists(file.path)) {dir.create(file.path)}
        file.path <- paste(file.path, file.name, sep = "\\")
        write.csv(file.data, file.path, row.names = TRUE) 
}

## Run the script
# Get user input from the console
files <- getUserInput()
file.path <- files[5]

# Get group names
# Assumes each group is named with "..._summary.csv"
group.name1 <- strsplit(files[1], "_summary.csv")
group.name2 <- strsplit(files[2], "_summary.csv")
group.name3 <- strsplit(files[3], "_summary.csv")
group.name4 <- strsplit(files[4], "_summary.csv")

# Read in data from each file and format it for analysis
data.group1  <- getDataTables(files[1],files[5])
data.group2  <- getDataTables(files[2],files[5])
data.group3  <- getDataTables(files[3],files[5])
data.group4  <- getDataTables(files[4],files[5])

# Get descriptive statistics for each group
stats.group1 <- getStats(data.group1)
stats.group2 <- getStats(data.group2)
stats.group3 <- getStats(data.group3)
stats.group4 <- getStats(data.group4)

# Write decriptive statistics to file
write.desc1 <- writeDescriptive(stats.group1, group.name1, file.path) 
write.desc2 <- writeDescriptive(stats.group2, group.name2, file.path) 
write.desc3 <- writeDescriptive(stats.group3, group.name3, file.path) 
write.desc4 <- writeDescriptive(stats.group4, group.name4, file.path) 

# Combine data for comparison between exposure/treatment groups
data.compare1 <- rbind(data.group1, data.group2)
data.compare2 <- rbind(data.group1, data.group3)
data.compare3 <- rbind(data.group1, data.group4)
data.compare4 <- rbind(data.group2, data.group3)
data.compare5 <- rbind(data.group2, data.group4)
data.compare6 <- rbind(data.group3, data.group4)

# Get one-way anova reports for each possible comparison
anova.compare1 <- getANOVA(data.compare1)
anova.compare2 <- getANOVA(data.compare2)
anova.compare3 <- getANOVA(data.compare3)
anova.compare4 <- getANOVA(data.compare4)
anova.compare5 <- getANOVA(data.compare5)
anova.compare6 <- getANOVA(data.compare6)

# Assemble ANOVA output file
file.compare1 <- paste(group.name1,group.name2, sep = "_vs_")
file.compare2 <- paste(group.name1,group.name3, sep = "_vs_")
file.compare3 <- paste(group.name1,group.name4, sep = "_vs_")
file.compare4 <- paste(group.name2,group.name3, sep = "_vs_")
file.compare5 <- paste(group.name2,group.name4, sep = "_vs_")
file.compare6 <- paste(group.name3,group.name4, sep = "_vs_")

# Write anova data to file
writeANOVA (anova.compare1, file.compare1, file.path)
writeANOVA (anova.compare2, file.compare2, file.path)
writeANOVA (anova.compare3, file.compare3, file.path)
writeANOVA (anova.compare4, file.compare4, file.path)
writeANOVA (anova.compare5, file.compare5, file.path)
writeANOVA (anova.compare6, file.compare6, file.path)

# Get confidence intervals for each comparison
conf.group1 <- getConfidence(data.compare1)
conf.group2 <- getConfidence(data.compare2)
conf.group3 <- getConfidence(data.compare3)
conf.group4 <- getConfidence(data.compare4)
conf.group5 <- getConfidence(data.compare5)
conf.group6 <- getConfidence(data.compare6)

# Write confidence intervals to file
writeConfidence (conf.group1, file.compare1, file.path)
writeConfidence (conf.group2, file.compare2, file.path)
writeConfidence (conf.group3, file.compare3, file.path)
writeConfidence (conf.group5, file.compare4, file.path)
writeConfidence (conf.group1, file.compare5, file.path)
writeConfidence (conf.group1, file.compare6, file.path)

