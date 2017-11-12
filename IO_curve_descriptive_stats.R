# Input/Output Curve Descriptive Statistics
# Author: Haley Speed, PhD
# University of Texas Southwestern Medical Center, Dallas

#-------------------------------------------------------------------------------------------------------#
# Example Spreadsheet Setup                                                                             #
#-------------------------------------------------------------------------------------------------------#
#                                                                                                       #
#        |   A  |   B   |     C     |     D    |     E   |                                              #
#       __________________________________________________                                              #
#       1|  id  | group | intensity |   slope  |    fv   |                                              #
#       2| A1.1 |  WT   |     0     | -0.00323 | -0.0043 |                                              #
#       3| A1.1 |  WT   |     50    | -0.03421 | -0.0211 |                                              #
#       4| A1.1 |  WT   |     100   | -0.52671 | -0.4521 |                                              #
#       2| A2.2 |  KI   |     0     | -0.00256 | -0.0023 |                                              #
#       3| A2.2 |  KI   |     50    | -0.07532 | -0.0215 |                                              #
#       4| A2.2 |  KI   |     100   | -0.64367 | -0.5212 |                                              #
#                                                                                                       #
#-------------------------------------------------------------------------------------------------------#
                        


#-------------------------------------------------------------------------------------------------------#
# Things the user should change                                                                         #
#-------------------------------------------------------------------------------------------------------#

# Assign filename and path (replace all backslashes with double forward slashes in path)
fileName  <- "C://Users"

# Import data from csv spreadsheet
data <- read.csv(filename, header = TRUE)

# Import data from xlsx spreadsheet
#library(openxlsx)
#data <- read.xlsx(filename, sheet = 1, colNames = TRUE) 

# Identify columns in spreadsheet
col.id        <- 1    # Slice Number
col.group     <- 2    # Group Name
col.intensity <- 3    # Intensity
col.slope     <- 4    # Slope
col.fv        <- 5    # Fiber Volley



#-------------------------------------------------------------------------------------------------------#
# Do NOT MODIFY                                                                     			#
#-------------------------------------------------------------------------------------------------------#
library(dplyr)
library(rcompanion)

# Get names of groups
groupNames <- unique(data[,col.group])
intensity <- unique(data[col.intensity])

# Create empty output table
output <- data.frame(group = character(), n = numeric(), intensity = factor(), slope.mean = numeric(), slope.sd = numeric(), 
                     slope.se = numeric(), fv.mean = numeric(), fv.sd = numeric(), fv.se = numeric())


# For each group, determine basic stats
for (x in groupNames) {
        # Filter rows containing data for the current group
        temp.data <- filter(data, group == x)

        # Get the sample size
        n <- tapply(as.numeric(temp.data[,col.group]),temp.data[,col.intensity], length)
        
        # 
        slope.mean  <- tapply(as.numeric(temp.data[,col.slope]),temp.data[,col.intensity], mean)
        slope.sd  <- tapply(as.numeric(temp.data[,col.slope]),temp.data[,col.intensity], sd)
        slope.se <- slope.sd/sqrt(n)
        
        fv.mean  <- tapply(as.numeric(temp.data[,col.fv]),temp.data[,col.intensity], mean)
        fv.sd  <- tapply(as.numeric(temp.data[,col.fv]),temp.data[,col.intensity], sd)
        fv.se <- fv.sd/sqrt(n)
        
        slope.conf <-groupwiseMean(slope ~ intensity, data=temp.data, conf=0.95, digits=4) 
        fv.conf <-groupwiseMean(fv ~ intensity, data=temp.data, conf=0.95, digits=4) 
        
        row <- cbind(group = x, n = n, intensity = intensity, slope.mean = slope.mean, 
                     slope.sd = slope.sd, slope.se = slope.se, slope.conf5 = slope.conf$Trad.lower, 
                     slope.conf95 = slope.conf$Trad.upper, fv.mean = fv.mean, fv.sd = fv.sd, 
                     fv.se = fv.se, fv.conf5 = fv.conf$Trad.lower, fv.conf95 = fv.conf$Trad.upper)
        output <- rbind(output, row)
}

# Create output folder
fileDir <- dirname(filename)
outDir <- paste(fileDir, "analyzed", sep = "//")

# Create output file name
fileBase <- basename(fileName)
fileBase <- strsplit(fileBase, ".csv")
#filebase <- strsplit(filebase, ".xlsx")
fileBase <- paste(fileBase, "descriptive.csv", sep = "_")
outFile <- paste(outDir, fileBase, sep = "//")

# Create Analyzed directory
if (!dir.exists(outDir)) { dir.create(outDir)}


# Write analyzed data to file
write.csv (output, outFile,  row.names = FALSE) 
