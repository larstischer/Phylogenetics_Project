
# set directory
setwd("C:/Users/lars-/Desktop/Uni/Master/2. Semester/Phylogenetics/Phylogenetics project")

temp <- read.csv("data/T.song.csv")
cd <- read.csv("data/scion_output.csv")
cd <- cd[order(cd$age),]
sr <- read.csv("data/Sr87.macarthur.csv", sep = ";")

## Temperature
# fixing the table because there was a wrong entry
temp[19626,]$Stage = "Induan"  

# Subset to Mesozoic and Cenozoic
postPal_T <- subset(temp, temp$Age <= 251.89)

stAge <- postPal_T[, c("Stage","Age")]


# mean temperature per stage 
stages <- unique(postPal_T$Stage)
meanT <- NULL

for(i in 1:length(stages)){
  currentStage <- stages[i]
  stageSub <- postPal_T[postPal_T$Stage == currentStage, ]
  stageT <- mean(stageSub$Temp, na.rm = TRUE)
  meanT <- c(meanT, stageT)
}
names(meanT) <- stages
meanT
length(meanT)

## CO2
# Subset
library(dplyr)
postPal_CO2 <- subset(cd, cd$age <= 251.89)

ten <- subset(postPal_CO2, postPal_CO2$age > 0 & postPal_CO2$age <= 10)

a <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250)
i <- 1
p <- colnames(postPal_CO2)
sub <- data.frame(Age = a[2:26], CO2 = NA )
#p <- colnames(postPal_CO2)
#sub <- as.data.frame(sub)

for(i in 1:(length(a)-1)){
  ten <- subset(postPal_CO2,postPal_CO2$age > a[i] & postPal_CO2$age <= a[i+1] )
  int <- mean(ten$CO2.PPM., na.rm = T)
  sub[i,2] <- int
}

## Sr
# Subset to Mesozoic and Cenozoic
postPal_Sr <- subset(sr, sr$age <= 251.89)

ten_Sr <- subset(postPal_Sr, postPal_Sr$age > 0 & postPal_Sr$age <= 10)

a <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250)
i <- 1
s <- colnames(postPal_Sr)
sub2 <- data.frame(Age = a[2:26], Sr = NA )
#p <- colnames(postPal_CO2)
#sub <- as.data.frame(sub)

for(i in 1:(length(a)-1)){
  ten_Sr <- subset(postPal_Sr,postPal_Sr$age > a[i] & postPal_Sr$age <= a[i+1] )
  mint <- mean(ten_Sr$Mean, na.rm = T)
  sub2[i,2] <- mint
}

## Overview tables
lookup <- sub[, c("Age", "CO2")]
lookup$Sr <- sub2[, c("Sr")]

write.csv2(lookup, file = "data/lookup.csv", row.names = FALSE)

lookup2 <- data.frame(stages)
lookup2$Temperature <- meanT

write.csv2(lookup2, file = "data/lookup2.csv", row.names = FALSE)
