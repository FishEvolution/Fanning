#### Setup ####

library(rfishbase)
library(morphospace, lib="C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R")
library(dplyr)
library(Momocs)
library(geomorph)
library(ape)
library(phangorn)
library(OUwie, lib="C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R")
library(picante)
library(regclass)
library(MCMCglmm)
library(terra)
library(raster)
library(sf)
library(Morpho)
library(motmot)
library(caper)

#### Swimming Method Data Collection ####

Swimming <- read.csv("Species.csv")

SwimmingMethod <- swimming(species_list=Swimming$FishBase.Species)

write.csv(SwimmingMethod, 'SwimmingMethod.csv')

#### Environment Data Collection ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Maps")

library(terra)

Groups <- read.csv("Groups.csv", header = F)

Species <- read.csv("Species.csv", header = F)

# Loop 1

N <- dim(Groups)[1]

for(i in 1:N){
  
  Name <- Groups$V1[i]
  Name.zip <- Groups$V2[i]
  
  unzip(Name.zip, exdir = Name)
  
  print(i)
}

# Loop 2

N <- dim(Species)[1]

Fish <- vect(Groups$V3[11])
  
for(i in 1:N){
    
  FishSpecies <- Fish[Fish$sci_name == Species$V1[i],]
    
  if (nrow(FishSpecies) > 0) {
    writeVector(FishSpecies, Species$V2[i])
  }
    
  print(i)
    
}

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Maps/Species Files")

x <- list.files()
x
Temperature <- x[765]
Oxygen <- x[553]

Temp <- raster(Temperature)
Oxy <- raster(Oxygen)
Chlor <- raster(Chlorophyll)

ListMaps<-read.csv("RangeMaps.csv", header=FALSE)
ListMaps<-as.vector(ListMaps)

N <- length(ListMaps)
TempAnswersMean <- rep(0,N)
OxyAnswersMean <- rep(0,N)
TempAnswersSD <- rep(0,N)
OxyAnswersSD <- rep(0,N)

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Temp,Shape),Shape)
  answer<-mean(getValues(y2),na.rm=TRUE) # Take SD, Simulate from Normal Distribution. MCMCglmm - 10 values per species
  
  TempAnswersMean[i]<-mean(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(TempAnswersMean, file=paste("ZZTempAnswersMean.csv"))

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Oxy,Shape),Shape)
  answer<-mean(getValues(y2),na.rm=TRUE)
  
  OxyAnswersMean[i]<-mean(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(OxyAnswersMean, file=paste("ZZOxyAnswersMean.csv"))

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Temp,Shape),Shape)
  answer<-SD(getValues(y2),na.rm=TRUE)
  
  TempAnswersSD[i]<-SD(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(ChlorAnswersSD, file=paste("ZZTempAnswersSD.csv"))

for(i in 1:N){
  Shape <- st_read(ListMaps[i])
  Shape<-Shape[1,1]
  y2<-mask(crop(Oxy,Shape),Shape)
  answer<-SD(getValues(y2),na.rm=TRUE)
  
  OxyAnswersSD[i]<-SD(getValues(y2),na.rm=TRUE)
  print(i)
}

write.csv(OxyAnswersSD, file=paste("ZZOxyAnswersSD.csv"))

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests")

EnvironmentalData <- read.csv("EnvironmentalData.csv")

N = 167

EnvironmentalDataRanges <- matrix(0, nrow = 20, ncol = 167)

for(i in 1:N){
  
  NumbersT <- rnorm(n = 10, mean = EnvironmentalData$TemperatureMean[i], sd = EnvironmentalData$TemperatureSD[i])
  NumbersO <- rnorm(n = 10, mean = EnvironmentalData$OxygenMean[i], sd = EnvironmentalData$OxygenSD[i])
  
  EnvironmentalDataRanges[1:10,i] <- NumbersT
  EnvironmentalDataRanges[11:20,i] <- NumbersO
  
  print(i)
  
}

write.csv(EnvironmentalDataRanges, "EnvironmentalDataRanges.csv")

#### Consensus Tree ####

# Sample Size 162

Phylogeny <- read.tree("Actinopterygii.trees")

GroupNames <- read.csv("GroupNames.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Phylogeny[[1]]$tip.label)

for(k in 1:100)
{
  Phylogeny[[k]] <- keep.tip(Phylogeny[[k]], Phylogeny[[k]]$tip.label[match(GroupNames$V1, Phylogeny[[k]]$tip.label)])
  print(k)
}

write.nexus(Phylogeny, file="Phylogeny.nex")

# Sample Size 244

Phylogeny <- read.tree("Actinopterygii.trees")

GroupNames <- read.csv("GroupNames1.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Phylogeny[[1]]$tip.label)

for(k in 1:100)
{
  Phylogeny[[k]] <- keep.tip(Phylogeny[[k]], Phylogeny[[k]]$tip.label[match(GroupNames$V1, Phylogeny[[k]]$tip.label)])
  print(k)
}

write.nexus(Phylogeny, file="Phylogeny1.nex")

#### Morphospace ####

source("Morphospace.R")

## Fanning - Sample Size = 162

# Import the harmonic array and data to generate a complete morphospace.

HarmonicArray <- read.csv("HarmonicArray.csv", header = FALSE)
coefmat <- import.will(HarmonicArray)
FishData <- read.csv("FishData.csv")
rownames(coefmat) <- FishData$Name
index <- rownames(coefmat) %in% FishData$Name 
FishShapes <- coefmat[index,] 

msp <- mspace(FishShapes, mag = 1, nh = 12, nv = 10, bg.models = "gray",
              invax=c(1), size.models = 3) %>%
  proj_shapes(FishShapes, pch = 16)

write.csv(msp$ordination$x, 'PCValues.csv')

# Graph of PC1 Abundance

hist(FishData$PC1, breaks = 20)

# Calculate area of morphospace occupation for fanning and non-fanning fish.

str(FishData)
FishData$Group = factor(FishData$Group)
fac = factor(FishData$Group)
xy = cbind(FishData$PC1,FishData$PC2)

# Non-fanning

taxclass.indexNF = FishData$Group == "NonFanning"
NFSubData = FishData[taxclass.indexNF,]

col = "salmon3"
alpha = .5
lty = 1
NFx <- xy[fac == levels(fac)[2], 1]
NFy <- xy[fac == levels(fac)[2], 2]
NFHullp <- grDevices::chull(x = NFx, y = NFy)
graphics::polygon(NFx[NFHullp], NFy[NFHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

NFHull.coo <- cbind(NFx[NFHullp], NFy[NFHullp]) 
Obs.NFArea = Momocs::coo_area(coo = NFHull.coo)
Obs.NFArea 

# Fanning

taxclass.indexF = FishData$Group
FSubData = FishData[taxclass.indexF,]

col = "skyblue2"
alpha = .5
lty = 1
Fx <- xy[fac == levels(fac)[1], 1]
Fy <- xy[fac == levels(fac)[1], 2]
FHullp <- grDevices::chull(x = Fx, y = Fy)
graphics::polygon(Fx[FHullp], Fy[FHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

FHull.coo <- cbind(Fx[FHullp], Fy[FHullp]) 
Obs.FArea = Momocs::coo_area(coo = FHull.coo)
Obs.FArea 

# Area - Fanning and Non-Fanning

taxclass.index = FishData$All == "All"
SubData = FishData[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 81, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1,Random.SubData$PC2)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.FArea, col = 2, lwd = 2, lty = 2)
abline(v = Obs.NFArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.FArea) / (perm + 1)

sum(Random.Areas < Obs.NFArea) / (perm + 1)

# Area - Overlap

taxclass.Oindex = FishData$All == "All"
OSubData = FishData[taxclass.Oindex,]

perm = 9999
Random.OAreas = NULL
for(i in 1:perm) {
  
  Random.OIndex = sample(x = 1:nrow(OSubData), size = 109, replace = FALSE)
  
  Random.OSubData = OSubData[Random.OIndex,]
  
  Random.Oxy = cbind(Random.OSubData$PC1,Random.OSubData$PC2)
  Random.Ox <- Random.Oxy[, 1]
  Random.Oy <- Random.Oxy[, 2]
  Random.OHullp <- grDevices::chull(x = Random.Ox, y = Random.Oy)
  graphics::polygon(Random.Ox[Random.OHullp], Random.Oy[Random.OHullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.OArea = coo_area(coo = Random.Oxy[Random.OHullp,])
  Random.OAreas[i] = Random.OArea
}

hist(Random.OAreas, breaks = 100)
abline(v = 0.10367743, col = 2, lwd = 2, lty = 2)

sum(Random.OAreas < 0.10367743) / (perm + 1)

# PC1 Shape Transformation

Extshps1 <- extract_shapes(mspace = msp, axis = 1, nshapes = 2)
Extshps_Neg1 <- inv_efourier(Extshps1$shapes[1,])
Extshps_Pos1 <- inv_efourier(Extshps1$shapes[2,])
plot(Extshps_Pos1, asp = 5, type = "l", col = "red")
points(Extshps_Neg1, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg1)



# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg2)

## Pectoral Swimming  - Sample Size = 162

# Import the harmonic array and data to generate a complete morphospace.

HarmonicArrayS <- read.csv("HarmonicArrayS.csv", header = FALSE)
coefmat <- import.will(HarmonicArrayS)
FishData <- read.csv("FishData.csv")
FishData <- FishData[-(which(FishData$PectoralSwimming %in% "-")),]
rownames(coefmat) <- FishData$Name
index <- rownames(coefmat) %in% FishData$Name 
FishShapes <- coefmat[index,] 

msp <- mspace(FishShapes, mag = 1, nh = 12, nv = 10, bg.models = "gray",
              invax=c(1), size.models = 3) %>%
  proj_shapes(FishShapes, pch = 16)

write.csv(msp$ordination$x, 'PCValues.csv')

# Calculate area of morphospace occupation for pectoral swimming and non-pectoral fish.

str(FishData)
FishData$PectoralSwimming = factor(FishData$PectoralSwimming)
fac = factor(FishData$PectoralSwimming)
xy = cbind(FishData$PC1S,FishData$PC2S)

# Non-Pectoral Swimming

taxclass.indexNS = FishData$PectoralSwimming
FSubData = FishData[taxclass.indexNS,]

col = "orchid2"
alpha = .5
lty = 1
NSx <- xy[fac == levels(fac)[1], 1]
NSy <- xy[fac == levels(fac)[1], 2]
NSHullp <- grDevices::chull(x = NSx, y = NSy)
graphics::polygon(NSx[NSHullp], NSy[NSHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

NSHull.coo <- cbind(NSx[NSHullp], NSy[NSHullp]) 
Obs.NSArea = Momocs::coo_area(coo = NSHull.coo)
Obs.NSArea 

# Pectoral Swimming

taxclass.indexS = FishData$PectoralSwimming == "0"
FSubData = FishData[taxclass.indexS,]

col = "palegreen2"
alpha = .5
lty = 1
Sx <- xy[fac == levels(fac)[2], 1]
Sy <- xy[fac == levels(fac)[2], 2]
SHullp <- grDevices::chull(x = Sx, y = Sy)
graphics::polygon(Sx[SHullp], Sy[SHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

SHull.coo <- cbind(Sx[SHullp], Sy[SHullp]) 
Obs.SArea = Momocs::coo_area(coo = SHull.coo)
Obs.SArea 

# Area - Pectoral and Non-Pectoral

taxclass.index = FishData$All == "All"
SubData = FishData[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 32, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1S,Random.SubData$PC2S)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.SArea, col = 2, lwd = 2, lty = 2)
abline(v = Obs.NSArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.SArea) / (perm + 1)

sum(Random.Areas < Obs.NSArea) / (perm + 1)

# Area - Overlap

taxclass.Oindex = FishData$All == "All"
OSubData = FishData[taxclass.Oindex,]

perm = 9999
Random.OAreas = NULL
for(i in 1:perm) {
  
  Random.OIndex = sample(x = 1:nrow(OSubData), size = 22, replace = FALSE)
  
  Random.OSubData = OSubData[Random.OIndex,]
  
  Random.Oxy = cbind(Random.OSubData$PC1S,Random.OSubData$PC2S)
  Random.Ox <- Random.Oxy[, 1]
  Random.Oy <- Random.Oxy[, 2]
  Random.OHullp <- grDevices::chull(x = Random.Ox, y = Random.Oy)
  graphics::polygon(Random.Ox[Random.OHullp], Random.Oy[Random.OHullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.OArea = coo_area(coo = Random.Oxy[Random.OHullp,])
  Random.OAreas[i] = Random.OArea
}

hist(Random.OAreas, breaks = 100)
abline(v = 0.031143575, col = 2, lwd = 2, lty = 2)

sum(Random.OAreas < 0.031143575) / (perm + 1)

# PC1 Shape Transformation

Extshps1 <- extract_shapes(mspace = msp, axis = 1, nshapes = 2)
Extshps_Neg1 <- inv_efourier(Extshps1$shapes[1,])
Extshps_Pos1 <- inv_efourier(Extshps1$shapes[2,])
plot(Extshps_Pos1, asp = 5, type = "l", col = "red")
points(Extshps_Neg1, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg2)

## Fanning 1  - Sample Size = 244

# Import the harmonic array and data to generate a complete morphospace.

HarmonicArray1 <- read.csv("HarmonicArray1.csv", header = FALSE)
coefmat <- import.will(HarmonicArray1)
FishData1 <- read.csv("FishData1.csv")
rownames(coefmat) <- FishData1$Name
index <- rownames(coefmat) %in% FishData1$Name 
FishShapes <- coefmat[index,] 

msp <- mspace(FishShapes, mag = 1, nh = 12, nv = 10, bg.models = "gray",
              invax=c(1), size.models = 3) %>%
  proj_shapes(FishShapes, pch = 16)

write.csv(msp$ordination$x, 'PCValues.csv')

# Graph of PC1 Abundance

hist(FishData1$PC1, breaks = 161)

# Calculate area of morphospace occupation for fanning and non-fanning fish.

str(FishData1)
FishData1$Group = factor(FishData1$Group)
fac = factor(FishData1$Group)
xy = cbind(FishData1$PC1,FishData1$PC2)

# Non-fanning 1

taxclass.indexNF = FishData1$Group == "NonFanning"
NFSubData = FishData1[taxclass.indexNF,]

col = "salmon3"
alpha = .5
lty = 1
NFx <- xy[fac == levels(fac)[2], 1]
NFy <- xy[fac == levels(fac)[2], 2]
NFHullp <- grDevices::chull(x = NFx, y = NFy)
graphics::polygon(NFx[NFHullp], NFy[NFHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

NFHull.coo <- cbind(NFx[NFHullp], NFy[NFHullp]) 
Obs.NFArea = Momocs::coo_area(coo = NFHull.coo)
Obs.NFArea 

# Fanning 1

taxclass.indexF = FishData1$Group
FSubData = FishData1[taxclass.indexF,]

col = "skyblue2"
alpha = .5
lty = 1
Fx <- xy[fac == levels(fac)[1], 1]
Fy <- xy[fac == levels(fac)[1], 2]
FHullp <- grDevices::chull(x = Fx, y = Fy)
graphics::polygon(Fx[FHullp], Fy[FHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

FHull.coo <- cbind(Fx[FHullp], Fy[FHullp]) 
Obs.FArea = Momocs::coo_area(coo = FHull.coo)
Obs.FArea 

# Area - Non-Fanning

taxclass.index = FishData1$All == "All"
SubData = FishData1[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 163, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1,Random.SubData$PC2)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.NFArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.NFArea) / (perm + 1)

# Area - Fanning

taxclass.index = FishData1$All == "All"
SubData = FishData1[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 81, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1,Random.SubData$PC2)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.FArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.FArea) / (perm + 1)

# Area - Overlap

taxclass.Oindex = FishData1$All == "All"
OSubData = FishData1[taxclass.Oindex,]

perm = 9999
Random.OAreas = NULL
for(i in 1:perm) {
  
  Random.OIndex = sample(x = 1:nrow(OSubData), size = 149, replace = FALSE)
  
  Random.OSubData = OSubData[Random.OIndex,]
  
  Random.Oxy = cbind(Random.OSubData$PC1,Random.OSubData$PC2)
  Random.Ox <- Random.Oxy[, 1]
  Random.Oy <- Random.Oxy[, 2]
  Random.OHullp <- grDevices::chull(x = Random.Ox, y = Random.Oy)
  graphics::polygon(Random.Ox[Random.OHullp], Random.Oy[Random.OHullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.OArea = coo_area(coo = Random.Oxy[Random.OHullp,])
  Random.OAreas[i] = Random.OArea
}

hist(Random.OAreas, breaks = 100)
abline(v = 0.128251902, col = 2, lwd = 2, lty = 2)

sum(Random.OAreas < 0.128251902) / (perm + 1)

# PC1 Shape Transformation

Extshps1 <- extract_shapes(mspace = msp, axis = 1, nshapes = 2)
Extshps_Neg1 <- inv_efourier(Extshps1$shapes[1,])
Extshps_Pos1 <- inv_efourier(Extshps1$shapes[2,])
plot(Extshps_Pos1, asp = 5, type = "l", col = "red")
points(Extshps_Neg1, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg2)

# Pectoral Swimming 1   - Sample Size = 244

# Import the harmonic array and data to generate a complete morphospace.

HarmonicArrayS1 <- read.csv("HarmonicArrayS1.csv", header = FALSE)
coefmat <- import.will(HarmonicArrayS1)
FishData1 <- read.csv("FishData1.csv")
FishData1 <- FishData1[-(which(FishData1$PectoralSwimming %in% "-")),]
rownames(coefmat) <- FishData1$Name
index <- rownames(coefmat) %in% FishData1$Name 
FishShapes <- coefmat[index,] 

msp <- mspace(FishShapes, mag = 1, nh = 12, nv = 10, bg.models = "gray",
              , size.models = 3) %>%
  proj_shapes(FishShapes, pch = 16)

write.csv(msp$ordination$x, 'PCValues.csv')

# Calculate area of morphospace occupation for pectoral swimming and non-pectoral fish.

str(FishData1)
FishData1$PectoralSwimming = factor(FishData1$PectoralSwimming)
fac = factor(FishData1$PectoralSwimming)
xy = cbind(FishData1$PC1S,FishData1$PC2S)

# Non-Pectoral Swimming 1

taxclass.indexNS = FishData1$PectoralSwimming
FSubData = FishData1[taxclass.indexNS,]

col = "orchid2"
alpha = .5
lty = 1
NSx <- xy[fac == levels(fac)[1], 1]
NSy <- xy[fac == levels(fac)[1], 2]
NSHullp <- grDevices::chull(x = NSx, y = NSy)
graphics::polygon(NSx[NSHullp], NSy[NSHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

NSHull.coo <- cbind(NSx[NSHullp], NSy[NSHullp]) 
Obs.NSArea = Momocs::coo_area(coo = NSHull.coo)
Obs.NSArea 

# Pectoral Swimming

taxclass.indexS = FishData$PectoralSwimming == "0"
FSubData = FishData[taxclass.indexS,]

col = "palegreen2"
alpha = .5
lty = 1
Sx <- xy[fac == levels(fac)[2], 1]
Sy <- xy[fac == levels(fac)[2], 2]
SHullp <- grDevices::chull(x = Sx, y = Sy)
graphics::polygon(Sx[SHullp], Sy[SHullp], border = col, 
                  col = grDevices::adjustcolor(col, alpha.f = alpha), lty = lty)

SHull.coo <- cbind(Sx[SHullp], Sy[SHullp]) 
Obs.SArea = Momocs::coo_area(coo = SHull.coo)
Obs.SArea 

# Area - Non-Pectoral Swimming 1

taxclass.index = FishData1$All == "All"
SubData = FishData1[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 112, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1S,Random.SubData$PC2S)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.NSArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.NSArea) / (perm + 1)

# Area - Pectoral Swimming 1

taxclass.index = FishData1$All == "All"
SubData = FishData1[taxclass.index,]

perm = 9999
Random.Areas = NULL
for(i in 1:perm) {
  
  Random.Index = sample(x = 1:nrow(SubData), size = 34, replace = FALSE)
  
  Random.SubData = SubData[Random.Index,]
  
  Random.xy = cbind(Random.SubData$PC1S,Random.SubData$PC2S)
  Random.x <- Random.xy[, 1]
  Random.y <- Random.xy[, 2]
  Random.Hullp <- grDevices::chull(x = Random.x, y = Random.y)
  graphics::polygon(Random.x[Random.Hullp], Random.y[Random.Hullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.Area = coo_area(coo = Random.xy[Random.Hullp,])
  Random.Areas[i] = Random.Area
}

hist(Random.Areas, breaks = 100)
abline(v = Obs.SArea, col = 2, lwd = 2, lty = 2)

sum(Random.Areas < Obs.SArea) / (perm + 1)

# Area - Overlap

taxclass.Oindex = FishData1$All == "All"
OSubData = FishData1[taxclass.Oindex,]

perm = 9999
Random.OAreas = NULL
for(i in 1:perm) {
  
  Random.OIndex = sample(x = 1:nrow(OSubData), size = 55, replace = FALSE)
  
  Random.OSubData = OSubData[Random.OIndex,]
  
  Random.Oxy = cbind(Random.OSubData$PC1,Random.OSubData$PC2)
  Random.Ox <- Random.Oxy[, 1]
  Random.Oy <- Random.Oxy[, 2]
  Random.OHullp <- grDevices::chull(x = Random.Ox, y = Random.Oy)
  graphics::polygon(Random.Ox[Random.OHullp], Random.Oy[Random.OHullp], border = col, 
                    col = grDevices::adjustcolor(col, alpha.f = alpha), 
                    lty = lty)
  
  Random.OArea = coo_area(coo = Random.Oxy[Random.OHullp,])
  Random.OAreas[i] = Random.OArea
}

hist(Random.OAreas, breaks = 100)
abline(v = 0.073639351, col = 2, lwd = 2, lty = 2)

sum(Random.OAreas < 0.073639351) / (perm + 1)

# PC1 Shape Transformation

Extshps1 <- extract_shapes(mspace = msp, axis = 1, nshapes = 2)
Extshps_Neg1 <- inv_efourier(Extshps1$shapes[1,])
Extshps_Pos1 <- inv_efourier(Extshps1$shapes[2,])
plot(Extshps_Pos1, asp = 5, type = "l", col = "red")
points(Extshps_Neg1, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, 
                gridPars = gridPar(pt.bg = "green3", pt.size = 1.25),
                method = "vector", mag = 0.5)
points(Extshps_Neg2)

#### OUwie - Adaptive Landscape ####

Phylogeny <- read.tree("Actinopterygii.trees")

## Fanning - Sample Size = 162

# Model Testing

OUwieData <- read.csv("OUwie.csv")
GroupNames <- read.csv("GroupNames.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Fanning<-rep("Fanning", dim(OUwieData)[1])
Fanning[OUwieData$Fanning == "Fanning"] <- "Fanning"
Fanning[OUwieData$Fanning == "NonFanning"] <- "NonFanning"
Fanning <- as.factor(Fanning)
names(Fanning) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Fanning, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Fanning, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Fanning)
  anc <- rep("NonFanning", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Fanning"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("BM1"),algorithm=c("invert"))
  answers[k,2]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("BMS"))
  answers[k,3]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OU1"))
  answers[k,4]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMV"))
  answers[k,6]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMA"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,11]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,12]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMA"))
  answers[k,13]<-test$AICc
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputs.csv"))

# Fanning Peak Values

OUwieData <- read.csv("OUwie.csv")
GroupNames <- read.csv("GroupNames.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Fanning<-rep("Fanning", dim(OUwieData)[1])
Fanning[OUwieData$Fanning == "Fanning"] <- "Fanning"
Fanning[OUwieData$Fanning == "NonFanning"] <- "NonFanning"
Fanning <- as.factor(Fanning)
names(Fanning) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Fanning, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Fanning, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Fanning)
  anc <- rep("NonFanning", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Fanning"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,2]<-test$AICc
  answers[k,3]<-test$theta[1,1]
  answers[k,4]<-test$theta[2,1]
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  answers[k,6]<-test$theta[1,1]
  answers[k,7]<-test$theta[2,1]
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsPC.csv"))

## Swimming - Sample Size = 162

OUwieData <- read.csv("OUwieS.csv")
GroupNames <- read.csv("GroupNamesS.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Swimming<-rep("Swimming", dim(OUwieData)[1])
Swimming[OUwieData$Swimming == "Swimming"] <- "Swimming"
Swimming[OUwieData$Swimming == "NonSwimming"] <- "NonSwimming"
Swimming <- as.factor(Swimming)
names(Swimming) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Swimming, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Swimming, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Swimming)
  anc <- rep("NonSwimming", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Swimming"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("BM1"),algorithm=c("invert"))
  answers[k,2]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("BMS"))
  answers[k,3]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OU1"))
  answers[k,4]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMV"))
  answers[k,6]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMA"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,11]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,12]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMA"))
  answers[k,13]<-test$AICc
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsS.csv"))

# Swimming Peak Values

OUwieData <- read.csv("OUwieS.csv")
GroupNames <- read.csv("GroupNamesS.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Swimming<-rep("Swimming", dim(OUwieData)[1])
Swimming[OUwieData$Swimming == "Swimming"] <- "Swimming"
Swimming[OUwieData$Swimming == "NonSwimming"] <- "NonSwimming"
Swimming <- as.factor(Swimming)
names(Swimming) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Swimming, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Swimming, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Swimming)
  anc <- rep("Swimming", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "NonSwimming"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,2]<-test$AICc
  answers[k,3]<-test$theta[1,1]
  answers[k,4]<-test$theta[2,1]
  
  test<-OUwie(Tree,trait2.data,model=c("OUMA"))
  answers[k,5]<-test$AICc
  answers[k,6]<-test$theta[1,1]
  answers[k,7]<-test$theta[2,1]
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsSPC.csv"))

## Fanning 1 - Sample Size = 244

OUwieData <- read.csv("OUwie1.csv")
GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Fanning<-rep("Fanning", dim(OUwieData)[1])
Fanning[OUwieData$Fanning == "Fanning"] <- "Fanning"
Fanning[OUwieData$Fanning == "NonFanning"] <- "NonFanning"
Fanning <- as.factor(Fanning)
names(Fanning) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Fanning, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Fanning, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Fanning)
  anc <- rep("NonFanning", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Fanning"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("BM1"),algorithm=c("invert"))
  answers[k,2]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("BMS"))
  answers[k,3]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OU1"))
  answers[k,4]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMV"))
  answers[k,6]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMA"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,11]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,12]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMA"))
  answers[k,13]<-test$AICc
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputs1.csv"))

# Fanning 1 Peak Values

OUwieData <- read.csv("OUwie1.csv")
GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Fanning<-rep("Fanning", dim(OUwieData)[1])
Fanning[OUwieData$Fanning == "Fanning"] <- "Fanning"
Fanning[OUwieData$Fanning == "NonFanning"] <- "NonFanning"
Fanning <- as.factor(Fanning)
names(Fanning) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Fanning, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Fanning, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Fanning)
  anc <- rep("NonFanning", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Fanning"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("OUMV"))
  answers[k,2]<-test$AICc
  answers[k,3]<-test$theta[1,1]
  answers[k,4]<-test$theta[2,1]
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  answers[k,6]<-test$theta[1,1]
  answers[k,7]<-test$theta[2,1]
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputs1PC.csv"))

## Swimming 1 - Sample Size = 244

OUwieData <- read.csv("OUwieS1.csv")
GroupNames <- read.csv("GroupNamesS1.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Swimming<-rep("Swimming", dim(OUwieData)[1])
Swimming[OUwieData$Swimming == "Swimming"] <- "Swimming"
Swimming[OUwieData$Swimming == "NonSwimming"] <- "NonSwimming"
Swimming <- as.factor(Swimming)
names(Swimming) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Swimming, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Swimming, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Swimming)
  anc <- rep("NonSwimming", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Swimming"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("BM1"),algorithm=c("invert"))
  answers[k,2]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("BMS"))
  answers[k,3]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OU1"))
  answers[k,4]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMV"))
  answers[k,6]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUMA"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,11]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,12]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMA"))
  answers[k,13]<-test$AICc
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsS1.csv"))

# Swimming 1 Peak Values

OUwieData <- read.csv("OUwieS1.csv")
GroupNames <- read.csv("GroupNamesS1.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1
names(trait1)<-OUwieData$Genus_species

trait2<-OUwieData$PC2
names(trait2)<-OUwieData$Genus_species

Swimming<-rep("Swimming", dim(OUwieData)[1])
Swimming[OUwieData$Swimming == "Swimming"] <- "Swimming"
Swimming[OUwieData$Swimming == "NonSwimming"] <- "NonSwimming"
Swimming <- as.factor(Swimming)
names(Swimming) <- OUwieData$Genus_species

for(k in 1:100){
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Swimming, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Swimming, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Swimming)
  anc <- rep("NonSwimming", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Swimming"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Swimming=Swimming, X = trait2)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,2]<-test$AICc
  answers[k,3]<-test$theta[1,1]
  answers[k,4]<-test$theta[2,1]
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,5]<-test$AICc
  answers[k,6]<-test$theta[1,1]
  answers[k,7]<-test$theta[2,1]
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsS1PC.csv"))
#### Intraspecific Variation ####

source("Morpho.R")

HarmonicArray <- read.csv("HarmonicArrayVariation.csv", header = FALSE)
coefmat <- import.will(HarmonicArray)
FishData <- read.csv("IntraspecificVariation.csv")
rownames(coefmat) <- FishData$Species
index <- rownames(coefmat) %in% FishData$Species 
FishShapes <- coefmat[index,] 

msp <- mspace(FishShapes, mag = 1, nh = 12, nv = 10, bg.models = "gray",
              invax=c(1), size.models = 3) %>%
  proj_shapes(FishShapes, pch = 16)

write.csv(msp$ordination$x, 'PCValues.csv')

# OUwie

Environment <- read.csv("IntraspecificVariation.csv")

A <- subset(Environment, Sample == "DI")
B <- subset(Environment, Sample == "D")
C <- subset(Environment, Sample == "CK")
D <- subset(Environment, Sample == "CH")
E <- subset(Environment, Sample == "HS")
F <- subset(Environment, Sample == "A")
G <- subset(Environment, Sample == "DR")
H <- subset(Environment, Sample == "AM")
I <- subset(Environment, Sample == "FL")
J <- subset(Environment, Sample == "BY")

sdA <- sd(A$PC1)
sdB <- sd(B$PC1)
sdC <- sd(C$PC1)
sdD <- sd(D$PC1)
sdE <- sd(E$PC1)
sdF <- sd(F$PC1)
sdG <- sd(G$PC1)
sdH <- sd(H$PC1)
sdI <- sd(I$PC1)
sdJ <- sd(J$PC1)

(sdA + sdB + sdC + sdD + sdE + sdF + sdG + sdH + sdI + sdJ)/10

# 0.03309501

sdA <- sd(A$PC2)
sdB <- sd(B$PC2)
sdC <- sd(C$PC2)
sdD <- sd(D$PC2)
sdE <- sd(E$PC2)
sdF <- sd(F$PC2)
sdG <- sd(G$PC2)
sdH <- sd(H$PC2)
sdI <- sd(I$PC2)
sdJ <- sd(J$PC2)

(sdA + sdB + sdC + sdD + sdE + sdF + sdG + sdH + sdI + sdJ)/10

# 0.0203181

N = 244

PCRanges <- matrix(0, nrow = 20, ncol = 244)

for(i in 1:N){
  
  NumbersPC1 <- rnorm(n = 10, mean = Environment$PC1[i], sd = 0.03309501)
  NumbersPC2 <- rnorm(n = 10, mean = Environment$PC2[i], sd = 0.0203181)
  
  PCRanges[1:10,i] <- NumbersPC1
  PCRanges[11:20,i] <- NumbersPC2
  
  print(i)
  
}

write.csv(PCRanges, "PCRanges.csv")

# OUWie Test

# Model Testing

OUwieData <- read.csv("OUwieVariation.csv")
GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

N <- 100

answers <- data.frame()

t100 <- Phylogeny[1:100]
Tree <- t100[[1]]

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])

trait1 <- OUwieData$PC1A
names(trait1)<-OUwieData$Genus_species

trait2 <- OUwieData$PC1B
names(trait2)<-OUwieData$Genus_species

trait3 <- OUwieData$PC1C
names(trait3)<-OUwieData$Genus_species

trait4 <- OUwieData$PC1D
names(trait4)<-OUwieData$Genus_species

trait5 <- OUwieData$PC1E
names(trait5)<-OUwieData$Genus_species

trait6 <- OUwieData$PC1F
names(trait6)<-OUwieData$Genus_species

trait7 <- OUwieData$PC1G
names(trait7)<-OUwieData$Genus_species

trait8 <- OUwieData$PC1H
names(trait8)<-OUwieData$Genus_species

trait9 <- OUwieData$PC1I
names(trait9)<-OUwieData$Genus_species

trait10 <- OUwieData$PC1J
names(trait10)<-OUwieData$Genus_species

trait11 <- OUwieData$PC2A
names(trait11)<-OUwieData$Genus_species

trait12 <- OUwieData$PC2B
names(trait12)<-OUwieData$Genus_species

trait13 <- OUwieData$PC2C
names(trait13)<-OUwieData$Genus_species

trait14 <- OUwieData$PC2D
names(trait14)<-OUwieData$Genus_species

trait15 <- OUwieData$PC2E
names(trait15)<-OUwieData$Genus_species

trait16 <- OUwieData$PC2F
names(trait16)<-OUwieData$Genus_species

trait17 <- OUwieData$PC2G
names(trait17)<-OUwieData$Genus_species

trait18 <- OUwieData$PC2H
names(trait18)<-OUwieData$Genus_species

trait19 <- OUwieData$PC2I
names(trait19)<-OUwieData$Genus_species

trait20 <- OUwieData$PC2J
names(trait20)<-OUwieData$Genus_species

Fanning<-rep("Fanning", dim(OUwieData)[1])
Fanning[OUwieData$Fanning == "Fanning"] <- "Fanning"
Fanning[OUwieData$Fanning == "NonFanning"] <- "NonFanning"
Fanning <- as.factor(Fanning)
names(Fanning) <- OUwieData$Genus_species

for(k in 1:100){ 
  
  Tree <- t100[[k]]
  Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  fitER <- ace(Fanning, Tree, model = "ER", type = "discrete")
  fitARD <- ace(Fanning, Tree, model = "ARD", type = "discrete")
  
  states <- fitER
  
  if(fitARD$loglik - fitER$loglik > 3.8){states <- fitARD} 
  
  N <- length(Fanning)
  anc <- rep("NonFanning", N-1)
  names(anc) <- c((N+1):(2*N-1))
  anc[which(states$lik.anc[,1] > 0.5)] <- "Fanning"
  
  Tree$node.label <- anc
  
  trait1.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait1)
  
  trait2.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait2)
  
  trait3.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait3)
  
  trait4.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait4)
  
  trait5.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait5)
  
  trait6.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait6)
  
  trait7.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait7)
  
  trait8.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait8)
  
  trait9.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait9)
  
  trait10.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait10)
  
  trait11.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait11)
  
  trait12.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait12)
  
  trait13.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait13)
  
  trait14.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait14)
  
  trait15.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait15)
  
  trait16.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait16)
  
  trait17.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait17)
  
  trait18.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait18)
  
  trait19.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait19)
  
  trait20.data <- data.frame(Genus_species = OUwieData$Genus_species, Fanning=Fanning, X = trait20)
  
  answers[k,1] <- k
  
  test<-OUwie(Tree,trait1.data,model=c("BM1"))
  answers[k,2]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OU1"))
  answers[k,3]<-test$AICc
  
  test<-OUwie(Tree,trait1.data,model=c("OUM"))
  answers[k,4]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,5]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,6]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait3.data,model=c("BM1"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait3.data,model=c("OU1"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait3.data,model=c("OUM"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait4.data,model=c("BM1"))
  answers[k,11]<-test$AICc
  
  test<-OUwie(Tree,trait4.data,model=c("OU1"))
  answers[k,12]<-test$AICc
  
  test<-OUwie(Tree,trait4.data,model=c("OUM"))
  answers[k,13]<-test$AICc
  
  test<-OUwie(Tree,trait5.data,model=c("BM1"))
  answers[k,14]<-test$AICc
  
  test<-OUwie(Tree,trait5.data,model=c("OU1"))
  answers[k,15]<-test$AICc
  
  test<-OUwie(Tree,trait5.data,model=c("OUM"))
  answers[k,16]<-test$AICc
  
  test<-OUwie(Tree,trait6.data,model=c("BM1"))
  answers[k,17]<-test$AICc
  
  test<-OUwie(Tree,trait6.data,model=c("OU1"))
  answers[k,18]<-test$AICc
  
  test<-OUwie(Tree,trait6.data,model=c("OUM"))
  answers[k,19]<-test$AICc
  
  test<-OUwie(Tree,trait7.data,model=c("BM1"))
  answers[k,20]<-test$AICc
  
  test<-OUwie(Tree,trait7.data,model=c("OU1"))
  answers[k,21]<-test$AICc
  
  test<-OUwie(Tree,trait7.data,model=c("OUM"))
  answers[k,22]<-test$AICc
  
  test<-OUwie(Tree,trait8.data,model=c("BM1"))
  answers[k,23]<-test$AICc
  
  test<-OUwie(Tree,trait8.data,model=c("OU1"))
  answers[k,24]<-test$AICc
  
  test<-OUwie(Tree,trait8.data,model=c("OUM"))
  answers[k,25]<-test$AICc
  
  test<-OUwie(Tree,trait9.data,model=c("BM1"))
  answers[k,26]<-test$AICc
  
  test<-OUwie(Tree,trait9.data,model=c("OU1"))
  answers[k,27]<-test$AICc
  
  test<-OUwie(Tree,trait9.data,model=c("OUM"))
  answers[k,28]<-test$AICc
  
  test<-OUwie(Tree,trait10.data,model=c("BM1"))
  answers[k,29]<-test$AICc
  
  test<-OUwie(Tree,trait10.data,model=c("OU1"))
  answers[k,30]<-test$AICc
  
  test<-OUwie(Tree,trait10.data,model=c("OUM"))
  answers[k,31]<-test$AICc
  
  test<-OUwie(Tree,trait11.data,model=c("BM1"))
  answers[k,32]<-test$AICc
  
  test<-OUwie(Tree,trait11.data,model=c("OU1"))
  answers[k,33]<-test$AICc
  
  test<-OUwie(Tree,trait11.data,model=c("OUM"))
  answers[k,34]<-test$AICc
  
  test<-OUwie(Tree,trait12.data,model=c("BM1"))
  answers[k,35]<-test$AICc
  
  test<-OUwie(Tree,trait12.data,model=c("OU1"))
  answers[k,36]<-test$AICc
  
  test<-OUwie(Tree,trait12.data,model=c("OUM"))
  answers[k,37]<-test$AICc
  
  test<-OUwie(Tree,trait13.data,model=c("BM1"))
  answers[k,38]<-test$AICc
  
  test<-OUwie(Tree,trait13.data,model=c("OU1"))
  answers[k,39]<-test$AICc
  
  test<-OUwie(Tree,trait13.data,model=c("OUM"))
  answers[k,40]<-test$AICc
  
  test<-OUwie(Tree,trait14.data,model=c("BM1"))
  answers[k,41]<-test$AICc
  
  test<-OUwie(Tree,trait14.data,model=c("OU1"))
  answers[k,42]<-test$AICc
  
  test<-OUwie(Tree,trait14.data,model=c("OUM"))
  answers[k,43]<-test$AICc
  
  test<-OUwie(Tree,trait15.data,model=c("BM1"))
  answers[k,44]<-test$AICc
  
  test<-OUwie(Tree,trait15.data,model=c("OU1"))
  answers[k,45]<-test$AICc
  
  test<-OUwie(Tree,trait15.data,model=c("OUM"))
  answers[k,46]<-test$AICc
  
  test<-OUwie(Tree,trait16.data,model=c("BM1"))
  answers[k,47]<-test$AICc
  
  test<-OUwie(Tree,trait16.data,model=c("OU1"))
  answers[k,48]<-test$AICc
  
  test<-OUwie(Tree,trait16.data,model=c("OUM"))
  answers[k,49]<-test$AICc
  
  test<-OUwie(Tree,trait17.data,model=c("BM1"))
  answers[k,50]<-test$AICc
  
  test<-OUwie(Tree,trait17.data,model=c("OU1"))
  answers[k,51]<-test$AICc
  
  test<-OUwie(Tree,trait17.data,model=c("OUM"))
  answers[k,52]<-test$AICc
  
  test<-OUwie(Tree,trait18.data,model=c("BM1"))
  answers[k,53]<-test$AICc
  
  test<-OUwie(Tree,trait18.data,model=c("OU1"))
  answers[k,54]<-test$AICc
  
  test<-OUwie(Tree,trait18.data,model=c("OUM"))
  answers[k,55]<-test$AICc
  
  test<-OUwie(Tree,trait19.data,model=c("BM1"))
  answers[k,56]<-test$AICc
  
  test<-OUwie(Tree,trait19.data,model=c("OU1"))
  answers[k,57]<-test$AICc
  
  test<-OUwie(Tree,trait19.data,model=c("OUM"))
  answers[k,58]<-test$AICc
  
  test<-OUwie(Tree,trait20.data,model=c("BM1"))
  answers[k,59]<-test$AICc
  
  test<-OUwie(Tree,trait20.data,model=c("OU1"))
  answers[k,60]<-test$AICc
  
  test<-OUwie(Tree,trait20.data,model=c("OUM"))
  answers[k,61]<-test$AICc
  
  print(k)
}

write.csv(answers, file=paste("OuwieOutputsVariation.csv"))

#### KMult ####

# Fanning - Sample Size = 162

GroupNames <- read.csv("GroupNames.csv", header=FALSE)

LoadHarmonicArray <- function(filename){
  Array2D <- read.csv(file = filename, header = FALSE)
  nRow <- dim(Array2D)[1]
  nTaxa <- dim(Array2D)[2]
  nTaxa = nTaxa/2
  Array3D <- array(data = NA, c(nRow,2,nTaxa))
  for (i in 1:nTaxa){
    for (j in 1:nRow){
      Array3D[j,1,i] = Array2D[j, i*2 - 1]
      Array3D[j,2,i] = Array2D[j, i*2]
    }
  }
  return(Array3D)
}

Array <- LoadHarmonicArray("HarmonicArray.csv")

labels <- unlist(read.csv("GroupNames.csv", header=FALSE))
dimnames(Array)[[3]] <- labels

PSP.shape <- physignal(Array, Tree, iter=9999)
PSP.shape

## Fanning 1 - Sample Size = 244

GroupNames1 <- read.csv("GroupNames1.csv", header=FALSE)

LoadHarmonicArray <- function(filename){
  Array2D <- read.csv(file = filename, header = FALSE)
  nRow <- dim(Array2D)[1]
  nTaxa <- dim(Array2D)[2]
  nTaxa = nTaxa/2
  Array3D <- array(data = NA, c(nRow,2,nTaxa))
  for (i in 1:nTaxa){
    for (j in 1:nRow){
      Array3D[j,1,i] = Array2D[j, i*2 - 1]
      Array3D[j,2,i] = Array2D[j, i*2]
    }
  }
  return(Array3D)
}

Array <- LoadHarmonicArray("HarmonicArray1.csv")

labels <- unlist(read.csv("GroupNames1.csv", header=FALSE))
dimnames(Array)[[3]] <- labels

PSP.shape <- physignal(Array, Tree1, iter=9999)
PSP.shape

#### Ancestral State Reconstruction and Analyses #### 

## Fanning and Swimming - Sample Size = 162

GroupNames <- read.csv("GroupNames.csv", header=FALSE)

AncestralStatesContinuous <- read.csv("AncestralStatesContinuous.csv", row.names=1)
AncestralStatesPC1 <- as.matrix(AncestralStatesContinuous)[,1]
AncestralStatesPC2 <- as.matrix(AncestralStatesContinuous)[,2]

AncestralStatesDiscrete <- read.csv("AncestralStatesDiscrete.csv", row.names=1, header=FALSE)
AncestralStatesFanning <- as.matrix(AncestralStatesDiscrete)[,1]

AncestralStatesSwimming <- read.csv("AncestralStatesDiscrete.csv", row.names=1, header=FALSE)
if(sum(is.na(AncestralStatesSwimming$V3)>0)){AncestralStatesSwimming<-AncestralStatesSwimming[-which(is.na(AncestralStatesSwimming$V3)),]} 
AncestralStatesSwimming <- as.matrix(AncestralStatesSwimming)[,2]

# Phylogenetic Signal

phylosignal(AncestralStatesPC1, Tree, reps = 9999)
phylosignal(AncestralStatesPC2, Tree, reps = 9999)
phylosignal(AncestralStatesFanning, Tree, reps = 9999)
phylosignal(AncestralStatesSwimming, Tree, reps = 9999)

# PC1 Ancestral States

fitContinuous(Tree, AncestralStatesPC1, SE = 0, model = c("BM"))
fitContinuous(Tree, AncestralStatesPC1, SE = 0, model = c("OU"))
fitContinuous(Tree, AncestralStatesPC1, SE = 0, model = c("EB"))

FitPC1 <- fastAnc(Tree, AncestralStatesPC1, vars=TRUE, CI=TRUE)
FitPC1

FitPC1$CI[1,]
range(AncestralStatesPC1)

obj <- contMap(Tree, AncestralStatesPC1, plot=FALSE)
obj <- setMap(obj,viridisLite::viridis(n=8))
plot(obj, lwd = 3, fsize = 1.5, legend = 0.7*max(nodeHeights(Tree)))

# PC2 Ancestral States

fitContinuous(Tree, AncestralStatesPC2, SE = 0, model = c("BM"))
fitContinuous(Tree, AncestralStatesPC2, SE = 0, model = c("OU"))
fitContinuous(Tree, AncestralStatesPC2, SE = 0, model = c("EB"))

FitPC2 <- fastAnc(Tree, AncestralStatesPC2, vars=TRUE, CI=TRUE)
FitPC2

FitPC2$CI[1,]
range(AncestralStatesPC2)

obj <- contMap(Tree, AncestralStatesPC2, plot=FALSE)
obj <- setMap(obj,viridisLite::viridis(n=8))
plot(obj, lwd = 3, fsize = 1.5, legend = 0.7*max(nodeHeights(Tree)))

# Fanning Ancestral States

fitDiscrete(Tree, AncestralStatesFanning, model = c("ER"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("SYM"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("ARD"))

FitFanning <- ace(AncestralStatesDiscrete$V2, Tree, model = "ER", type = "discrete")
FitFanning
round(FitFanning$lik.anc, 3)

cols<-palette(c("black", "salmon3"))

plot.phylo(Tree, type = "fan", cex = 0.4, label.offset = 20, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V2))],
               sort(unique(AncestralStatesDiscrete$V2)))
par(lwd = 0.1)
nodelabels(node = 1:Tree$Nnode+Ntip(Tree), pie = FitFanning$lik.anc, 
           piecol = cols, cex = 0.15)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V2, sort(unique(AncestralStatesDiscrete$V2))), 
          piecol = cols, cex = 0.15)

# Swimming Ancestral States

fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("ER"))
fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("SYM"))
fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("ARD"))

AncestralStatesSwimming <- read.csv("AncestralStatesDiscrete.csv", row.names=1, header=FALSE)
if(sum(is.na(AncestralStatesSwimming$V3)>0)){AncestralStatesSwimming<-AncestralStatesSwimming[-which(is.na(AncestralStatesSwimming$V3)),]} 
AncestralStatesSwimming <- as.matrix(AncestralStatesSwimming)[,2]

FitSwimming <- ace(AncestralStatesDiscrete$V3, Tree, model = "SYM", type = "discrete")
FitSwimming
round(FitSwimming$lik.anc, 3)

cols<-palette(c("black", "salmon3", "white"))

plot.phylo(Tree, type = "fan", cex = 0.4, label.offset = 20, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V3))],
               sort(unique(AncestralStatesDiscrete$V3)))
par(lwd = 0.1)
nodelabels(node = 1:Tree$Nnode+Ntip(Tree), pie = FitSwimming$lik.anc, 
           piecol = cols, cex = 0.15)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V3, sort(unique(AncestralStatesDiscrete$V3))), 
          piecol = cols, cex = 0.15)

## Fanning and Swimming 1 - Sample Size = 244

GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

AncestralStatesContinuous <- read.csv("AncestralStatesContinuous1.csv", row.names=1, header=FALSE)
AncestralStatesPC1 <- as.matrix(AncestralStatesContinuous)[,1]
AncestralStatesPC2 <- as.matrix(AncestralStatesContinuous)[,2]

AncestralStatesDiscrete <- read.csv("AncestralStatesDiscrete1.csv", row.names=1, header=FALSE)
AncestralStatesFanning <- as.matrix(AncestralStatesDiscrete)[,1]

AncestralStatesSwimming <- read.csv("AncestralStatesDiscrete1.csv", row.names=1, header=FALSE)
if(sum(is.na(AncestralStatesSwimming$V3)>0)){AncestralStatesSwimming<-AncestralStatesSwimming[-which(is.na(AncestralStatesSwimming$V3)),]} 
AncestralStatesSwimming <- as.matrix(AncestralStatesSwimming)[,2]

# Phylogenetic Signal

phylosignal(AncestralStatesPC1, Tree1, reps = 9999)
phylosignal(AncestralStatesPC2, Tree1, reps = 9999)
phylosignal(AncestralStatesFanning, Tree1, reps = 999)
phylosignal(AncestralStatesSwimming, Tree1, reps = 999)

# PC1 Ancestral States

fitContinuous(Tree1, AncestralStatesPC1, SE = 0, model = c("BM"))
fitContinuous(Tree1, AncestralStatesPC1, SE = 0, model = c("OU"))
fitContinuous(Tree1, AncestralStatesPC1, SE = 0, model = c("EB"))

FitPC1 <- fastAnc(Tree1, AncestralStatesPC1, vars=TRUE, CI=TRUE)
FitPC1

FitPC1$CI[1,]
range(AncestralStatesPC1)

obj <- contMap(Tree1, AncestralStatesPC1, plot=FALSE)
obj <- setMap(obj,viridisLite::viridis(n=8))
plot(obj, lwd = 3, fsize = 1.5, legend = 0.7*max(nodeHeights(Tree1)))

# PC2 Ancestral States

fitContinuous(Tree1, AncestralStatesPC2, SE = 0, model = c("BM"))
fitContinuous(Tree1, AncestralStatesPC2, SE = 0, model = c("OU"))
fitContinuous(Tree1, AncestralStatesPC2, SE = 0, model = c("EB"))

FitPC2 <- fastAnc(Tree1, AncestralStatesPC2, vars=TRUE, CI=TRUE)
FitPC2

FitPC2$CI[1,]
range(AncestralStatesPC2)

obj <- contMap(Tree1, AncestralStatesPC2, plot=FALSE)
obj <- setMap(obj,viridisLite::viridis(n=8))
plot(obj, lwd = 3, fsize = 1.5, legend = 0.7*max(nodeHeights(Tree1)))

# Fanning Ancestral States

fitDiscrete(Tree1, AncestralStatesFanning, model = c("ER"))
fitDiscrete(Tree1, AncestralStatesFanning, model = c("SYM"))
fitDiscrete(Tree1, AncestralStatesFanning, model = c("ARD"))

FitFanning <- ace(AncestralStatesDiscrete$V2, Tree1, model = "ER", type = "discrete")
FitFanning
round(FitFanning$lik.anc, 3)

palette(c("black", "salmon3", "white"))
palette()

plot.phylo(Tree1, type = "fan", cex = 0.4, label.offset = 20, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V2))],
               sort(unique(AncestralStatesDiscrete$V2)))
par(lwd = 0.1)
nodelabels(node = 1:Tree1$Nnode+Ntip(Tree1), pie = FitFanning$lik.anc, 
           piecol = cols, cex = 0.15)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V2, sort(unique(AncestralStatesDiscrete$V2))), 
          piecol = cols, cex = 0.15)

# Swimming Ancestral States

fitDiscrete(Tree1, AncestralStatesSwimming, SE = 0, model = c("ER"))
fitDiscrete(Tree1, AncestralStatesSwimming, SE = 0, model = c("SYM"))
fitDiscrete(Tree1, AncestralStatesSwimming, SE = 0, model = c("ARD"))

FitSwimming <- ace(AncestralStatesDiscrete$V3, Tree1, model = "ARD", type = "discrete", use.eigen = FALSE, use.expm = TRUE)
FitSwimming
round(FitSwimming$lik.anc, 3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

palette(c("black", "salmon3", "transparent"))
palette()

plot.phylo(Tree1, type = "fan", cex = 0.4, label.offset = 20, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V3))],
               sort(unique(AncestralStatesDiscrete$V3)))
par(lwd = 0.1)
nodelabels(node = 1:Tree1$Nnode+Ntip(Tree1), pie = FitSwimming$lik.anc, 
           piecol = cols, cex = 0.15)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V3, sort(unique(AncestralStatesDiscrete$V3))), 
          piecol = cols, cex = 0.15)

#### BayesTraits ####

# Fanning - Sample Size = 162

Tree <- Phylogeny

BayesTraits <- read.csv("BayesTraits.csv", header=FALSE)

GroupNames <- read.csv("GroupNames.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Tree[[1]]$tip.label)

for(k in 1:100)
{
  Tree[[k]] <- keep.tip(Tree[[k]], Tree[[k]]$tip.label[match(GroupNames$V1, Tree[[k]]$tip.label)])
  print(k)
}

write.nexus(Tree, file="Tree.trees")

# Fanning - Sample Size = 244

Tree1 <- Phylogeny

BayesTraits <- read.csv("BayesTraits1.csv", header=FALSE)

GroupNames <- read.csv("GroupNames1.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% Tree1[[1]]$tip.label)

for(k in 1:100)
{
  Tree1[[k]] <- keep.tip(Tree1[[k]], Tree1[[k]]$tip.label[match(GroupNames$V1, Tree1[[k]]$tip.label)])
  print(k)
}

write.nexus(Tree1, file="Tree1.trees")

# Swimming - Sample Size = 162

TreeS <- Phylogeny

BayesTraits <- read.csv("BayesTraitsSPC.csv", header=FALSE)

GroupNames <- read.csv("GroupNamesS.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% TreeS[[1]]$tip.label)

for(k in 1:100)
{
  TreeS[[k]] <- keep.tip(TreeS[[k]], TreeS[[k]]$tip.label[match(GroupNames$V1, TreeS[[k]]$tip.label)])
  print(k)
}

write.nexus(TreeS, file="TreeS.trees")

# Swimming - Sample Size = 244

TreeS1 <- Phylogeny

BayesTraits <- read.csv("BayesTraitsSPC1.csv", header=FALSE)

GroupNames <- read.csv("GroupNamesS1.csv", header = F)

GroupNames <- subset(GroupNames, V1 %in% TreeS1[[1]]$tip.label)

for(k in 1:100)
{
  TreeS1[[k]] <- keep.tip(TreeS1[[k]], TreeS1[[k]]$tip.label[match(GroupNames$V1, TreeS1[[k]]$tip.label)])
  print(k)
}

write.nexus(TreeS1, file="TreeS1.trees")

#### Temperature and Oxygen PCA ####

PCA <- read.csv("PCA.csv")

Data.PCA <- princomp(PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

#### Environment Depth Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Depth)>0)){Environment<-Environment[-which(is.na(Environment$Depth)),]} 

Environment$zDepth<-scale(Environment$Depth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zDepth,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentDepth.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zDepth,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentDepth.Rdata")
  
}

save(Final.Model, file="EnvironmentDepth.Rdata")

#### Environment Oxygen Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen)),]} 

Environment$zOxygen<-scale(Environment$Oxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentOxygen.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentOxygen.Rdata")
  
}

save(Final.Model, file="EnvironmentOxygen.Rdata")

#### Environment Temperature Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[2]]

Environment <- read.csv("Environment.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Temperature)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zTemp<-scale(Environment$Temperature)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTemp.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTemp.Rdata")
  
}
save(Final.Model, file="EnvironmentTemp.Rdata")

#### Environment Temperature + Oxygen PC Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[2]]

Environment <- read.csv("Environment.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$PC1)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zPC1<-scale(Environment$PC1)
Environment$zPC2<-scale(Environment$PC2)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempOxyPC.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempOxyPC.Rdata")
  
}

save(Final.Model, file="EnvironmentTempOxyPC.Rdata")

#### Environment 1 Depth Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Depth)>0)){Environment<-Environment[-which(is.na(Environment$Depth)),]} 

Environment$zDepth<-scale(Environment$Depth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zDepth,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentDepth1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zDepth,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentDepth1.Rdata")
  
}

save(Final.Model, file="EnvironmentDepth1.Rdata")

#### Environment 1 Oxygen Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen)),]} 

Environment$zOxygen<-scale(Environment$Oxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentOxygen1.Rdata")
  
}

save(Final.Model, file="EnvironmentOxygen1.Rdata")

#### Environment 1 Temperature Modelling ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Temperature)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zTemp<-scale(Environment$Temperature)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTemp1.Rdata")
  
}
save(Final.Model, file="EnvironmentTemp1.Rdata")

#### Environment 1 Temperature + Oxygen PC Modelling Fanning ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("Environment1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$PC1)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zPC1<-scale(Environment$PC1)
Environment$zPC2<-scale(Environment$PC2)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempOxyPC1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempOxyPC1.Rdata")
  
}

save(Final.Model, file="EnvironmentTempOxyPC1.Rdata")

#### Environment Depth Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Depth)>0)){Environment<-Environment[-which(is.na(Environment$Depth)),]} 

Environment$zDepth<-scale(Environment$Depth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zDepth,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentDepthS.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zDepth,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentDepthS.Rdata")
  
}

save(Final.Model, file="EnvironmentDepthS.Rdata")

#### Environment Oxygen Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen)),]} 

Environment$zOxygen<-scale(Environment$Oxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zOxygen,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentOxygenS.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zOxygen,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentOxygenS.Rdata")
  
}

save(Final.Model, file="EnvironmentOxygenS.Rdata")

#### Environment Temperature Modelling Swimming  ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[2]]

Environment <- read.csv("EnvironmentS.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Temperature)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zTemp<-scale(Environment$Temperature)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zTemp,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempS.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zTemp,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempS.Rdata")
  
}
save(Final.Model, file="EnvironmentTempS.Rdata")

#### Environment Temperature + Oxygen PC Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$PC1)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zPC1<-scale(Environment$PC1)
Environment$zPC2<-scale(Environment$PC2)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempOxyPCS.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempOxyPCS.Rdata")
  
}

save(Final.Model, file="EnvironmentTempOxyPCS.Rdata")

#### Environment 1 Depth Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Depth)>0)){Environment<-Environment[-which(is.na(Environment$Depth)),]} 

Environment$zDepth<-scale(Environment$Depth)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zDepth,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zDepth,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentDepthS1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zDepth,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentDepthS1.Rdata")
  
}

save(Final.Model, file="EnvironmentDepthS1.Rdata")

#### Environment 1 Oxygen Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen)),]} 

Environment$zOxygen<-scale(Environment$Oxygen)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zOxygen,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentOxygenS1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zOxygen,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentOxygenS1.Rdata")
  
}

save(Final.Model, file="EnvironmentOxygenS1.Rdata")

#### Environment 1 Temperature Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[2]]

Environment <- read.csv("EnvironmentS1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Temperature)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zTemp<-scale(Environment$Temperature)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Swimming)~zTemp,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempS1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Swimming)~zTemp,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempS1.Rdata")
  
}
save(Final.Model, file="EnvironmentTempS1.Rdata")

#### Environment 1 Temperature + Oxygen PC Modelling Swimming ####

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentS1.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$PC1)>0)){Environment<-Environment[-which(is.na(Environment$Temperature)),]}

Environment$zPC1<-scale(Environment$PC1)
Environment$zPC2<-scale(Environment$PC2)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 3),
                             V = gelman.prior(~zPC1 + zPC2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentTempOxyPCS1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zPC1 + zPC2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentTempOxyPCS1.Rdata")
  
}

save(Final.Model, file="EnvironmentTempOxyPCS1.Rdata")

#### Environment Variability Oxygen Modelling Fanning ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/Data Repository/Generalised Linear Mixed Models/Variability Tests")

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentalData.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen1)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen1)),]} 

Environment$zOxygen1<-scale(Environment$Oxygen1)
Environment$zOxygen2<-scale(Environment$Oxygen2)
Environment$zOxygen3<-scale(Environment$Oxygen3)
Environment$zOxygen4<-scale(Environment$Oxygen4)
Environment$zOxygen5<-scale(Environment$Oxygen5)
Environment$zOxygen6<-scale(Environment$Oxygen6)
Environment$zOxygen7<-scale(Environment$Oxygen7)
Environment$zOxygen8<-scale(Environment$Oxygen8)
Environment$zOxygen9<-scale(Environment$Oxygen9)
Environment$zOxygen10<-scale(Environment$Oxygen10)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen1,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen1,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen1,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen1.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen1.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen2.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen2.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen2.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen3,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen3,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen3.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen3,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen3.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen3.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen4,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen4,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen4.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen4,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen4.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen4.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen5,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen5,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen5.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen5,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen5.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen5.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen6,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen6,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen6.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen6,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen6.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen6.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen7,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen7,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen7.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen7,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen7.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen7.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen8,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen8,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen8.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen8,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen8.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen8.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen9,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen9,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen9.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen9,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen9.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen9.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen10,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zOxygen10,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityOxygen10.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zOxygen10,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityOxygen10.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityOxygen10.Rdata")

summary(Final.Model)

# Summary Model

####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests")

load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen1.Rdata")
Final.Model1 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen2.Rdata")
Final.Model2 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen3.Rdata")
Final.Model3 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen4.Rdata")
Final.Model4 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen5.Rdata")
Final.Model5 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen6.Rdata")
Final.Model6 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen7.Rdata")
Final.Model7 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen8.Rdata")
Final.Model8 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen9.Rdata")
Final.Model9 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityOxygen10.Rdata")
Final.Model10 <- Final.Model

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentalData.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen1)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen1)),]} 

Environment$zOxygen<-scale(Environment$Oxygen1)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zOxygen,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Final.Model<-MCMCglmm(as.factor(Fanning)~zOxygen,
                      random = ~Species, 
                      ginverse = list(Species = animalA), 
                      prior = gelmanprior, 
                      verbose = TRUE, 
                      family = "categorical", 
                      data = Environment,
                      nitt = 10000,
                      thin = 1,
                      burnin = 0,
                      pl = TRUE, 
                      pr = TRUE, 
                      slice = TRUE) 

Final.Model$VCV[((1-1)*1000+1):(1*1000), ] <- Final.Model1$VCV[1:1000,] 
Final.Model$Sol[((1-1)*1000+1):(1*1000), ] <- Final.Model1$Sol[1:1000,] 
Final.Model$VCV[((2-1)*1000+1):(2*1000), ] <- Final.Model2$VCV[1:1000,] 
Final.Model$Sol[((2-1)*1000+1):(2*1000), ] <- Final.Model2$Sol[1:1000,] 
Final.Model$VCV[((3-1)*1000+1):(3*1000), ] <- Final.Model3$VCV[1:1000,] 
Final.Model$Sol[((3-1)*1000+1):(3*1000), ] <- Final.Model3$Sol[1:1000,] 
Final.Model$VCV[((4-1)*1000+1):(4*1000), ] <- Final.Model4$VCV[1:1000,] 
Final.Model$Sol[((4-1)*1000+1):(4*1000), ] <- Final.Model4$Sol[1:1000,] 
Final.Model$VCV[((5-1)*1000+1):(5*1000), ] <- Final.Model5$VCV[1:1000,] 
Final.Model$Sol[((5-1)*1000+1):(5*1000), ] <- Final.Model5$Sol[1:1000,] 
Final.Model$VCV[((6-1)*1000+1):(6*1000), ] <- Final.Model6$VCV[1:1000,] 
Final.Model$Sol[((6-1)*1000+1):(6*1000), ] <- Final.Model6$Sol[1:1000,] 
Final.Model$VCV[((7-1)*1000+1):(7*1000), ] <- Final.Model7$VCV[1:1000,] 
Final.Model$Sol[((7-1)*1000+1):(7*1000), ] <- Final.Model7$Sol[1:1000,] 
Final.Model$VCV[((8-1)*1000+1):(8*1000), ] <- Final.Model8$VCV[1:1000,] 
Final.Model$Sol[((8-1)*1000+1):(8*1000), ] <- Final.Model8$Sol[1:1000,] 
Final.Model$VCV[((9-1)*1000+1):(9*1000), ] <- Final.Model9$VCV[1:1000,] 
Final.Model$Sol[((9-1)*1000+1):(9*1000), ] <- Final.Model9$Sol[1:1000,] 
Final.Model$VCV[((10-1)*1000+1):(10*1000), ] <- Final.Model10$VCV[1:1000,] 
Final.Model$Sol[((10-1)*1000+1):(10*1000), ] <- Final.Model10$Sol[1:1000,] 

save(Final.Model, file = "EnvironmentVariabilityOxygen.Rdata")

summary(Final.Model)

#### Environment Variability Temperature Modelling Fanning ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/Data Repository/Generalised Linear Mixed Models/Variability Tests")

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentalData.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Temperature1)>0)){Environment<-Environment[-which(is.na(Environment$Temperature1)),]} 

Environment$zTemperature1<-scale(Environment$Temperature1)
Environment$zTemperature2<-scale(Environment$Temperature2)
Environment$zTemperature3<-scale(Environment$Temperature3)
Environment$zTemperature4<-scale(Environment$Temperature4)
Environment$zTemperature5<-scale(Environment$Temperature5)
Environment$zTemperature6<-scale(Environment$Temperature6)
Environment$zTemperature7<-scale(Environment$Temperature7)
Environment$zTemperature8<-scale(Environment$Temperature8)
Environment$zTemperature9<-scale(Environment$Temperature9)
Environment$zTemperature10<-scale(Environment$Temperature10)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp1,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp1,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp1.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp1,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp1.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp1.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp2,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp2,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp2.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp2,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp2.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp2.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp3,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp3,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp3.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp3,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp3.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp3.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp4,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp4,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp4.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp4,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp4.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp4.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp5,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp5,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp5.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp5,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp5.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp5.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp6,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp6,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp6.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp6,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp6.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp6.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp7,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp7,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp7.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp7,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp7.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp7.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp8,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp8,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp8.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp8,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp8.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp8.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp9,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp9,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp9.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp9,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp9.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp9.Rdata")

summary(Final.Model)

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp10,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Model<-MCMCglmm(as.factor(Fanning)~zTemp10,
                random = ~Species, 
                ginverse = list(Species = animalA), 
                prior = gelmanprior, 
                verbose = TRUE, 
                family = "categorical", 
                data = Environment,
                nitt = 11000,
                thin = 10,
                burnin = 1000,
                pl = TRUE, 
                pr = TRUE, 
                slice = TRUE) 

Final.Model<-Model
Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 

nsamp.l <- nrow(Model$VCV)
start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l, "Species"]))

save(Final.Model, file = "EnvironmentVariabilityTemp10.Rdata")

# Real Model

for(i in 1:100){ 
  Tree <- T100[[i]] 
  Tree <- collapse.singles(Tree)
  Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
  
  animalA <- inverseA(Tree)$Ainv 
  
  Model <- MCMCglmm(as.factor(Fanning)~zTemp10,
                    random = ~Species, 
                    ginverse=list(Species = animalA), 
                    prior = gelmanprior, 
                    verbose = FALSE, 
                    family = "categorical",  
                    start = start1.l,
                    data = Environment,
                    nitt = 22000,
                    thin = 2000,
                    burnin = 2000,
                    pl = TRUE,
                    pr = TRUE,
                    slice = TRUE)
  
  print(i)
  
  Final.Model$VCV[((i-1)*10+1):(i*10), ] <- Model$VCV[1:10,] 
  Final.Model$Sol[((i-1)*10+1):(i*10), ] <- Model$Sol[1:10,] 
  Final.Model$Liab[((i-1)*10+1):(i*10), ] <- Model$Liab[1:10,] 
  
  nsamp.l <- nrow(Model$VCV)
  start1.l = list(R = Model$VCV[nsamp.l, "units"], G = list(G1 = Model$VCV[nsamp.l,"Species"]))
  
  save(Final.Model, file = "EnvironmentVariabilityTemp10.Rdata")
  
}

save(Final.Model, file="EnvironmentVariabilityTemp10.Rdata")

summary(Final.Model)

# Summary Model

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests")

load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp1.Rdata")
Final.Model1 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp2.Rdata")
Final.Model2 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp3.Rdata")
Final.Model3 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp4.Rdata")
Final.Model4 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp5.Rdata")
Final.Model5 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp6.Rdata")
Final.Model6 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp7.Rdata")
Final.Model7 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp8.Rdata")
Final.Model8 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp9.Rdata")
Final.Model9 <- Final.Model
load("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Variability Tests/EnvironmentVariabilityTemp10.Rdata")
Final.Model10 <- Final.Model

Trees <- read.tree("Actinopterygii.trees")
T100 <- Trees 
Tree <- T100[[1]]

Environment <- read.csv("EnvironmentalData.csv")

if(sum(is.na(Environment$Species)>0)){Environment <- Environment[-which(is.na(Environment$Species)), ]}
if(sum(is.na(Environment$Oxygen1)>0)){Environment<-Environment[-which(is.na(Environment$Oxygen1)),]} 

Environment$zTemp<-scale(Environment$Temperature1)

T100 <- lapply(T100, drop.tip, tip=setdiff(Tree$tip.label, Environment$Species))

# Dummy Model

i=1

Tree <- T100[[i]]  
Tree <- collapse.singles(Tree)
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)
animalA <- inverseA(Tree)$Ainv

prior.gaussian <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

gelmanprior <- list(B = list(mu = rep(0, 2),
                             V = gelman.prior(~zTemp,
                                              data = Environment, scale=1+pi^2/3)), 
                    R = list(V = 1, fix = 1), G = list(G1 = list(V = 1E-10, nu = -1)))

Final.Model<-MCMCglmm(as.factor(Fanning)~zTemp,
                      random = ~Species, 
                      ginverse = list(Species = animalA), 
                      prior = gelmanprior, 
                      verbose = TRUE, 
                      family = "categorical", 
                      data = Environment,
                      nitt = 10000,
                      thin = 1,
                      burnin = 0,
                      pl = TRUE, 
                      pr = TRUE, 
                      slice = TRUE) 

Final.Model$VCV[((1-1)*1000+1):(1*1000), ] <- Final.Model1$VCV[1:1000,] 
Final.Model$Sol[((1-1)*1000+1):(1*1000), ] <- Final.Model1$Sol[1:1000,] 
Final.Model$VCV[((2-1)*1000+1):(2*1000), ] <- Final.Model2$VCV[1:1000,] 
Final.Model$Sol[((2-1)*1000+1):(2*1000), ] <- Final.Model2$Sol[1:1000,] 
Final.Model$VCV[((3-1)*1000+1):(3*1000), ] <- Final.Model3$VCV[1:1000,] 
Final.Model$Sol[((3-1)*1000+1):(3*1000), ] <- Final.Model3$Sol[1:1000,] 
Final.Model$VCV[((4-1)*1000+1):(4*1000), ] <- Final.Model4$VCV[1:1000,] 
Final.Model$Sol[((4-1)*1000+1):(4*1000), ] <- Final.Model4$Sol[1:1000,] 
Final.Model$VCV[((5-1)*1000+1):(5*1000), ] <- Final.Model5$VCV[1:1000,] 
Final.Model$Sol[((5-1)*1000+1):(5*1000), ] <- Final.Model5$Sol[1:1000,] 
Final.Model$VCV[((6-1)*1000+1):(6*1000), ] <- Final.Model6$VCV[1:1000,] 
Final.Model$Sol[((6-1)*1000+1):(6*1000), ] <- Final.Model6$Sol[1:1000,] 
Final.Model$VCV[((7-1)*1000+1):(7*1000), ] <- Final.Model7$VCV[1:1000,] 
Final.Model$Sol[((7-1)*1000+1):(7*1000), ] <- Final.Model7$Sol[1:1000,] 
Final.Model$VCV[((8-1)*1000+1):(8*1000), ] <- Final.Model8$VCV[1:1000,] 
Final.Model$Sol[((8-1)*1000+1):(8*1000), ] <- Final.Model8$Sol[1:1000,] 
Final.Model$VCV[((9-1)*1000+1):(9*1000), ] <- Final.Model9$VCV[1:1000,] 
Final.Model$Sol[((9-1)*1000+1):(9*1000), ] <- Final.Model9$Sol[1:1000,] 
Final.Model$VCV[((10-1)*1000+1):(10*1000), ] <- Final.Model10$VCV[1:1000,] 
Final.Model$Sol[((10-1)*1000+1):(10*1000), ] <- Final.Model10$Sol[1:1000,] 

save(Final.Model, file = "EnvironmentVariabilityTemp.Rdata")

summary(Final.Model)

#### OU Tree ####

GroupNames <- read.csv("FishData1.csv")

ConsensusTree <- read.nexus("ConsensusTree.nex")

ConsensusTree <- nnls.tree(method = "ultrametric", cophenetic(ConsensusTree), ConsensusTree)

min(ConsensusTree$edge.length)

plot(ConsensusTree)

GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

Tree <- drop.tip(ConsensusTree, ConsensusTree$tip.label[-match(GroupNames$V1, ConsensusTree$tip.label)])
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)

AncestralStatesContinuous <- read.csv("AncestralStatesContinuous1.csv", header=FALSE, row.names=1)
AncestralStatesPC1 <- as.matrix(AncestralStatesContinuous)[,1]
AncestralStatesPC2 <- as.matrix(AncestralStatesContinuous)[,2]

fit <- fitContinuous(ConsensusTree, AncestralStatesPC1, SE = 0, model = c("OU"))
alpha <- fit$opt$alpha

OUConsensusTree <- transformPhylo(ConsensusTree, model = "OU", alpha=alpha)

plotTree(OUConsensusTree)

# PGLS

Environment <- read.csv("Environment1.csv")

FishEnvironment <- caper::comparative.data(OUConsensusTree, Environment, Species, vcv=TRUE, vcv.dim=3)

TemperatureModel <- caper::pgls(Temperature ~ Fanning, FishEnvironment, lambda='ML')
summary(TemperatureModel)

OxygenModel <- caper::pgls(Oxygen ~ Fanning, FishEnvironment, lambda='ML')
summary(OxygenModel)

DepthModel <- caper::pgls(Depth ~ Fanning, FishEnvironment, lambda='ML')
summary(DepthModel)

TemperatureModelS <- caper::pgls(Temperature ~ Swimming, FishEnvironment, lambda='ML')
summary(TemperatureModelS)

OxygenModelS <- caper::pgls(Oxygen ~ Swimming, FishEnvironment, lambda='ML')
summary(OxygenModelS)

DepthModelS <- caper::pgls(Depth ~ Swimming, FishEnvironment, lambda='ML')
summary(DepthModelS)