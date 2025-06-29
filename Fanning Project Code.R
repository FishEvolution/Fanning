#### Setup ####

setwd("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R")

library(rfishbase)
library(morphospace, lib="C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R")
library(dplyr)
library(Momocs)
library(geomorph)
library(ape)
library(phangorn)
library(OUwie, lib="C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R")
library(picante)
library(ggtree) #Only if Non-Circular Used
library(ggplot2) #Only if Non-Circular Used

#### Swimming Method Data Collection ####

Swimming <- read.csv("Species.csv")

SwimmingMethod <- swimming(species_list=Swimming$FishBase.Species)

write.csv(SwimmingMethod, 'SwimmingMethod.csv')

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
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, method = "vector", mag = 0.5)
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
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, method = "vector", mag = 0.5)
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
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, method = "vector", mag = 0.5)
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
plotRefToTarget(M1 = Extshps_Neg1, Extshps_Pos1, method = "vector", mag = 0.5)
points(Extshps_Neg1)

# PC2 Shape Transformation

Extshps2 <- extract_shapes(mspace = msp, axis = 2, nshapes = 2)
Extshps_Neg2 <- inv_efourier(Extshps2$shapes[1,])
Extshps_Pos2 <- inv_efourier(Extshps2$shapes[2,])
plot(Extshps_Pos2, asp = 1, type = "l", col = "red")
points(Extshps_Neg2, type = "l", col = "blue")
dev.new(width=1, height=1)
plotRefToTarget(M1 = Extshps_Neg2, Extshps_Pos2, method = "vector", mag = 0.5)
points(Extshps_Neg2)

#### Adaptive Landscape ####

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
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,11]<-test$AICc
  
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
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,11]<-test$AICc
  
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
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,5]<-test$AICc
  answers[k,6]<-test$theta[1,1]
  
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
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,11]<-test$AICc
  
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
  
  test<-OUwie(Tree,trait2.data,model=c("BM1"))
  answers[k,7]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("BMS"))
  answers[k,8]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OU1"))
  answers[k,9]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUM"))
  answers[k,10]<-test$AICc
  
  test<-OUwie(Tree,trait2.data,model=c("OUMV"))
  answers[k,11]<-test$AICc
  
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

#### KMult ####

# Import Tree

Phylogeny <- read.tree("Actinopterygii.trees")
T100 <- Phylogeny 
Tree <- T100[[1]]

## Fanning - Sample Size = 162

GroupNames <- read.csv("GroupNames.csv", header=FALSE)
GroupNames$V1 %in% Tree$tip.label

PrunedTree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
write.tree(PrunedTree)

write.csv(PrunedTree$tip.label, "TreeNames.csv")

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

PSP.shape <- physignal(Array, PrunedTree, iter=9999)
PSP.shape

## Swimming - Sample Size = 162

GroupNames <- read.csv("GroupNamesS.csv", header=FALSE)
GroupNames$V1 %in% Tree$tip.label

PrunedTree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
write.tree(PrunedTree)

write.csv(PrunedTree$tip.label, "TreeNamesS.csv")

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

Array <- LoadHarmonicArray("HarmonicArrayS.csv")

labels <- unlist(read.csv("GroupNamesS.csv", header=FALSE))
dimnames(Array)[[3]] <- labels

PSP.shape <- physignal(Array, PrunedTree, iter=9999)
PSP.shape

## Fanning 1 - Sample Size = 244

GroupNames1 <- read.csv("GroupNames1.csv", header=FALSE)
GroupNames1$V1 %in% Tree$tip.label

PrunedTree1 <- drop.tip(Tree, Tree$tip.label[-match(GroupNames1$V1, Tree$tip.label)])
write.tree(PrunedTree1)

write.csv(PrunedTree1$tip.label, "TreeNames1.csv")

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

PSP.shape <- physignal(Array, PrunedTree1, iter=9999)
PSP.shape

## Swimming 1 - Sample Size = 244

GroupNames <- read.csv("GroupNamesS1.csv", header=FALSE)
GroupNames$V1 %in% Tree$tip.label

PrunedTree1 <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
write.tree(PrunedTree1)

write.csv(PrunedTree1$tip.label, "TreeNames1S.csv")

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

Array <- LoadHarmonicArray("HarmonicArrayS1.csv")

labels <- unlist(read.csv("GroupNamesS1.csv", header=FALSE))
dimnames(Array)[[3]] <- labels

PSP.shape <- physignal(Array, PrunedTree1, iter=9999)
PSP.shape

#### Ancestral State Reconstruction and Analyses #### 

# Setup

Phylogeny <- read.tree("Actinopterygii.trees")
T100 <- Phylogeny 
Tree <- T100[[1]]

## Fanning and Swimming - Sample Size = 162

GroupNames <- read.csv("GroupNames.csv", header=FALSE)

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)

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

phenogram(Tree, AncestralStatesPC1, fsize=0.6, spread.costs = c(1,0))

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

phenogram(Tree,AncestralStatesPC2,fsize=0.6,spread.costs=c(1,0))

# Fanning Ancestral States

fitDiscrete(Tree, AncestralStatesFanning, model = c("ER"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("SYM"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("ARD"))

Nodes <- read.csv("Nodes.csv", header=FALSE)

N <- 34 

Answers <- rep(0,N)

for(i in 1:N){
  Species <- c(Nodes$V1[i], Nodes$V2[i])
  if (Species[1] == Species[2]) {
    Species <- c(Nodes$V1[i])
  }
  Node <- MRCA(Tree, Species)
  Answers[i] <- Node
  print(i)
}

write.csv(Answers, file=paste("Answers.csv"))

FitFanning <- ace(AncestralStatesDiscrete$V2, Tree, model = "ER", type = "discrete")
FitFanning
round(FitFanning$lik.anc, 3)

ancstats <- as.data.frame(FitFanning$lik.anc)
ancstats$node <- 1:Tree$Nnode+Ntip(Tree)

cols<-palette(c("black", "salmon3"))

pies <- nodepie(ancstats, cols = 1:2, outline.color = "black")
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

Nodes <- read.csv("Nodes.csv", header=FALSE)

clade_df <- data.frame(node = Nodes$V4, clade = Nodes$V3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggTree <- ggtree(Tree) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(50, 300, 50, 50),bgcolor="transparent",fgcolor="transparent") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggTree <- ggTree +
  geom_cladelab(data = clade_df,mapping = aes(node = node, label = clade),fontsize = 7, barsize = 2.5, extend = 0.3, offset = 4, align = TRUE)

ggTree <- ggTree + geom_strip("Girella_nigricans", "Scorpaena_scrofa", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Oncorhynchus_kisutch", "Cyclothone_microdon", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Acanthurus_nigroris", "Miichthys_miiuy", barsize=  2.5, extend = 0.3, offset = 4)

ggTree <- ggTree + geom_inset(pies, width = 0.035, height = 0.035) 

ggTree <- gheatmap(ggTree, data.frame(factor(AncestralStatesFanning)), width=0.02, low="black", high="salmon3", colnames = FALSE, color = "black") +
  scale_fill_manual(values=c("black","salmon3")) +
  theme(legend.position = 'none')

ggTree

ggsave("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Tree.pdf", width = 50, height = 80, units = "cm", limitsize = FALSE)

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

ancstats <- as.data.frame(FitSwimming$lik.anc)
ancstats$node <- 1:Tree$Nnode+Ntip(Tree)

cols<-palette(c("black", "salmon3", "white"))

pies <- nodepie(ancstats, cols = 1:3, outline.color = "black")
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

Nodes <- read.csv("Nodes.csv", header=FALSE)

clade_df <- data.frame(node = Nodes$V4, clade = Nodes$V3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggTree <- ggtree(Tree) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(50, 300, 50, 50),bgcolor="transparent",fgcolor="transparent") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggTree <- ggTree +
  geom_cladelab(data = clade_df,mapping = aes(node = node, label = clade),fontsize = 7, barsize = 2.5, extend = 0.3, offset = 4, align = TRUE)

ggTree <- ggTree + geom_strip("Girella_nigricans", "Scorpaena_scrofa", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Oncorhynchus_kisutch", "Cyclothone_microdon", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Acanthurus_nigroris", "Miichthys_miiuy", barsize=  2.5, extend = 0.3, offset = 4)

ggTree <- ggTree + geom_inset(pies, width = 0.035, height = 0.035) 

ggTree <- gheatmap(ggTree, data.frame(factor(AncestralStatesSwimming)), width=0.02, low="black", high="white", colnames = FALSE, color = "black") +
  scale_fill_manual(values=c("black","salmon3","white")) +
  theme(legend.position = 'none')

ggTree

## Fanning and Swimming 1 - Sample Size = 244

GroupNames <- read.csv("GroupNames1.csv", header=FALSE)

Tree <- drop.tip(Tree, Tree$tip.label[-match(GroupNames$V1, Tree$tip.label)])
Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)

AncestralStatesContinuous <- read.csv("AncestralStatesContinuous1.csv", row.names=1, header=FALSE)
AncestralStatesPC1 <- as.matrix(AncestralStatesContinuous)[,1]
AncestralStatesPC2 <- as.matrix(AncestralStatesContinuous)[,2]

AncestralStatesDiscrete <- read.csv("AncestralStatesDiscrete1.csv", row.names=1, header=FALSE)
AncestralStatesFanning <- as.matrix(AncestralStatesDiscrete)[,1]

AncestralStatesSwimming <- read.csv("AncestralStatesDiscrete1.csv", row.names=1, header=FALSE)
if(sum(is.na(AncestralStatesSwimming$V3)>0)){AncestralStatesSwimming<-AncestralStatesSwimming[-which(is.na(AncestralStatesSwimming$V3)),]} 
AncestralStatesSwimming <- as.matrix(AncestralStatesSwimming)[,2]

# Phylogenetic Signal

phylosignal(AncestralStatesPC1, Tree, reps = 9999)
phylosignal(AncestralStatesPC2, Tree, reps = 9999)
phylosignal(AncestralStatesFanning, Tree, reps = 999)
phylosignal(AncestralStatesSwimming, Tree, reps = 999)

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

phenogram(Tree, AncestralStatesPC1, fsize=0.6, spread.costs = c(1,0))

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

phenogram(Tree,AncestralStatesPC2,fsize=0.6,spread.costs=c(1,0))

# Fanning Ancestral States

fitDiscrete(Tree, AncestralStatesFanning, model = c("ER"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("SYM"))
fitDiscrete(Tree, AncestralStatesFanning, model = c("ARD"))

FitFanning <- ace(AncestralStatesDiscrete$V2, Tree, model = "ER", type = "discrete")
FitFanning
round(FitFanning$lik.anc, 3)

palette(c("black", "salmon3", "white"))
palette()

plot.phylo(Tree, type = "fan", cex = 0.7, label.offset = 10, 
           no.margin = TRUE, x.lim = 100)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V2))],
               sort(unique(AncestralStatesDiscrete$V2)))
par(lwd = 0.1)
nodelabels(node = 1:Tree$Nnode+Ntip(Tree), pie = FitFanning$lik.anc, 
           piecol = cols, cex = 0.2)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V2, sort(unique(AncestralStatesDiscrete$V2))), 
          piecol = cols, cex = 0.2)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9*par()$usr[1],
                  y = -max(nodeHeights(Tree)), fsize=0.8)

Nodes <- read.csv("Nodes1.csv", header=FALSE)

N <- 40 

Answers <- rep(0,N)

for(i in 1:N){
  Species <- c(Nodes$V1[i], Nodes$V2[i])
  if (Species[1] == Species[2]) {
    Species <- c(Nodes$V1[i])
  }
  Node <- MRCA(Tree, Species)
  Answers[i] <- Node
  print(i)
}

write.csv(Answers, file=paste("Answers1.csv"))

FitFanning <- ace(AncestralStatesDiscrete$V2, Tree, model = "ER", type = "discrete")
FitFanning
round(FitFanning$lik.anc, 3)

ancstats <- as.data.frame(FitFanning$lik.anc)
ancstats$node <- 1:Tree$Nnode+Ntip(Tree)

cols<-palette(c("black", "salmon3"))

pies <- nodepie(ancstats, cols = 1:2, outline.color = "black")
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

Nodes <- read.csv("Nodes1.csv", header=FALSE)

clade_df <- data.frame(node = Nodes$V4, clade = Nodes$V3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggTree <- ggtree(Tree) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(50, 300, 50, 50),bgcolor="transparent",fgcolor="transparent") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggTree <- ggTree +
  geom_cladelab(data = clade_df,mapping = aes(node = node, label = clade),fontsize = 7, barsize = 2.5, extend = 0.3, offset = 4, align = TRUE)

ggTree <- ggTree + geom_strip("Acanthurus_nigroris", "Lutjanus_russellii", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Istiophorus_platypterus", "Sphyraena_barracuda", barsize = 2.5, extend = 0.3, offset = 4)

ggTree <- ggTree + geom_inset(pies, width = 0.035, height = 0.035) 

ggTree <- gheatmap(ggTree, data.frame(factor(AncestralStatesFanning)), width=0.02, low="black", high="salmon3", colnames = FALSE, color = "black") +
  scale_fill_manual(values=c("black","salmon3")) +
  theme(legend.position = 'none')

ggTree

ggsave("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Tree.pdf", width = 50, height = 80, units = "cm", limitsize = FALSE)

# Swimming Ancestral States

fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("ER"))
fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("SYM"))
fitDiscrete(Tree, AncestralStatesSwimming, SE = 0, model = c("ARD"))

FitSwimming <- ace(AncestralStatesDiscrete$V3, Tree, model = "ARD", type = "discrete", use.eigen = FALSE, use.expm = TRUE)
FitSwimming
round(FitSwimming$lik.anc, 3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

palette(c("black", "salmon3", "transparent"))
palette()

plot.phylo(Tree, type = "fan", cex = 0.7, label.offset = 10, 
           no.margin = TRUE)

cols<-setNames(palette()[1:length(unique(AncestralStatesDiscrete$V3))],
               sort(unique(AncestralStatesDiscrete$V3)))
par(lwd = 0.1)
nodelabels(node = 1:Tree$Nnode+Ntip(Tree), pie = FitSwimming$lik.anc, 
           piecol = cols, cex = 0.2)
tiplabels(pie=to.matrix(AncestralStatesDiscrete$V3, sort(unique(AncestralStatesDiscrete$V3))), 
          piecol = cols, cex = 0.2)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9*par()$usr[1],
                  y = -max(nodeHeights(Tree)), fsize=0.8)


FitSwimming <- ace(AncestralStatesDiscrete$V3, Tree, model = "ARD", type = "discrete", use.eigen = FALSE, use.expm = TRUE)
FitSwimming
round(FitSwimming$lik.anc, 3)

ancstats <- as.data.frame(FitSwimming$lik.anc)
ancstats$node <- 1:Tree$Nnode+Ntip(Tree)

cols<-palette(c("black", "salmon3", "white"))

pies <- nodepie(ancstats, cols = 1:3, outline.color = "black")
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols))

Nodes <- read.csv("Nodes1.csv", header=FALSE)

clade_df <- data.frame(node = Nodes$V4, clade = Nodes$V3)

mycol <- rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

ggTree <- ggtree(Tree) + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(50, 300, 50, 50),bgcolor="transparent",fgcolor="transparent") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggTree <- ggTree +
  geom_cladelab(data = clade_df,mapping = aes(node = node, label = clade),fontsize = 7, barsize = 2.5, extend = 0.3, offset = 4, align = TRUE)

ggTree <- ggTree + geom_strip("Acanthurus_nigroris", "Lutjanus_russellii", barsize = 2.5, extend = 0.3, offset = 4) +
  geom_strip("Istiophorus_platypterus", "Sphyraena_barracuda", barsize = 2.5, extend = 0.3, offset = 4)

ggTree <- ggTree + geom_inset(pies, width = 0.035, height = 0.035) 

ggTree <- gheatmap(ggTree, data.frame(factor(AncestralStatesSwimming)), width=0.02, low="black", high="salmon3", colnames = FALSE, color = "black") +
  scale_fill_manual(values=c("black","salmon3","white")) +
  theme(legend.position = 'none')

ggTree

ggsave("C:/Users/tmunr/OneDrive/Aberdeen/Projects/Fanning Fin Project/R/Tree.pdf", width = 50, height = 80, units = "cm", limitsize = FALSE)

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

PCA <- read.csv("PCA.csv")

Data.PCA <- princomp(PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

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

PCA <- read.csv("PCA.csv")

Data.PCA <- princomp(PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

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

PCA <- read.csv("PCA.csv")

Data.PCA <- princomp(PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

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

PCA <- read.csv("PCA.csv")

Data.PCA <- princomp(PCA)

write.csv(Data.PCA$scores, file = "PCAScores.csv")

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
