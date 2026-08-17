#### Morphospace ####
#### Setup ####

library(Morpho)
library(morphospace)
library(dplyr)
library(Momocs)
library(geomorph)

#### Morphospace Analyses ####

source("Morpho.R")

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
