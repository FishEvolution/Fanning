#### Ancestral State Reconstruction and Analyses #### 
#### Setup ####

library(ape)
library(phangorn)
library(picante)
library(geiger)
library(geomorph)

#### Consensus Trees ####

# Sample Size = 162

Tree <- read.nexus("ConsensusTree.nex")

Tree <- nnls.tree(method = "ultrametric", cophenetic(Tree), Tree)

min(Tree$edge.length)

plot(Tree)

# Sample Size = 244

Tree1 <- read.nexus("ConsensusTree1.nex")

Tree1 <- nnls.tree(method = "ultrametric", cophenetic(Tree1), Tree1)

min(Tree1$edge.length)

plot(Tree1)

#### Fanning and Swimming - Sample Size = 162 ####

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

#### Fanning and Swimming 1 - Sample Size = 244 ####

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