#### OUwie ####
#### Setup ####

library(ape)
library(OUwie)
library(phangorn)
library(Morpho)
library(morphospace)
library(dplyr)
library(rptR)

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
