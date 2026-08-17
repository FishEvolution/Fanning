#### OU Tree ####
#### Setup ####

library(ape)
library(phangorn)
library(geiger)
library(motmot)
library(caper)

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