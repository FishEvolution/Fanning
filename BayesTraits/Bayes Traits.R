#### Bayes Traits ####
#### Setup ####

library(ape)

Phylogeny <- read.tree("Actinopterygii.trees")

#### Tree Files ####

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