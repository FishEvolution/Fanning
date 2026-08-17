#### Generalised Linear Mixed Models ####
#### Setup ####

library(ape)
library(regclass)
library(phangorn)
library(MCMCglmm)

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
