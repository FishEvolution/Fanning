# Fanning

READ ME

Data Overview.xlsx - Contains all the original data for Fanning and Swimming Method, Image Sources and Environmental Data used.

Fanning Project Code.R - Overall R script file for all analyses. Separate scripts for each individual analyses are found within their respect folders.

Tree File - Due to size it is obtainable from Rabosky et al. 2018 (Rabosky DL, Chang J, Cowman PF, Sallan L, Friedman M, Kaschner K, Garilao C, Near TJ, Coll M, Alfaro ME. An inverse latitudinal gradient in speciation rate for marine fishes. Nature. 2018 Jul;559(7714):392-5.)

Below is a description of how each part of the methods were carried out including folder and file names used in each analysis. Further descriptions of the files are found within their respective folder in their own ReadMe Files.

A) Folder - Morphospace 

Images (All Images are found as links in Data Overview.xlsx under the Image column.

1. Collect fish images that are facing the left, or flipping them if they are facing the right.
2. Produce an outline in Inkscape by tracing the pectoral fin using the pen tool to create a spiral path.
3. Convert the image to a binary image and then to an outline in ImageJ. You may need to invert the image to get a black outline on a white background.
4. Save these outlines as .tif files in the MATLAB folder.

Morphospace

1. Generate a taxon data set (TDS) and a morphospace in MATLAB using Fish.csv (Column 1 = Image.tif, Column 2 and 3 = all 0s, Column 4 = 1 for Non-Fanning, 2 for Fanning) and FishNames.csv (1st cell = Non-Fanning and 2nd = Fanning).
2. Generate a second TDS and morphospace using FishSwimming.csv and FishSwimmingNames.csv (1st cell = Non-Fanning, 2nd = Fanning, 3rd = Non-Fanning and Pectoral Swimming and 4th = Fanning and Pectoral Swimming) to allow for a comparison.
3. Perform sensitivity analyses on number of landmarks and number of samples for the first morphospace.
4. Produce a harmonic array from the TDS produced for the first morphospace.
5. Import the harmonic array into R (HarmonicArray.csv) and generate the PC values these values into FishData.csv (used for Fanning and Swimming).
6. Recreate morphospace using the package "morphospace".
7. Calculate the area of occupation of fanning and non-fanning fish in morphospace and compare this to the size of a random area of 81 species using the package "Momocs".
8. Determine the area of overlap between the fanning and non-fanning fish by measuring the area in InkScape and comparing this with the area of fanning and non-fanning using the known values calculated in R. 
9. Determine the shape transformation for PC1 and PC2 from negative to positive values using the packages "morphospace" and "geomorph".

B) Folder - OUwie

1. Create an OUwie.csv file containing the "Genus_species" names, the "Fanning" states and the "PC1" and "PC2" values.
2. Import the multiple trees file (Actinopterygii.trees)
3. Run the code to test for multiple evolutionary models using the packages "OUwie" and "phytools"
4. Collect the AICc in the file OuwieOutputs.csv.
5. Re-run the code on only the best model for each PC axes and collect the optimality values in the file OuwiePCOutputs to determine the optimality peaks for each PC axis.

C) Folder - Ancestral State Reconstruction

1. Create two ancestral state files with no headings, one called AncestralStateContinuous.csv which contains all the species names and the PC1 and PC2 values, and one called AncestralStateDiscrete.csv which contains all the species names and the Fanning and Swimming trait data, with traits represented as 1s (Non-Fanning and Non-Pectoral Swimming), 2s (Fanning and Pectoral Swimming) and any missing data put as a 3.
2. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree for the analyses.
3. Prune the Tree to only include species found in the GroupNames.csv file using the package "phangorn".
4. Ensure the tree is ultrametric.
5. Create a consensus tree using BEAST.
6. Determine phylogenetic signal for the traits using the package "picante".
7. Use FitContinuous for PC values and FitDiscrete for fanning and swimming behaviour using the package "geiger" to test the different evolutionary models.
8. Perform an ancestral state reconstruction using the above models using the package "ape" and plot these on a phylogeny.

Kmult

1. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
2. For all analyses prune the tree(s) needed using the appropriate GroupNames.csv file.
3. Ensure names in the phylogenetic tree, MATLAB Fish.csv file and R GroupNames.csv file match and are in the same order as the tree using the "Check Names" section.
4. Test for phylogenetic signal (Kmult) using the package "geomorph".

D) Folder - BayesTraits

1. Copy species names, fanning behaviour (coded as 0 for Non-Fanning and 1 for Fanning) and swimming behaviour (coded as 0 for Non-Pectoral Swimming and 1 for Pectoral Swimming) from FishData.csv into a BayesTraits.csv and .txt file.
2. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape".
3. Produce a Pruned Tree of only species in the BayesTrait.csv file using the package "ape" and export it.
4. Remove all text from after the word #NEXUS until the BEGIN TREES; line from the .trees file.
5. Place the .txt file and .trees file in the BayesTraits folder.
6. Open the command window and the set directory to the BayesTraits folder using the cd command.
7. Run the following code: BayesTraitsV4.exe Tree.trees FanningSwimming.txt for the discrete dependent model or BayesTraitsV4.exe Tree.trees FanningSwimmingI.txt for discrete independent model.
8. Choose from the following options: 
	2 / 3 (2 and 3 to compare dependent and independent models respectively)
	2 (For MCMC)
	PriorAll exp 10 (If numbers are very large)
	EqualTrees 10000 (If one tree is being focused on)
	Stones 100 1000 (For comparing dependent and independent models)
	run 
9. Repeat steps 7 for Fanning behaviour and PC1 value. using PC.txt for the dependent model and PCI for the independent model, turning PC1 into a binary value of 0 being < 0.01 and 1 being > 0.01.
10. Choose from the following options: 
	2 / 3 (2 and 3 to compare dependent and independent models respectively)
	2 (For MCMC)
	PriorAll exp 10 (If numbers are very large)
	EqualTrees 10000 (If one tree is being focused on)
	Stones 100 1000 (For comparing dependent and independent models)
	run 
11. Models:
	BayesTraits.csv = Fanning and Swimming, sample size = 162
	BayesTraits1.csv = Fanning and Swimming, sample size = 244
	BayesTraitsFPC.csv = Fanning and PC1, sample size = 162
	BayesTraitsFPC1.csv = Fanning and PC1, sample size = 244
	BayesTraitsSPC.csv = Swimming and PC1, sample size = 64
	BayesTraitsSPC1.csv = Swimming and PC1, sample size = 146

E) Folder - Generalised Linear Mixed Models

1. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
2. Import Environment data (Environment.csv) and remove any rows that do not contain data for the variable being tested (Temperature, Oxygen Concentration and Depth separately).
3. Perform a MCMCglmm analysis on the data and save the model as "Environment[Variable].Rdata"
4. Also perform a VIF analysis on temperature and oxygen concentration to determine collinearity using the VIF function in the package "regclass". 

F) Folder - PGLS

1. Import the Consensus tree (ConsensusTree.nex).
2. Transform the tree using an Ornstein-Uhlenbeck transformation.
3. Important the Environmental data for all 244 species (Envrionment1.csv).
4. Perform a PGLS analysis using the pgls function in the package "caper" for temperature, oxygen and depth.
