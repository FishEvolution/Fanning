READ ME - Fanning

Data Overview - Contains all the original data for Fanning and Swimming Method, Image Sources and Environmental Data used.

Below is a description of how each part of the methods were carried out include file names used in each analysis.

For all project the end of the file represents which part of the project it is used in:

None of the below - Fanning Method on Original Sample Size of 162
1 - Fanning Method on Increase Sample Size of 244
S - Swimming Method on Original Sample Size of 162
1S - Swimming Method on Increase Sample Size of 244

A) Images

1. Collect fish images that are facing the left, or flipping them if they are facing the right.
2. Produce an outline in Inkscape by tracing the pectoral fin using the pen tool to create a spiral path.
3. Convert the image to a binary image and then to an outline in ImageJ. You may need to invert the image to get a black outline on a white background.
4. Save these outlines as .tif files in the MATLAB folder.

B) Morphospace

1. Generate a taxon data set (TDS) and a morphospace in MATLAB using Fish.csv (Column 1 = Image.tif, Column 2 and 3 = all 0s, Column 4 = 1 for Non-Fanning, 2 for Fanning) and FishNames.csv (1st cell = Non-Fanning and 2nd = Fanning).
2. Generate a second TDS and morphospace using FishSwimming.csv and FishSwimmingNames.csv (1st cell = Non-Fanning, 2nd = Fanning, 3rd = Non-Fanning and Pectoral Swimming and 4th = Fanning and Pectoral Swimming) to allow for a comparison.
3. Perform sensitivity analyses on number of landmarks and number of samples for the first morphospace.
4. Produce a harmonic array from the TDS produced for the first morphospace.
5. Import the harmonic array into R (HarmonicArray.csv) and generate the PC values these values into FishData.csv (used for Fanning and Swimming).
6. Recreate morphospace using the package "morphospace".
7. Calculate the area of occupation of fanning and non-fanning fish in morphospace and compare this to the size of a random area of 81 species using the package "Momocs".
8. Determine the area of overlap between the fanning and non-fanning fish by measuring the area in InkScape and comparing this with the area of fanning and non-fanning using the known values calculated in R. 
9. Determine the shape transformation for PC1 and PC2 from negative to positive values using the packages "morphospace" and "geomorph".

C) Adaptive Landscape

1. Create an OUwie.csv file containing the "Genus_species" names, the "Fanning" states and the "PC1" and "PC2" values.
2. Import the multiple trees file (Actinopterygii.trees)
3. Run the code to test for multiple evolutionary models using the packages "OUwie" and "phytools"
4. Collect the AICc in the file OuwieOutputs.csv.
5. Re-run the code on only the best model for each PC axes and collect the optimality values in the file OuwiePCOutputs to determine the optimality peaks for each PC axis.

D) Phylomorphospace

1. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
2. For all analyses prune the tree(s) needed using the appropriate GroupNames.csv file.
3. Ensure names in the phylogenetic tree, MATLAB Fish.csv file and R GroupNames.csv file match and are in the same order as the tree using the "Check Names" section.
4. Test for phylogenetic signal (Kmult) using the package "geomorph".

E) Ancestral State Reconstruction and Analyses

1. Create two ancestral state files with no headings, one called AncestralStateContinuous.csv which contains all the species names and the PC1 and PC2 values, and one called AncestralStateDiscrete.csv which contains all the species names and the Fanning and Swimming trait data, with traits represented as 1s (Non-Fanning and Non-Pectoral Swimming), 2s (Fanning and Pectoral Swimming) and any missing data put as a 3.
2. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree for the analyses.
3. Prune the Tree to only include species found in the GroupNames.csv file using the package "phangorn".
4. Ensure the tree is ultrametric.
5. Determine phylogenetic signal for the traits using the package "picante".
6. Use FitContinuous for PC values and FitDiscrete for fanning and swimming behaviour using the package "geiger" to test the different evolutionary models.
7. Perform an ancestral state reconstruction using the above models using the package "ape" and plot these on a phylogeny.

F) BayesTraits

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

G) Bayesian Phylogenetic Generalised Linear Mixed Models

1. Import the 100 phylogenetic trees (Actinopterygii.trees) into R using the package "ape" and use one tree at random for the analyses.
2. Import Environment data (Environment.csv) and remove any rows that do not contain data for the variable being tested (Temperature, Oxygen Concentration and Depth separately).
3. Perform a MCMCglmm analysis on the data and save the model as "Environment[Variable].Rdata"
4. Also perform a VIF analysis on temperature and oxygen concentration to determine collinearity using the VIF function in the package "regclass". 
