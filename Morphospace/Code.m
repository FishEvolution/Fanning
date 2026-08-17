% Fanning - Sample Size = 162

% Generate Morphospace of Fish Fins Showing Fanning and Non-Fanning

TDS = taxaSetFromFile("Fish.csv","FishNames.csv",150);

SS = TDS.theoMorph(1,2,500,0);

% Harmonic array generator for theoretical morphospace

TDS.GenerateHarmonicArray('HarmonicArray.csv')

% Robustness for morphospace

TDS.Robustness(1000,5) 

% Sensitivity analysis on one example wing picture. Multiple pictures were used.

LandmarkSensitivity("Acanthopagrus_australis.tif",10:10:100);

% Swimming - Sample Size = 162

% Generate Morphospace of Fish Fins Showing Fanning and Non-Fanning

TDSS = taxaSetFromFile("FishS.csv","FishNamesS.csv",150);

SSS = TDSS.theoMorph(1,2,500,0);

% Harmonic array generator for theoretical morphospace

TDSS.GenerateHarmonicArray('HarmonicArrayS.csv')

% Fanning - Sample Size = 244

% Generate Morphospace of Fish Fins Showing Fanning and Non-Fanning

TDS1 = taxaSetFromFile("Fish1.csv","FishNames.csv",150);

SS1 = TDS1.theoMorph(1,2,500,0);

% Harmonic array generator for theoretical morphospace

TDS1.GenerateHarmonicArray('HarmonicArray1.csv')

% Swimming - Sample Size = 244

% Generate Morphospace of Fish Fins Showing Fanning and Non-Fanning

TDSS1 = taxaSetFromFile("FishS1.csv","FishNamesS.csv",150);

SSS1 = TDSS1.theoMorph(1,2,500,0);

% Harmonic array generator for theoretical morphospace

TDSS1.GenerateHarmonicArray('HarmonicArrayS1.csv')

% Intraspecific Variation

% Generate Morphospace of Fish Fins Showing Fanning and Non-Fanning

TDSV = taxaSetFromFile("FishTest.csv","FishNamesTest.csv",150);

SSV = TDS.theoMorph(1,2,500,0);

% Harmonic array generator for theoretical morphospace

TDSV.GenerateHarmonicArray('HarmonicArrayVariation.csv')