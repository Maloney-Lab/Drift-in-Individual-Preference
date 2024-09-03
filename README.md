# Drift in Individual Preference
 Supporting code and analysis for Maloney et al paper investigating the change in behavioral preferences in Drosophila melanogaster and evaluating the potential ecological fitness effects of randomly shifting behavioral preferences.

 

| Data  | This folder contains raw data used in this paper|
|------|-------|
`serotonin_mutant_longform.csv`, `serotonin_drug_longform.csv` |	Handedness Data compiled across all experiments for each day, with number of turns and number of right turns.
Raw Circling Data | This folder contains the .mat files for all of the centroid data for each experiments in the circling data. This is the most compact form of the data (14.1 GB), but is transformed into analyzed data with the scripts `exportdatabyfly.m`
`AllChunks.csv`	| The real world environment chunks used for real world environment data. Used in `realworlddata.py`
`allsimresults.nc` |	Simulation results for all real world environmental simulations. Used in `AnalyzeRealWorldDataSims.ipynb`
`MatureAge_Run7_lin_matchinterestingemv` | Data from phase sweep (lin space) used for mature age analysis (Fig 3D. See `matureage.ipynb`
`MatureAge_Run3` |	Data from phase sweep (log space) used for mature age analysis. See `matureage.ipynb`


The following scripts and notebooks (MATLAB 2023b, Python 3.10) are available on [github](https://github.com/Maloney-Lab/Drift-in-Individual-Preference). 

**Scripts**	| **This folder contains the scripts used to analyze the data and run simulations**
---|----
`exportdatabyfly.m` |	script for converting matlab continuous circling data to python. Dependency: `AngleArrays.m`
`AngleArrays.m`|	helper function for analyzing circling data
`stanhelpers.py`|	Helper functions for Bayesian analysis with STAN in python. Used by `SerotoninStan.ipynb` and `DGRPstan.ipynb`
`loadcontinuousmatlabfiles.py`|	Helper functions for analyzing continuous data and transitioning from matlab to python
`continuousanalysis.py`|	Helper functions for analyzing continuous data
`dmodel6_AR_transformed.stan`|	STAN code for Bayesian Autoregressive Model. Used by `SerotoninStan.ipynb` and `DGRPStan.ipynb`
`simpledrift.py` |	Core simulation code for running individual simulations over time
`matrixmaker.py` |	Code for running multiple parameters at once to make phase space plots of parameters
`realworlddata.py` |	Code for running simulations on realworld data on the cluster
`realworlddata.sh` |	Helper shell script for above to run on server
`frequency_environment_server.py` |	Script for generating simulations on server across multiple environmental amplitudes and frequencies
`frequency_environment.sh` |	Helper shell script for above to run on server
`matureage_server.py` |	Script for generating simulations on server across multiple ages of reproductive maturities and frequencies
`matureage.sh` |	Helper shell script for above to run on server
**Notebooks** | **This folder contains notebooks used for data analysis and creating figures. Note: Paths may need to be edited for depend files and folders listed, as the data has been separated due to file size restrictions for repositories.**
`SerotoninStan.ipynb` |	Notebook for analyzing data from Serotonin Experiments in STAN and generating figures in paper. Dependencies: `serotonin_mutant_longform.csv`, `serotonin_drug_longform.csv`, `stanhelpers.py`
`AngleAnalysis.ipynb`| Notebook used for analyzing continuous circling data and respective figures/power analysis
`Analyze_Centroids.ipynb` |	Notebook used for lowpass filtered data (Fig 1B)
5`HT_Supplemental_Figures.ipynb` |	Notebook used for r values over time (Fig S1 H-M)
`DGRPStan.ipynb` | Code for analyzing DGRP data. Dependencies are `stanhelpers.py` and dmodel6`_AR_transformed.stan`
`freqvsenvmean.ipynb` |	Code for frequency vs envmean analysis. Note, the actual depencies (simulations based on randomly generated data) are too large (>100GB) to be included, however this file is maintained for clarity of sourcing. Data is created with `frequency_environment.sh` and `frequency_environment_server.py`
`AnalyzeRealWorldDataSims.ipynb` | Code for analyzing results of real world simulation data; depencies `allsimresults.nc`
