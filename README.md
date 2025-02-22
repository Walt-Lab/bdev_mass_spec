# OLINK Proteomic Analysis to Identify Potential Extracellular Vesicle-Associated Proteins

Analysis of OLINK proteomic data to identify proteins that may be associated with brain-derived extracellular vesicles.

### Key Features
- A dataset containing information concenrning 5416 unique proteins, collected via the OLINK HT panel using frationated human cerebrospinal fluid.
- Read OLINK parquet files and identify proteins that may be associated with extracellular vesicles using relative protein abundances in fractionated human cerebrospinal fluid.
- Overlay lists of proteins that may be associated with extracellular vesicles with single-cell RNA sequencing data and subcellular localization analysis to determine if a particular protein could be a potential cell-type specific immunocapture or validation target.

## Modules

### config.py
- Contains several global variables.

### raw_data_preprocessing.py
- Converts the raw parquet file produced by OLINK into a tidy dataframe.
- Generates graphs to display the median fractionation pattern of a protein of interest.
- Calculates the EV Association Score of a protein of interest.

Required Packages 
- matplotlib.axes
- matplotlib.pyplot
- pandas

Required Documentation
- config.py

### olink_fractionation.py
- Uses fractionation patterns reported by Olink to identify proteins that may be associated with extracellular vesicles.

Required Packages
- pandas

### specificity_functions.py
- Calculates various statistical measures of specificity, including tau score, tissue specificity index, gini coefficient, Shannon entropy, specificity measure, and zscore.

Required Packages
- numpy
- pandas
- scipy

### brainrnaseq_specificity.py
- Uses data collected and made available by BrainRNA-Seq to determine proteins that are specific to a cell type of interest.

Required Packages
- requests
- numpy
- pandas
- io
- pathlib

Required Documentation
- config
- specificity_functions

### gtex_specificity.py
- Uses data made available by GTEx to determine proteins that are specific to the brain.

Required Packages
- pandas

Required Documentation
- specificity_functions
- config

### deeptmhmm_localization.py
- Uses the DeepTMHMM deep learning model to identify the most likely subcellular localization of proteins of interest.

Required Packages
- biolib
- gzip
- os
- pathlib
- pandas

### identify_targets.py
- Uses specified fractionation, localization, and cell-type specificity criteria to identify protein targets.

Required Packages
- pandas
- typing
- pathlib

Required Documentation
- raw_data_preprocessing.py
- olink_fractionation.py
- brainrnaseq_specificity.py
- deeptmhmm_localization.py
- config.py