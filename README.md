# OLINK Proteomic Analysis to Identify Potential Extracellular Vesicle-Associated Proteins

Analysis of OLINK proteomic data to identify proteins that may be associated with brain-derived extracellular vesicles.

### Key Features
- A dataset containing information concenrning 5416 unique proteins, collected via the OLINK HT panel using frationated human cerebrospinal fluid.
- Read OLINK parquet files and identify proteins that may be associated with extracellular vesicles using relative protein abundances in fractionated human cerebrospinal fluid.
- Overlay lists of proteins that may be associated with extracellular vesicles with single-cell RNA sequencing data and subcellular localization analysis to determine if a particular protein could be a potential cell-type specific immunocapture or validation target.

## Modules
### raw_data_preprocessing.py
- Converts the raw parquet file produced by OLINK into a tidy dataframe.
- Generates graphs to display the median fractionation pattern of a protein of interest.
- Calculates the EV Association Score of a protein of interest.

Required Packages 
- matplotlib.pyplot
- pandas
### olink_fractionation.py
- Uses fractionation patterns to identify proteins that may be associated with extracellular vesicles.

Required Packages
- pandas
### brainrnaseq_specificity.py
- Uses data collected and made available by BrainRNA-Seq to determine proteins that are specific to a cell type of interest.

Required Packages
- requests
- scipy
- numpy
- pandas
- io
### deeptmhmm_localization.py
- Uses the DeepTMHMM deep learning model to identify the most likely subcellular localization of proteins of interest.

Required Packages
- biolib
- gzip
- os
- pathlib
- pandas