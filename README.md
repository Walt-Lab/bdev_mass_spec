# Mass spectrometry analysis of cell-type specific brain-derived extracellular vesicle pulldowns

## Modules

### config.py
- Contains several global variables.

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