# gtex_specificity

ORGANS = {
    "Neural": ["Brain", "Pituitary", "Nerve"],
    "Heart": ["Heart"],
    "Small_Intestine": ["Small_Intestine"],
    "Colon": ["Colon"],
    "Liver": ["Liver"],
    "Pancreas": ["Pancreas"],
    "Esophagus": ["Esophagus"],
    "Stomach": ["Stomach"],
    "Kidney": ["Kidney"],
    "Adipose": ["Adipose"],
    "Artery": ["Artery"],
    "Skin": ["Skin"],
    "Muscle": ["Muscle"],
    "Cervix": ["Cervix"],
    "Lung": ["Lung"],
    "Spleen": ["Spleen"],
    "Testis": ["Testis"],
    "Ovary": ["Ovary"],
    "Prostate": ["Prostate"],
    "Thyroid": ["Thyroid"],
    "Bladder": ["Bladder"],
    "Uterus": ["Uterus"],
    "Vagina": ["Vagina"],
    "Breast": ["Breast"],
    "Cells": ["Cells"],
    "Fallopian_Tube": ["Fallopian_Tube"],
    "Minor_Salivary_Gland": ["Minor_Salivary_Gland"],
    "Whole_Blood": ["Whole_Blood"],
    "Adrenal_Gland": ["Adrenal_Gland"],
}

# brainrnaseq_specificity

hgnc_ids = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

CELL_TYPES = {
    "microglia": "microglla",  # mirrors inconsistency in raw data file
    "astrocyte": "mature",
    "oligodendrocyte": "oligodendrocyte",
    "neuron": "neuron",
    "endothelial": "endothelial",
}

hgnc_ids = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
