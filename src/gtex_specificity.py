import pandas as pd

import requests

from io import StringIO

from config import ORGANS
from specificity_functions import calculate_enrichment

def gtex_specificity(raw_file, specificity_metric = "tau"):
    df = pd.read_csv(raw_file, compression='gzip', sep='\t', skiprows=2)
    
    grouped_columns = {key: [] for key in ORGANS.keys()}

    # Group columns by organ name
    for col in df.columns:
        for group, organ_list in ORGANS.items():
            if any(organ in col for organ in organ_list):
                grouped_columns[group].append(col)

    # Create a new DataFrame with grouped data (e.g., mean expression per organ)
    grouped_data = {}

    for group, cols in grouped_columns.items():
        if cols:  
            grouped_data[group] = df[cols].median(axis=1)

    # Convert the grouped data into a new DataFrame
    grouped_df = pd.DataFrame(grouped_data)
    grouped_df.insert(0, "Name", df["Name"])  

    grouped_df.set_index("Name", inplace=True)
    specificity_scores = grouped_df.apply(lambda row: calculate_enrichment(row, specificity_metric), axis = 1)
    specificity_scores = specificity_scores.reset_index()

    specificity_scores["ensembl"] = specificity_scores["Name"].str.split('.').str[0]
    specificity_scores.drop(["Name"], axis = 1, inplace=True)

    specificity_scores.columns = ["body_tau_score", "ensembl_gene_id"]

    return specificity_scores