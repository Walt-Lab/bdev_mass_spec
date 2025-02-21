import pandas as pd

from config import ORGANS
from specificity_functions import calculate_enrichment


def gtex_specificity(raw_file, specificity_metric="tau") -> pd.DataFrame:
    """
    Computes tissue specificity scores from GTEx RNA-seq data.

    Parameters
    ----------
    raw_file: str
        Path to the GTEx RNA-seq data file (gzipped TSV format).
    specificity_metric: {'tau', 'tsi', 'gini', 'hg'}
        Metric used to calculate enrichment (default: "tau"). Options:
        - 'tau': Tau specificity score
        - 'tsi': Tissue Specificity Index
        - 'gini': Gini coefficient
        - 'hg': Shannon entropy of gene expression distribution

    Returns
    ----------
    DataFrame
        A DataFrame containing columns for tissue specificity scores and Ensembl gene IDs.
    """
    df = pd.read_csv(raw_file, compression="gzip", sep="\t", skiprows=2)

    grouped_columns = {key: [] for key in ORGANS.keys()}

    for col in df.columns:
        for group, organ_list in ORGANS.items():
            if any(organ in col for organ in organ_list):
                grouped_columns[group].append(col)

    grouped_data = {}

    for group, cols in grouped_columns.items():
        if cols:
            grouped_data[group] = df[cols].median(axis=1)

    grouped_df = pd.DataFrame(grouped_data)
    grouped_df.insert(0, "Name", df["Name"])

    grouped_df.set_index("Name", inplace=True)
    specificity_scores = grouped_df.apply(
        lambda row: calculate_enrichment(row, specificity_metric), axis=1
    )
    specificity_scores = specificity_scores.reset_index()

    specificity_scores["ensembl"] = specificity_scores["Name"].str.split(".").str[0]
    specificity_scores.drop(["Name"], axis=1, inplace=True)

    specificity_scores.columns = ["body_tau_score", "ensembl_gene_id"]

    return specificity_scores
