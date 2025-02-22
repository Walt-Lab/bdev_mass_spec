from io import StringIO
from pathlib import Path

import pandas as pd
import requests

from config import CELL_TYPES
from specificity_functions import calculate_enrichment


def calculate_mean(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the mean expression value across all columns for each row (gene).
    Adds a new column "Mean" containing the computed values.

    Parameters
    ----------
    df: pandas.DataFrame
        A DataFrame where rows represent genes and columns contain
        numeric expression values.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame with an additional "Mean" column containing the
        mean expression values.
    """

    return df.assign(Mean=df.mean(axis=1, numeric_only=True))


def process_hgnc_data(hgnc_ids: str) -> pd.DataFrame:
    hgnc_uniprot_mapping_data = pd.read_csv(
        StringIO(requests.get(hgnc_ids).text),
        sep="\t",
        usecols=[
            "hgnc_id",
            "uniprot_ids",
            "symbol",
            "name",
            "alias_symbol",
            "alias_name",
            "ensembl_gene_id",
        ],
    )
    hgnc_uniprot_mapping_data["uniprot_ids"] = hgnc_uniprot_mapping_data[
        "uniprot_ids"
    ].str.split("|")
    hgnc_uniprot_mapping_data = hgnc_uniprot_mapping_data.explode("uniprot_ids")
    return hgnc_uniprot_mapping_data.reset_index(drop=True)


def map_hgnc_ids(hgnc_ids: str, brain_rna_seq_raw_path: Path | str) -> pd.DataFrame:
    """
    Maps HGNC gene IDs from the Brain RNA-Seq dataset to UniProt protein IDs.
    This function integrates data from the Brain RNA-Seq dataset and
    HGNC ID mapping file.

    Parameters
    ----------
    hgnc_ids : str
        Path to a .txt file containing mappings from HGNC IDs to UniProt
        protein names.
    brain_rna_seq_raw_path : str
        Path to a .csv file containing raw RNA-Seq expression data.

    Returns
    -------
    pandas.DataFrame
        A DataFrame indexed by UniProt IDs, containing mapped gene
        expression data.

    Notes
    -----
    - The function merges the RNA-Seq dataset with HGNC mappings and
    removes duplicates.
    - Rows without UniProt IDs are dropped.
    References
    ----------
    Zhang et al. (2016) Purification and characterization of progenitor
    and mature human astrocytes reveals transcriptional and functional
    differences with mouse. Neuron 89(1):37-53. PMID: 26687838.
    """
    brain_rna_seq = pd.read_csv(brain_rna_seq_raw_path)
    hgnc_uniprot_mapping_data = process_hgnc_data(hgnc_ids)

    brain_rna_seq = (
        pd.merge(
            brain_rna_seq,
            hgnc_uniprot_mapping_data,
            left_on="id",
            right_on="hgnc_id",
            how="inner",
        )
        .dropna(subset=["uniprot_ids"])
        .drop_duplicates(subset=["uniprot_ids"])
        .set_index("uniprot_ids")
    )

    return brain_rna_seq


def mean_cell_type(brain_rna_seq_data: pd.DataFrame, cell_type: str) -> pd.DataFrame:
    """
    Computes the mean expression level of a specified cell type from the
    Brain RNA-Seq dataset.

    Parameters
    ----------
    brain_rna_seq_data : pandas.DataFrame
        A DataFrame where rows represent genes and columns contain
        expression levels for multiple cell types.
    cell_type : {'astrocyte', 'endothelial', 'microglia',
    'oligodendrocyte', 'neuron'}
        Cell type for which mean expression should be calculated.
        Available options:
        - 'astrocyte': Mature astrocytes
        - 'endothelial': Endothelial cells
        - 'microglia': Microglia cells
        - 'oligodendrocyte': Oligodendrocytes
        - 'neuron': Neurons

    Returns
    -------
    pandas.DataFrame
        A new DataFrame containing only the UniProt IDs and the mean
        expression level for the specified cell type.

    Notes
    -----
    - Uses the "CELL_TYPES" dictionary to account for naming
    inconsistencies in the raw dataset.
    """
    key = CELL_TYPES[str(cell_type)]
    cell_type_df = brain_rna_seq_data.filter(like=str(key), axis=1)
    means_df = calculate_mean(cell_type_df)
    return means_df.rename(columns={"uniprot_ids": "uniprot_ids", "Mean": key})


def create_enrichment_dataframe(brain_rna_seq_data: pd.DataFrame) -> pd.DataFrame:
    """
    Creates a summary DataFrame containing mean expression levels for each cell type.

    This function computes the mean expression of each gene for
    different brain cell types using the Brain RNA-Seq dataset. The
    resulting DataFrame is indexed by UniProt IDs and contains one
    column per cell type, with values representing the mean expression.

    Parameters
    ----------
    brain_rna_seq_data : pandas.DataFrame
        A DataFrame where rows represent genes (indexed by UniProt IDs)
        and columns contain cell-type-specific expression levels.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with UniProt IDs as the index and mean expression
        values for each cell type as columns. The column names are
        standardized to:
        - 'astrocyte': Mature astrocytes
        - 'endothelial': Endothelial cells
        - 'microglia': Microglia cells
        - 'oligodendrocyte': Oligodendrocytes
        - 'neuron': Neurons

    Notes
    -----
    - The function calls "mean_cell_type" for each cell type and merges
      the results.
    - The column names are adjusted using "CELL_TYPES" to match expected
      formatting.
    - The function ensures that all cell types are included in the final
      DataFrame.
    """

    cell_type_dfs = {
        cell_type: mean_cell_type(brain_rna_seq_data, cell_type)
        for cell_type in [
            "astrocyte",
            "endothelial",
            "microglia",
            "oligodendrocyte",
            "neuron",
        ]
    }

    all_cell_types = pd.concat(cell_type_dfs.values(), axis=1)
    all_cell_types_filtered = all_cell_types.filter([CELL_TYPES.values()])
    cell_type_dict_inverted = {v: k for k, v in CELL_TYPES.items()}
    return all_cell_types_filtered.rename(columns=cell_type_dict_inverted)


def cell_type_enrichment(
    expression_df: pd.DataFrame,
    cell_type: str,
    specificity_metric: str,
    specificity_cutoff: float,
) -> list[str]:
    """
    Identifies genes that meet a given specificity threshold for a particular cell type.

    Parameters
    ----------
    expression_df : pandas.DataFrame
        A DataFrame where rows represent genes and columns contain
        cell-type-specific expression data.
    cell_type : {'astrocyte', 'endothelial', 'microglia',
    'oligodendrocyte', 'neuron'}
        The cell type of interest.
    specificity_metric : {'enrichment', 'tsi', 'zscore', 'spm', 'tau',
    'gini', 'hg'}
        The metric used to determine specificity (see detailed
        descriptions in "calculate_enrichment").
    specificity_cutoff : float
        The minimum threshold value for considering a gene as specific
        to a given cell type.

    Returns
    -------
    list
        A list of UniProt IDs for genes that meet the specificity
        criteria.

    Notes
    -----
    - If "specificity_metric" is 'enrichment', genes must have an
      expression level **X times** higher than the highest non-target
      cell type.
    - If "specificity_metric" is 'zscore' or 'spm', genes must exceed
      "specificity_cutoff" directly.
    - If "specificity_metric" is one of {'tau', 'gini', 'hg', 'tsi'},
      genes must have the highest expression in the target cell type and
      exceed "specificity_cutoff".

    References
    ----------
    Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene
    expression tissue-specificity metrics. Brief Bioinform. 2017 Mar
    1;18(2):205-214. doi: 10.1093/bib/bbw008. PMID: 26891983; PMCID:
    PMC5444245.
    Schug J, Schuller WP, Kappen C, Salbaum JM, Bucan M, Stoeckert CJ
    Jr. Promoter features related to tissue specificity as measured by
    Shannon entropy. Genome Biol. 2005;6(4):R33. doi:
    10.1186/gb-2005-6-4-r33. Epub 2005 Mar 29. PMID: 15833120; PMCID:
    PMC1088961.
    Wright Muelas, M., Mughal, F., O’Hagan, S. et al. The role and
    robustness of the Gini coefficient as an unbiased tool for the
    selection of Gini genes for normalising expression profiling data.
    Sci Rep 9, 17960 (2019). https://doi.org/10.1038/s41598-019-54288-7.
    """
    cell_type_uniprot_ids = []

    if specificity_metric == "enrichment":
        other_cell_types = expression_df.drop(cell_type, axis=1)
        cell_type_targets = expression_df[[cell_type]]
        cell_type_targets["comparison"] = other_cell_types.apply(
            lambda row: row.max(), axis=1
        )
        if cell_type == "all":
            for index, row in cell_type_targets.iterrows():
                cell_type_uniprot_ids.append(
                    row[
                        (row[cell_type]) < (row["comparison"] * specificity_cutoff)
                    ].index
                )
        else:
            for index, row in cell_type_targets.iterrows():
                cell_type_uniprot_ids.append(
                    row[
                        (row[cell_type]) > (row["comparison"] * specificity_cutoff)
                    ].index
                )
    else:
        enrichment_values = expression_df.apply(
            lambda row: calculate_enrichment(row, specificity_metric), axis=1
        )
        if cell_type == "all":
            if specificity_metric in ["tau", "gini", "hg", "tsi"]:
                cell_type_uniprot_ids = enrichment_values[
                    enrichment_values < specificity_cutoff
                ].index.tolist()
            else:
                for index, row in enrichment_values.iterrows():
                    cell_type_uniprot_ids = enrichment_values[
                        (enrichment_values[cell_type]) > (specificity_cutoff)
                    ].index.tolist()
        else:
            if specificity_metric in ["tau", "gini", "hg", "tsi"]:
                cell_type_uniprots = enrichment_values[
                    enrichment_values > specificity_cutoff
                ].index.tolist()
                cell_type_df = expression_df.loc[cell_type_uniprots]
                for index, row in cell_type_df.iterrows():
                    max_column = row.idxmax()
                    if max_column == cell_type:
                        cell_type_uniprot_ids.append(index)
            if specificity_metric in ["zscore", "spm"]:
                cell_type_uniprot_ids = enrichment_values[
                    (enrichment_values[cell_type]) > (specificity_cutoff)
                ].index.tolist()
    return cell_type_uniprot_ids
