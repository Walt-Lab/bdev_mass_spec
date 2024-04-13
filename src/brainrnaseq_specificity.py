import requests
import scipy

import numpy as np
import pandas as pd

from io import StringIO

cell_type_dict = {
    "microglia" : "microglla",
    "astrocyte" : "mature",
    "oligodendrocyte" : "oligodendrocyte",
    "neuron" : "neuron",
    "endothelial" : "endothelial"
}

def calculate_mean(df):
    """
    Calculates the mean of all numeric values in a row of a dataframe, and assigns to a new column called "Mean".
    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe with columns containing numeric values.
    """
    return df.assign(Mean=df.mean(axis=1, numeric_only=True))


def map_hgnc_ids(brain_rna_seq_raw_path):
    """
    Maps the HGNC IDs in the Brain RNA-Seq file to UniProt IDs.
    Parameters
    ----------
    brain_rna_seq_path: csv file path
        Path to the "homo sapiens.csv" file, downloaded from brainrnaseq.org.
    References
    ----------
    Zhang et al. (2016) Purification and characterization of progenitor and mature human astrocytes reveals transcriptional and functional differences with mouse. Neuron 89(1):37-53. PMID: 26687838.
    """
    hgnc_ids = (
        "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
    )
    brain_rna_seq = pd.read_csv(brain_rna_seq_raw_path)

    hgnc_uniprot_mapping_data = pd.read_csv(
        (StringIO(requests.get(hgnc_ids).text)),
        sep="\t",
        usecols=["hgnc_id", "uniprot_ids"],
    )

    hgnc_uniprot_mapping_data["uniprot_ids"] = hgnc_uniprot_mapping_data[
        "uniprot_ids"
    ].str.split("|")
    hgnc_uniprot_mapping_data = hgnc_uniprot_mapping_data.explode("uniprot_ids")
    hgnc_uniprot_mapping_data = hgnc_uniprot_mapping_data.reset_index(drop=True)

    brain_rna_seq = pd.merge(
        brain_rna_seq,
        hgnc_uniprot_mapping_data,
        left_on="id",
        right_on="hgnc_id",
        how="inner",
    )
    brain_rna_seq.dropna(subset=["uniprot_ids"], inplace=True)
    brain_rna_seq.drop_duplicates(subset=["uniprot_ids"], inplace=True)
    brain_rna_seq.set_index("uniprot_ids", inplace = True)

    return brain_rna_seq


def mean_cell_type(brain_rna_seq_data, cell_type):
    """
    Returns only the mean of the data for the specified cell type, as well as the UniProt ID information in an additional column
    Parameters
    ----------
    brain_rna_seq_data : pandas.DataFrame
        Dataframe with a column called "uniprot_ids" (contains UniProt ID), and other columns containing cell-type specific Brain RNA Seq data
    cell_type : {'astrocyte', 'endothelial', 'microglia', 'oligodendrocyte', 'neuron'}
        Cell type of interest requested. Options:
         - 'astrocyte': mature astrocytes
         - 'endothelial': endothelial cells
         - 'microglia': microglia cells
         - 'oligodendrocyte': oligodendrocytes
         - 'neuron': neurons
    """
    
    key = cell_type_dict[str(cell_type)]
    cell_type_df = brain_rna_seq_data.filter(like=str(key), axis = 1)
    means_df = calculate_mean(cell_type_df)
    return means_df.rename(
        columns={"uniprot_ids": "uniprot_ids", "Mean": key}
    )


def calculate_enrichment(row, specificity_metric):
    """
    Uses the numeric values in a row of a dataframe (with the row being expression data for a specific gene of interest) to determine a gene's specificity Returns a numpy array with an index of indentifier for the gene and values of specificity scores.
    Parameters
    ----------
    'row' : pandas.DataFrame
        Row of a DataFrame containing expression data for a gene of interest with one column for each cell type/tissue being analyzed. Index should be identifier for each gene.
    'specificity_metric' : {'tau', 'tsi', 'gini', 'hg', 'spm', 'zscore'}
        Method to determine specificity requested. Options:
            - tau: Tau specificity score
            - tsi: Tissue Specificity Index
            - gini: Gini coefficient
            - hg: Entropy of a gene's expression distribution
            - spm: Specificity Metric
            - zscore: Z-score
    References
    ----------
    Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform. 2017 Mar 1;18(2):205-214. doi: 10.1093/bib/bbw008. PMID: 26891983; PMCID: PMC5444245.
    Schug J, Schuller WP, Kappen C, Salbaum JM, Bucan M, Stoeckert CJ Jr. Promoter features related to tissue specificity as measured by Shannon entropy. Genome Biol. 2005;6(4):R33. doi: 10.1186/gb-2005-6-4-r33. Epub 2005 Mar 29. PMID: 15833120; PMCID: PMC1088961.
    Wright Muelas, M., Mughal, F., O’Hagan, S. et al. The role and robustness of the Gini coefficient as an unbiased tool for the selection of Gini genes for normalising expression profiling data. Sci Rep 9, 17960 (2019). https://doi.org/10.1038/s41598-019-54288-7.
    """
    row_array = np.array(row)
    if specificity_metric == "tau":
        row_x = row_array / max(row_array)
        return np.sum(1 - row_x) / ((len(row_x)) - 1)
    if specificity_metric == "tsi":
        return max(row_array) / sum(row_array)
    if specificity_metric == "gini":
        sorted_types = np.sort(row_array)
        cumulative_fraction_types = np.cumsum(sorted_types) / np.sum(sorted_types)
        cumulative_fraction_total = np.linspace(0, 1, len(sorted_types))
        area_under_line_of_perfect_equality = scipy.integrate.simps(
            cumulative_fraction_total, cumulative_fraction_total
        )
        area_under_lorenz_curve = scipy.integrate.simps(
            cumulative_fraction_types, cumulative_fraction_total
        )
        return area_under_line_of_perfect_equality / (
            area_under_line_of_perfect_equality + area_under_lorenz_curve
        )
    if specificity_metric == "hg":
        row_sum = np.sum(row_array)
        p_sub_i = row_array / row_sum
        return -1 * (sum(p_sub_i * np.log2(p_sub_i)))
    if specificity_metric == "spm":
        squared_array = row_array**2
        sum_squared_array = np.sum(squared_array)
        spm_score = squared_array / sum_squared_array
        return pd.Series(spm_score, index=row.index)
    if specificity_metric == "zscore":
        mean_array = np.mean(row_array)
        std_array = np.std(row_array)
        zscore_values = (row_array - mean_array) / std_array
        return pd.Series(zscore_values, index=row.index)
    
def create_enrichment_dataframe(brain_rna_seq_data):
    """
    Returns a dataframe containing the mean expression for each protein in the BrainRNA-Seq dataset.
    Parameters
    ----------
    brain_rna_seq_data: pandas.DataFrame
        Dataframe created from the raw BrainRNA-Seq dataset. Index should contain a unique identifier for each gene.
    """
    astrocytes = mean_cell_type(brain_rna_seq_data, "astrocyte")
    endothelial = mean_cell_type(brain_rna_seq_data, "endothelial")
    microglia = mean_cell_type(brain_rna_seq_data, "microglia")
    oligodendrocytes = mean_cell_type(brain_rna_seq_data, "oligodendrocyte")
    neurons = mean_cell_type(brain_rna_seq_data, "neuron")

    all_cell_types = pd.merge(
        pd.merge(
            pd.merge(
                pd.merge(astrocytes, endothelial, left_index=True, right_index = True),
                microglia,
                left_index=True, right_index = True,
            ),
            oligodendrocytes,
            left_index=True, right_index = True,
        ),
        neurons,
        left_index=True, right_index = True,
    )
    all_cell_types = all_cell_types[["mature", "microglla", "endothelial", "oligodendrocyte", "neuron"]]
    cell_type_dict_inverted = {v: k for k, v in cell_type_dict.items()}
    return all_cell_types.rename(columns = cell_type_dict_inverted)

def cell_type_enrichment(
    expression_df,
    cell_type,
    specificity_metric,
    specificity_cutoff,
):
    """
    Returns a list of UniProt IDs corresponding to targets that meet specified criteria to determine cell-type specificity.
    Parameters
    ----------
    'brain_rna_seq_data' : pandas.DataFrame
        Dataframe with a column called "uniprot_ids" (contains UniProt ID), and other columns containing cell-type specific Brain RNA Seq data for the cell types listed under cell_type
    'cell_type' : {'astrocyte', 'endothelial', 'microglia', 'oligodendrocyte', 'neuron'}
        Cell type of interest requested. Options:
         - 'astrocyte': mature astrocytes
         - 'endothelial': endothelial cells
         - 'microglia': microglia cells
         - 'oligodendrocyte': oligodendrocytes
         - 'neuron': neurons
    'specificity_metric': {'tsi', 'zscore', 'spm', 'tau', 'gini', 'hg'}
        Individualized metric of determining cell type specificity requested. Options:
        - 'tsi': tissue specificity index
        - 'zscore': z-score
        - 'spm': specificity measure
        - 'tau': tau index
        - 'gini' : gini coefficient
        - 'hg' : entropy of a gene's expression distribution
    'specificity_cutoff' : numeric
        Numeric value representing the minimum value of the second enrichment cutoff.
    References
    ----------
    Kryuchkova-Mostacci N, Robinson-Rechavi M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform. 2017 Mar 1;18(2):205-214. doi: 10.1093/bib/bbw008. PMID: 26891983; PMCID: PMC5444245.
    Schug J, Schuller WP, Kappen C, Salbaum JM, Bucan M, Stoeckert CJ Jr. Promoter features related to tissue specificity as measured by Shannon entropy. Genome Biol. 2005;6(4):R33. doi: 10.1186/gb-2005-6-4-r33. Epub 2005 Mar 29. PMID: 15833120; PMCID: PMC1088961.
    Wright Muelas, M., Mughal, F., O’Hagan, S. et al. The role and robustness of the Gini coefficient as an unbiased tool for the selection of Gini genes for normalising expression profiling data. Sci Rep 9, 17960 (2019). https://doi.org/10.1038/s41598-019-54288-7.
    """

    cell_type_uniprot_ids = []

    if specificity_metric == "enrichment":
        other_cell_types = expression_df.drop(cell_type, axis=1)
        cell_type_targets = expression_df[[cell_type]]
        cell_type_targets["comparison"] = other_cell_types.apply(
            lambda row: row.max(), axis=1
        )
        if cell_type != "all":
            for index, row in cell_type_targets.iterrows():
                cell_type_uniprot_ids.append(
                    row[
                        (row[cell_type]) > (row["comparison"] * specificity_cutoff)
                    ].index
                )
        if cell_type == "all":
            for index, row in cell_type_targets.iterrows():
                cell_type_uniprot_ids.append(
                    row[
                        (row[cell_type]) < (row["comparison"] * specificity_cutoff)
                    ].index
                )
    else:
        enrichment_values = expression_df.apply(
            lambda row: calculate_enrichment(row, specificity_metric), axis=1
        )
        if cell_type != "all":
            if specificity_metric == "zscore" or specificity_metric == "spm":
                cell_type_uniprot_ids = enrichment_values[
                    (enrichment_values[cell_type]) > (specificity_cutoff)
                ].index.tolist()
            if (
                specificity_metric == "tau"
                or specificity_metric == "gini"
                or specificity_metric == "hg"
                or specificity_metric == "tsi"
            ):
                cell_type_uniprots = enrichment_values[
                    enrichment_values > specificity_cutoff
                ].index.tolist()
                cell_type_df = expression_df.loc[cell_type_uniprots]
                for index, row in cell_type_df.iterrows():
                    max_column = row.idxmax()
                    if max_column == cell_type:
                        cell_type_uniprot_ids.append(index)

        if cell_type == "all":
            if specificity_metric == "tau" or "gini" or "hg" or "tsi":
                cell_type_uniprot_ids = enrichment_values[
                    enrichment_values < specificity_cutoff
                ].index.tolist()
            else:
                for index, row in enrichment_values.iterrows():
                    cell_type_uniprot_ids = enrichment_values[
                        (enrichment_values[cell_type]) > (specificity_cutoff)
                    ].index.tolist()

    return cell_type_uniprot_ids