import sys

sys.path.append(
    "C:\\Users\\Wyss User\\AppData\\Local\\Packages\\PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0\\LocalCache\\local-packages\\Python311\\site-packages"
)

import biolib
import gzip
import math
import os
import re
import requests
import pathlib

import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests as r
import xml.etree.ElementTree as ET

from io import StringIO


brain_rna_seq_path = (
    "C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\fe-wp-dataset-124.csv"
)
raw_data = "C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\LCSET-29370_AShah_MNorman_Extended_NPX_2024-02-14.parquet"
uniprot_fasta_database = (
    "C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\uniprot_fasta_database.gz"
)
plate_layout_path = "C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\Plate Layout.xlsx"
assay_list_path = (
    "C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\ht_panel_assay_list.xlsx"
)


def parse_gz_file(file_path):
    """
    Creates a dictionary of UniProt IDs and their corresponding FASTA sequences using a .gz file.
    Parameters
    ----------
    file_path: path to a .gz file containing UniProt IDs and FASTA sequences.
    """
    protein_dict = {}
    current_uniprot_id = None
    current_sequence = ""
    with gzip.open(file_path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_uniprot_id is not None:
                    protein_dict[current_uniprot_id] = current_sequence
                    current_sequence = ""
                if "|" in line:
                    current_uniprot_id = line.split("|")[1].strip()
                else:
                    print(f"Skipping line without expected format: {line}")
                    current_uniprot_id = None
            else:
                current_sequence += line
        if current_uniprot_id is not None:
            protein_dict[current_uniprot_id] = current_sequence
    return protein_dict


def tmhmm_localization(assays, output_directory):
    """
    Uses DeepTMHMM to characterize the localization of each amino acid in a protein.
    Parameters
    ----------
    assays : pandas.DataFrame
        DataFrame with columns called 'UniProt ID' and 'Sequence'
    """
    deeptmhmm = biolib.load("DTU/DeepTMHMM")
    assay_list = []
    with open("query.fasta", "w") as fasta_file:
        for _, row in assays.iterrows():
            fasta_line = f">{row['UniProt ID']}\n{row['Sequence']}\n"
            fasta_file.write(fasta_line)
            assay_list.append(fasta_line)
    deeptmhmm_job = deeptmhmm.cli(args="--fasta query.fasta", machine="local")
    deeptmhmm_job.save_files(output_directory)


def identify_localization(assays, region, output_directory="ht_output"):
    """
    Identifies the localization of proteins using the output of DeepTMHMM.
    Parameters
    ----------
    assays : pandas.DataFrame
        DataFrame with columns called 'UniProt ID' and 'Sequence'
    region : {'TMhelix', 'inside', 'outside', 'internal', 'external'}
        Subcellular region requested. Options:
          - 'TMhelix': transmembrane proteins
          - 'inside': at least some of the protein is inside the cell/EV
          - 'outside': at least some of the protein is outside the cell/EV
          - 'internal': the protein is only found inside the cell, no transmembrane or outside domains
          - 'external': the protein is only found outside the cell, no transmembrane or inside domains
    output_directory: directory path
        Path to a directory in which the localization data will be stored.
    """
    output_directory_path = pathlib.Path(output_directory)
    if not os.path.exists(output_directory_path):
        os.makedirs(output_directory_path)
    if not os.path.exists(output_directory_path / "TMRs.gff3"):
        tmhmm_localization(assays, output_directory_path)
    localization_df = pd.read_csv(
        output_directory_path / "TMRs.gff3",
        sep="\t",
        comment="#",
        names=[
            "uniprot_id",
            "region_location",
            "region_start",
            "region_end",
            0,
            1,
            2,
            3,
        ],
    )
    localization_df = localization_df[localization_df["uniprot_id"] != "//"].dropna(
        axis=1
    )
    get_regional_uniprots = lambda region: set(
        localization_df[localization_df["region_location"] == region]["uniprot_id"]
    )
    if region == "internal":
        tm_uniprots = get_regional_uniprots("TMhelix")
        outside_uniprots = get_regional_uniprots("outside")
        inside_uniprots = get_regional_uniprots("inside")
        return inside_uniprots - tm_uniprots - outside_uniprots
    if region == "external":
        tm_uniprots = get_regional_uniprots("TMhelix")
        outside_uniprots = get_regional_uniprots("outside")
        inside_uniprots = get_regional_uniprots("inside")
        return outside_uniprots - inside_uniprots - tm_uniprots
    else:
        return get_regional_uniprots(region)
    

def clean_strings(strings):
    cleaned_strings = []
    for string in strings:
        cleaned_string = (
            str(string)
            .replace("(", "")
            .replace(")", "")
            .replace("'", "")
            .replace(",", "")
        )
        cleaned_strings.append(cleaned_string)
    return cleaned_strings


def clean_up_raw_data(raw_data, plate_layout_path):
    ht_data = pd.read_parquet(raw_data)
    plate_layout = pd.read_excel(plate_layout_path)
    data = pd.merge(ht_data, plate_layout, how="left", on="SampleID")
    # repeat assays: P32455, Q02750
    replicate_assays = data[["UniProt"]].value_counts()
    replicate_assays = replicate_assays[replicate_assays > 94]
    replicate_assay_list = clean_strings(replicate_assays.index.tolist())
    unique_data = data[~data["UniProt"].isin(replicate_assay_list)]
    ctrl_dict = {}
    for assay in list(unique_data["Assay"].unique()):
        df = unique_data[unique_data["Assay"] == assay]
        neg_ctrl = df[df["SampleType"] == "NEGATIVE_CONTROL"][
            "PCNormalizedNPX"
        ].median()
        ctrl_dict[assay] = neg_ctrl
    unique_data["Delta"] = unique_data.apply(
        lambda row: (
            row["PCNormalizedNPX"] - ctrl_dict[row["Assay"]]
            if row["PCNormalizedNPX"] >= ctrl_dict[row["Assay"]]
            else ctrl_dict[row["Assay"]]
        ),
        axis=1,
    )
    unique_data.loc[:, "Linear NPX"] = unique_data["PCNormalizedNPX"].map(lambda x: 2**x)
    tidy_data = unique_data[unique_data["SampleType"] == "SAMPLE"].pivot(
        columns="UniProt",
        index=["SampleID", "Health", "Sample", "CSF_sample"],
        values="Linear NPX",
    )
    return tidy_data


def calculate_mean(df):
    """
    Calculates the mean of all numeric values in a row of a dataframe, and assigns to a new column called "Mean".
    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe with columns containing numeric values.
    """
    return df.assign(Mean=df.mean(axis=1, numeric_only=True))


def map_hgnc_ids(brain_rna_seq_path):
    """
    Maps the HGNC IDs in the Brain RNA-Seq file to UniProt IDs.
    Parameters
    ----------
    brain_rna_seq_path: csv file path
        Path to the "homo sapiens.csv" file, downloaded from brainrnaseq.org.
    """
    hgnc_ids = (
        "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
    )
    brain_rna_seq = pd.read_csv(brain_rna_seq_path)

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
    if cell_type == "microglia":
        microglia_df = calculate_mean(
            brain_rna_seq_data[
                brain_rna_seq_data.filter(like="microglla").columns.append(
                    pd.Index(["uniprot_ids"])
                )
            ]
        )[["uniprot_ids", "Mean"]]
        return microglia_df.rename(
            columns={"uniprot_ids": "uniprot_ids", "Mean": "microglia"}
        )
    if cell_type == "astrocyte":
        astrocyte_df = calculate_mean(
            brain_rna_seq_data[
                brain_rna_seq_data.filter(like="mature").columns.append(
                    pd.Index(["uniprot_ids"])
                )
            ]
        )[["uniprot_ids", "Mean"]]
        return astrocyte_df.rename(
            columns={"uniprot_ids": "uniprot_ids", "Mean": "astrocyte"}
        )
    else:
        cell_type_df = calculate_mean(
            brain_rna_seq_data[
                brain_rna_seq_data.filter(like=cell_type).columns.append(
                    pd.Index(["uniprot_ids"])
                )


]
        )[["uniprot_ids", "Mean"]]
        return cell_type_df.rename(
            columns={"uniprot_ids": "uniprot_ids", "Mean": cell_type}
        )
    # make this into a dictionary like astrocyte : mature and then iterate over it


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
        area_under_line_of_perfect_equality = simps(
            cumulative_fraction_total, cumulative_fraction_total
        )
        area_under_lorenz_curve = simps(
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
    

def cell_type_enrichment(
    brain_rna_seq_data,
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
    astrocytes = mean_cell_type(brain_rna_seq_data, "astrocyte")
    endothelial = mean_cell_type(brain_rna_seq_data, "endothelial")
    microglia = mean_cell_type(brain_rna_seq_data, "microglia")
    oligodendrocytes = mean_cell_type(brain_rna_seq_data, "oligodendrocyte")
    neurons = mean_cell_type(brain_rna_seq_data, "neuron")

    astrocytes = astrocytes[~astrocytes["uniprot_ids"].isna()]
    endothelial = endothelial[~endothelial["uniprot_ids"].isna()]
    microglia = microglia[~microglia["uniprot_ids"].isna()]
    oligodendrocytes = oligodendrocytes[~oligodendrocytes["uniprot_ids"].isna()]
    neurons = neurons[~neurons["uniprot_ids"].isna()]

    all_cell_types = pd.merge(
        pd.merge(
            pd.merge(
                pd.merge(astrocytes, endothelial, on="uniprot_ids"),
                microglia,
                on="uniprot_ids",
            ),
            oligodendrocytes,
            on="uniprot_ids",
        ),
        neurons,
        on="uniprot_ids",
    )
    all_cell_types.set_index("uniprot_ids", inplace=True)

    cell_type_uniprot_ids = []

    if specificity_metric == "enrichment":
        other_cell_types = all_cell_types.drop(cell_type, axis=1)
        cell_type_targets = all_cell_types[[cell_type]]
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
        enrichment_values = all_cell_types.apply(
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
                cell_type_df = all_cell_types.loc[cell_type_uniprots]
                for index, row in cell_type_df.iterrows():
                    max_column = row.idxmax()
                    if max_column == cell_type:
                        cell_type_uniprot_ids.append(index)

        if cell_type == "all":
            if specificity_metric == "tau" or "gini" or "hg" or "tsi":
                df = enrichment_values[
                    enrichment_values < specificity_cutoff
                ].index.tolist()
            else:
                for index, row in enrichment_values.iterrows():
                    cell_type_uniprot_ids = enrichment_values[
                        (enrichment_values[cell_type]) > (specificity_cutoff)
                    ].index.tolist()

    return cell_type_uniprot_ids


def analyze_fractionation(
    tidy_dataframe,
    high_fractions,
    low_fractions,
    sample_health="all",
    mean_median_individual="individual",
):
    """
    Parameters
    ----------
    'tidy_dataframe' : pandas.DataFrame
        DataFrame with one column for each assay, one row for each sample, linearized NPX as the vlaues, and the following indices:
            - 'SampleID'
            - 'Health'
            - 'Sample'
            - 'CSF_sample'
    'high_fractions' : list of strings
        Fractions that should be higher than the fractions in the list of low fractions.
    'low_fractions' : list of strings
        Fractions that should be lower than the fractions in the list of high fractions.
    'sample_health' : {'all', 'ad', 'mci', 'mci_spectrum'}
        Health of the sample requested. Options:
            - 'healthy': only samples from healthy individuals
            - 'all': all different health groups
            - 'ad': samples from individuals diagnosed with Alzheimer's Disease (AD)
            - 'mci': samples from individuals diagnosed with mild cognitive imapirment that has not yet progressed to AD
            - 'mci_spectrum': samples from individuals diagnosed with mild cognitive impairment and samples from individuals that have been diagnosed with AD
    'mean_median_individual' : {'mean', 'median', 'individual', 'individual_median', 'individual_mean'}
        How the groups of samples should be analyzed. Options:
            - 'mean': the means of each fraction should be compared against each other
            - 'median': the medians of each fraction should be compared against each other
            - 'individual': the fractions of each sample should be compared against each other with no aggregation/grouping
            - 'individual_median': for each fraction, the median value of all samples will be compared
            - 'individual_mean': for each fraction, the mean value of all samples will be compared
        Default value: 'individual'
    """

    non_ppa_data = tidy_dataframe[
        tidy_dataframe.index.get_level_values("Sample").str.contains("SEC")
    ]
    if sample_health == "healthy":
        requested_health_data = non_ppa_data[
            non_ppa_data.index.get_level_values("Health").str.contains("Healthy")
        ]
    if sample_health == "all":
        requested_health_data = non_ppa_data
    if sample_health == "ad":
        requested_health_data = non_ppa_data[
            non_ppa_data.index.get_level_values("Health").str.contains("AD")
        ]
    if sample_health == "mci":
        requested_health_data = non_ppa_data[
            non_ppa_data.index.get_level_values("Health").str.contains("MCI")
        ]
    if sample_health == "mci_spectrum":
        requested_health_data = non_ppa_data[
            (non_ppa_data.index.get_level_values("Health").str.contains("AD"))
            | (non_ppa_data.index.get_level_values("Health").str.contains("MCI"))
        ]
    high_fractions_dataframes = {}
    for fraction in high_fractions:
        high_fractions_dataframes[fraction] = requested_health_data[
            requested_health_data.index.get_level_values("Sample").str.contains(
                fraction
            )
        ]
    high_fractions_df = pd.concat(high_fractions_dataframes.values(), axis=0)
    low_fractions_dataframes = {}
    for fraction in low_fractions:
        low_fractions_dataframes[fraction] = requested_health_data[
            requested_health_data.index.get_level_values("Sample").str.contains(
                fraction
            )
        ]
    low_fractions_df = pd.concat(low_fractions_dataframes.values(), axis=0)

    correct_fractionation = []

    for assay in list(non_ppa_data.columns):
        if mean_median_individual == "median":
            if not high_fractions_df[assay].isna().all():
                high_fractions = high_fractions_df[assay].median()

            if not low_fractions_df[assay].isna().all():
                low_fractions = low_fractions_df[assay].median()
            if high_fractions > low_fractions:
                correct_fractionation.append(assay)
        if mean_median_individual == "mean":
            if not high_fractions_df[assay].isna().all():
                high_fractions = high_fractions_df[assay].mean()
            if not low_fractions_df[assay].isna().all():
                low_fractions = low_fractions_df[assay].mean()
            if high_fractions > low_fractions:
                correct_fractionation.append(assay)
        if mean_median_individual == "individual":
            for sample in list(
                requested_health_data.index.get_level_values("CSF_sample").unique()
            ):
                high_sample_data = high_fractions_df[
                    high_fractions_df.index.get_level_values("CSF_sample") == sample
                ][assay].tolist()
                low_sample_data = low_fractions_df[
                    low_fractions_df.index.get_level_values("CSF_sample") == sample
                ][assay].tolist()
                high_fraction = min(high_sample_data)
                low_fraction = max(low_sample_data)
                if high_fraction > low_fraction:
                    correct_fractionation.append(assay)
        if mean_median_individual == "individual_median":
            high_fractions_values = []
            for fraction in high_fractions:
                if not high_fractions_df[assay].isna().all():
                    high_fractions_values.append(
                        high_fractions_df[
                            high_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].median()
                    )
            low_fractions_values = []
            for fraction in low_fractions:
                if not low_fractions_df[assay].isna().all():
                    low_fractions_values.append(
                        low_fractions_df[
                            low_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].median()
                    )
            if high_fractions_values:
                high_fraction = min(high_fractions_values)
            if low_fractions_values:
                low_fraction = max(low_fractions_values)
            if high_fraction > low_fraction:
                correct_fractionation.append(assay)
        if mean_median_individual == "individual_mean":
            high_fractions_values = []
            for fraction in high_fractions:
                if not high_fractions_df[assay].isna().all():
                    high_fractions_values.append(
                        high_fractions_df[
                            high_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].mean()
                    )
            low_fractions_values = []
            for fraction in low_fractions:
                if not low_fractions_df[assay].isna().all():
                    low_fractions_values.append(
                        low_fractions_df[
                            low_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].mean()
                    )
            if high_fractions_values:
                high_fraction = min(high_fractions_values)
            if low_fractions_values:
                low_fraction = max(low_fractions_values)
            if high_fraction > low_fraction:
                correct_fractionation.append(assay)
    return correct_fractionation


def identify_targets(
    correct_fractionation_uniprot_ids, cell_type_uniprot_ids, localization_uniprot_ids
):
    """
    Identifies targets that meet fractionation, cell-type specificity, and localization criteria.
    Parameters:
    ----------
    'correct_fractionation_uniprot_ids': list
        List of UniProt IDs that meet desired fractionation criteria.
    'cell_type_uniprot_ids': list
        List of UniProt IDs that meet desired cell-type specificity criteria.
    'localization_uniprot_ids': list
        List of UniProt IDs that meet desired localization criteria.
    """
    return (
        set(correct_fractionation_uniprot_ids)
        & set(cell_type_uniprot_ids)
        & set(localization_uniprot_ids)
    )


def get_targets(
    assays_path,
    uniprot_fasta_database,
    region,
    brain_rna_seq_path,
    cell_type,
    specificity_metric,
    specificity_cutoff,
    high_fractions,
    low_fractions,
    sample_health,
    mean_median_individual,
    raw_olink_data_file="none",
    plate_layout_dataframe="none",
    tidy_dataframe="none",
    output_directory="ht_output",
):
    """
    Identifies targets that meet specified fractionation, cell-type specificity, and localization criteria.
    Parameters:
    ----------
    'assays_path': .xlsx file path
        Path to a .xlsx file with a column containing the UniProt IDs of all the proteins in the Olink panel.
    'uniprot_fasta_database': .gz file path
        Path to a .gz file containing UniProt IDs and FASTA sequences.
    'region': {'TMhelix', 'inside', 'outside', 'internal', 'external'}
        Subcellular region requested. Options:
          - 'TMhelix': transmembrane proteins
          - 'inside': at least some of the protein is inside the cell/EV
          - 'outside': at least some of the protein is outside the cell/EV
          - 'internal': the protein is only found inside the cell, no transmembrane or outside domains
          - 'external': the protein is only found outside the cell, no transmembrane or inside domains
    'brain_rna_seq_path': csv file path
        Path to the "homo sapiens.csv" file, downloaded from brainrnaseq.org.
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
    'high_fractions' : list of strings
        Fractions that should be higher than the fractions in the list of low fractions.
    'low_fractions' : list of strings
        Fractions that should be lower than the fractions in the list of high fractions.
    'sample_health' : {'all', 'ad', 'mci', 'mci_spectrum'}
        Health of the sample requested. Options:
            - 'healthy': only samples from healthy individuals
            - 'all': all different health groups
            - 'ad': samples from individuals diagnosed with Alzheimer's Disease (AD)
            - 'mci': samples from individuals diagnosed with mild cognitive imapirment that has not yet progressed to AD
            - 'mci_spectrum': samples from individuals diagnosed with mild cognitive impairment and samples from individuals that have been diagnosed with AD
    'mean_median_individual' : {'mean', 'median', 'individual', 'individual_median', 'individual_mean'}
        How the groups of samples should be analyzed. Options:
            - 'mean': the means of each fraction should be compared against each other
            - 'median': the medians of each fraction should be compared against each other
            - 'individual': the fractions of each sample should be compared against each other with no aggregation/grouping
            - 'individual_median': for each fraction, the median value of all samples will be compared
            - 'individual_mean': for each fraction, the mean value of all samples will be compared
        Default value: 'individual'
    'raw_data_path': path to a .csv file, separated by semicolons
        Path to the file containing the raw Olink data.
    'plate_layout_dataframe': pandas.Dataframe
        Dataframe containing information to map the SampleID of a sample to its description.
    'tidy_dataframe' : pandas.DataFrame
        DataFrame with one column for each assay, one row for each sample, linearized NPX as the vlaues, and the following indices:
            - 'SampleID'
            - 'Health'
            - 'Sample'
            - 'CSF_sample'
    'output_directory': directory path
        Path to a directory in which the localization data will be stored.
    """
    if raw_olink_data_file != "none" and plate_layout_dataframe != "none":
        tidy_dataframe = clean_up_raw_data(
            raw_olink_data_file, plate_layout_dataframe
        )

    # tidy_dataframe.dropna(axis=1, how="all", inplace=True)

    fasta_sequences = parse_gz_file(uniprot_fasta_database)
    fasta_sequences.update(
        {
            "NTproBNP": "HPLGSPGSASDLETSGLQEQRNHLQGKLSELQVEQTSLEPLQESPRPTGVWKSREVATEGIRGHRKMVLYTLRAPR",
            "O43521-2": "MAKQPSDVSSECDREGRQLQPAERPPQLRPGAPTSLQTEPQDRSPAPMSCDKSTQTPSPPCQAFNHYLSAMASMRQAEPADMRPEIWIAQELRRIGDEFNAYYARRVFLNNYQAAEDHPRMVILRLLRYIVRLVWRMH",
            "Q13114-2": "MESSKKMDSPGALQTNPPLKLHTDRSAGTPVFVPEQGGYKEKFVKTVEDKYKCEKCHLVLCSPKQTECGHRFCESCMAALLSSSSPKCTACQESIVKDKVFKDNCCKREILALQIYCRNESRGCAEQLMLGHLLVHLKNDCHFEELPCVRPDCKEKVLRKDLRDHVEKACKYREATCSHCKSQVPMIALQVSLLQNESVEKNKSIQSLHNQICSFEIEIERQKEMLRNNESKILHLQRVIDSQAEKLKELDKEIRPFRQNWEEADSMKSSVESLQNRVTELESVDKSAGQVARNTGLLESQLSRHDQMLSVHDIRLADMDLRFQVLETASYNGVLIWKIRDYKRRKQEAVMGKTLSLYSQPFYTGYFGYKMCARVYLNGDGMGKGTHLSLFFVIMRGEYDALLPWPFKQKVTLMLMDQGSSRRHLGDAFKPDPNSSSFKKPTGEMNIASGCPVFVAQTVLENGTYIKDDTIFIKVIVDTSDLPDP",
            "O75882-2": "MVAAAAATEARLRRRTAATAALAGRSGGPHWDWDVTRAGRPGLGAGLRLPRLLSPPLRPRLLLLLLLLSPPLLLLLLPCEAEAAAAAAAVSGSAAAEAKECDRPCVNGGRCNPGTGQCVCPAGWVGEQCQHCGGRFRLTGSSGFVTDGPGNYKYKTKCTWLIEGQPNRIMRLRFNHFATECSWDHLYVYDGDSIYAPLVAAFSGLIVPERDGNETVPEVVATSGYALLHFFSDAAYNLTGFNITYSFDMCPNNCSGRGECKISNSSDTVECECSENWKGEACDIPHCTDNCGFPHRGICNSSDVRGCSCFSDWQGPGCSVPVPANQSFWTREEYSNLKLPRASHKAVVNGNIMWVVGGYMFNHSDYNMVLAYDLASREWLPLNRSVNNVVVRYGHSLALYKDKIYMYGGKIDSTGNVTNELRVFHIHNESWVLLTPKAKEQYAVVGHSAHIVTLKNGRVVMLVIFGHCPLYGYISNVQEYDLDKNTWSILHTQGALVQGGYGHSSVYDHRTRALYVHGGYKAFSANKYRLADDLYRYDVDTQMWTILKDSRFFRYLHTAVIVSGTMLVFGGNTHNDTSMSHGAKCFSSDFMAYDIACDRWSVLPRPDLHHDVNRFGHSAVLHNSTMYVFGGFNSLLLSDILVFTSEQCDAHRSEAACLAAGPGIRCVWNTGSSQCISWALATDEQEEKLKSECFSKRTLDHDRCDQHTDCYSCTANTNDCHWCNDHCVPRNHSCSEGQISIFRYENCPKDNPMYYCNKKTSCRSCALDQNCQWEPRNQECIALPENICGIGWHLVGNSCLKITTAKENYDNAKLFCRNHNALLASLTTQKKVEFVLKQLRIMQSSQSMSKLTLTPWVGLRKINVSYWCWEDMSPFTNSLLQWMPSEPSDAGFCGILSEPSTRGLKAATCINPLNGSVCERPANHSAKQCRTPCALRTACGDCTSGSSECMWCSNMKQCVDSNAYVASFPFGQCMEWYTMSTCPPENCSGYCTCSHCLEQPGCGWCTDPSNTGKGKCIEGSYKGPVKMPSQAPTGNFYPQPLLNSSMCLEDSRYNWSFIHCPACQCNGHSKCINQSICEKCENLTTGKHCETCISGFYGDPTNGGKCQPCKCNGHASLCNTNTGKCFCTTKGVKGDECQLCEVENRYQGNPLRGTCYYTLLIDYQFTFSLSQEDDRYYTAINFVATPDEQNRDLDMFINASKNFNLNITWAASFSAGTQAGEEMPVVSKTNIKEYKDSFSNEKFDFRNHPNITFFVYVSNFTWPIKIQVQTE",
            "Q8WXW3-4": "MSRKISKESKKVNISSSLESEDISLETTVPTDDISSSEEREGKVRITRQLIERKELLHNIQLLKIELSQKTMMIDNLKVDYLTKIEELEEKLNDALHQKQLLTLRLDNQLAFQQKDASKYQELMKQEMETILLRQKQLEETNLQLREKAGDVRRNLRDFELTEEQYIKLKAFPEDQLSIPEYVSVRFYELVNPLRKEICELQVKKNILAEELSTNKNQLKQLTEELAAMKQILVKMHSKHSENSLLLTKTEPKHVTENQKSKTLNVPKEHEDNIFTPKPTLFTKKEAPEWSKKQKMKT",
        }
    )
    assays = pd.read_excel(assays_path)
    assays["Sequence"] = assays["UniProt ID"].map(
        lambda x: fasta_sequences.get(x, "N/A")
    )
    localization_uniprot_ids = identify_localization(assays, region, output_directory)

    brain_rna_seq = map_hgnc_ids(brain_rna_seq_path)
    cell_type_uniprot_ids = cell_type_enrichment(
        brain_rna_seq, cell_type, specificity_metric, specificity_cutoff
    )

    correct_fractionation_uniprot_ids = analyze_fractionation(
        tidy_dataframe,
        high_fractions,
        low_fractions,
        sample_health,
        mean_median_individual,
    )

    return identify_targets(
        correct_fractionation_uniprot_ids,
        cell_type_uniprot_ids,
        localization_uniprot_ids,
    )