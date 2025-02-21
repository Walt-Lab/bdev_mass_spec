from typing import Optional
from pathlib import Path

import pandas as pd

from brainrnaseq_specificity import (
    cell_type_enrichment,
    create_enrichment_dataframe,
    map_hgnc_ids,
)
from config import MISSING_FASTA_SEQUENCES, hgnc_ids
from deeptmhmm_localization import identify_localization, parse_gz_file
from olink_fractionation import analyze_fractionation
from raw_data_preprocessing import calculate_fractionation_scores, clean_up_raw_data


def identify_targets(
    assay_list_path: str | Path,
    uniprot_fasta_database: str | Path,
    brain_rna_seq_raw_path: str | Path,
    region: str,
    cell_type: str,
    specificity_metric: str,
    specificity_cutoff: float,
    high_fractions: list[str],
    low_fractions: list[str],
    sample_health: str,
    mean_median_individual: str = "median",
    raw_olink_data_file: Optional[str | Path] = None,
    plate_layout_dataframe: Optional[pd.DataFrame] = None,
    tidy_dataframe: Optional[pd.DataFrame] = None,
    output_directory: str | Path = "ht_output",
) -> set[str]:
    """
    Identifies proteins that meet specified fractionation, cell-type
    specificity, and localization criteria.

    This function integrates **fractionation analysis**,
    **gene expression specificity**, and **subcellular localization**
    to identify proteins that meet all three criteria simultaneously.

    Parameters
    ----------
    assay_list_path : str
        Path to a .xlsx file containing a column of UniProt IDs for
        proteins in the Olink panel.
    uniprot_fasta_database : str
        Path to a .gz file containing UniProt IDs and FASTA sequences.
    brain_rna_seq_raw_path : str
        Path to the "homo sapiens.csv" file, downloaded from
        brainrnaseq.org.
    region : {'TMhelix', 'inside', 'outside', 'internal', 'external'}
        Subcellular localization category to filter proteins. Options:
          - 'TMhelix' : Transmembrane proteins
          - 'inside' : At least some of the protein is inside the
            cell/EV
          - 'outside' : At least some of the protein is outside the
            cell/EV
          - 'internal' : The protein is only found inside the cell, with
            no transmembrane or outside domains
          - 'external' : The protein is only found outside the cell,
            with no transmembrane or inside domains
    cell_type : {'astrocyte', 'endothelial', 'microglia',
    'oligodendrocyte', 'neuron'}
        Cell type for specificity filtering. Options:
         - 'astrocyte' : Mature astrocytes
         - 'endothelial' : Endothelial cells
         - 'microglia' : Microglia cells
         - 'oligodendrocyte' : Oligodendrocytes
         - 'neuron' : Neurons
    specificity_metric : {'tsi', 'zscore', 'spm', 'tau', 'gini', 'hg'}
        Metric for evaluating cell-type specificity. Options:
        - 'tsi' : Tissue Specificity Index
        - 'zscore' : Z-score
        - 'spm' : Specificity measure
        - 'tau' : Tau index
        - 'gini' : Gini coefficient
        - 'hg' : Entropy of a gene's expression distribution
    specificity_cutoff : float
        Minimum value required for a protein to pass the specificity
        threshold.
    high_fractions : list of str
        List of fraction names where expression is expected to be
        higher.
    low_fractions : list of str
        List of fraction names where expression is expected to be lower.
    sample_health : {'all', 'healthy', 'ad', 'mci', 'mci_spectrum'}
        Sample health category for filtering fractionation data.
        Options:
            - 'all' : Includes all available health groups
            - 'healthy' : Only samples from healthy individuals
            - 'ad' : Samples from individuals diagnosed with Alzheimer's
              Disease (AD)
            - 'mci' : Samples from individuals diagnosed with mild
              cognitive impairment (MCI)
            - 'mci_spectrum' : Samples from individuals with either
              MCI or AD
    mean_median_individual : {'mean', 'median', 'individual',
    'individual_median', 'individual_mean'}, optional
        Method for comparing sample fraction values. Options:
            - 'mean' : Uses the mean of each fraction
            - 'median' : Uses the median of each fraction
            - 'individual' : Compares individual samples without
              aggregation
            - 'individual_median' : Compares median values across
              samples per fraction
            - 'individual_mean' : Compares mean values across samples
              per fraction
        Default is ''median''.
    raw_olink_data_file : str, optional
        Path to a .csv file containing raw Olink data
        (semicolon-separated).
    plate_layout_dataframe : pandas.DataFrame, optional
        A DataFrame mapping SampleIDs to their descriptions.
    tidy_dataframe : pandas.DataFrame, optional
        A preprocessed DataFrame with the following structure:
            - One column per assay
            - One row per sample
            - Values: Linearized NPX values
            - Indexed by:
                - 'SampleID'
                - 'Health'
                - 'Sample'
                - 'CSF_sample'
    output_directory : str, optional
        Directory path to save localization data. Default is
        "ht_output".

    Returns
    -------
    set
        A set of UniProt IDs that meet the fractionation, cell-type
        specificity, and localization criteria.

    Notes
    -----
    - If 'tidy_dataframe' is not provided, it is generated from
      'raw_olink_data_file' and 'plate_layout_dataframe'.
    - The function relies on multiple sub-functions from:
        - 'raw_data_preprocessing' (for cleaning raw Olink data)
        - 'olink_fractionation' (for fractionation analysis)
        - 'brainrnaseq_specificity' (for gene expression specificity)
        - 'deeptmhmm_localization' (for protein localization predictions)
    - The function requires a pre-downloaded UniProt FASTA database.
    """

    if tidy_dataframe is None:
        tidy_dataframe = clean_up_raw_data(raw_olink_data_file, plate_layout_dataframe)

    fasta_sequences = parse_gz_file(uniprot_fasta_database)
    fasta_sequences.update(MISSING_FASTA_SEQUENCES)
    assays = pd.read_excel(assay_list_path)
    assays["Sequence"] = assays["UniProt ID"].map(
        lambda x: fasta_sequences.get(x, "N/A")
    )
    localization_uniprot_ids = identify_localization(assays, region, output_directory)

    brain_rna_seq = map_hgnc_ids(hgnc_ids, brain_rna_seq_raw_path)
    expression_df = create_enrichment_dataframe(brain_rna_seq)
    cell_type_uniprot_ids = cell_type_enrichment(
        expression_df, cell_type, specificity_metric, specificity_cutoff
    )

    correct_fractionation_uniprot_ids = analyze_fractionation(
        tidy_dataframe,
        high_fractions,
        low_fractions,
        sample_health,
        mean_median_individual,
    )

    return (
        set(correct_fractionation_uniprot_ids)
        & set(cell_type_uniprot_ids)
        & set(localization_uniprot_ids)
    )


def generate_protein_dataframe(
    low_tau: pd.DataFrame,
    fractionation_ids: list[str],
    localization_ids: list[str],
    localization_label: str,
    tidy_data: pd.DataFrame,
    high_fractions: list[str],
    low_fractions: list[str],
) -> pd.DataFrame:
    """
    Generates a filtered protein dataframe based on tau scores,
    fractionation patterns, and localization.

    Parameters
    ----------
    low_tau: pd.DataFrame
        DataFrame containing proteins with low tau scores.
    fractionation_ids: list
        List of UniProt IDs with the correct fractionation pattern.
    localization_ids: list
        List of UniProt IDs corresponding to a specific localization
        category.
    localization_label: str
        Label indicating the localization category (e.g., "internal",
        "transmembrane").
    tidy_data: DataFrame
        Dataset containing protein expression data.
    high_fraction: list of str
        List of high fraction numbers used for ratio calculation.
    low_fractions: list of str
        List of low fraction numbers used for ratio calculation.

    Returns
    ----------
    DataFrame
        A processed DataFrame containing filtered proteins with
        calculated fractionation scores and localization labels.
    """
    filtered_proteins = low_tau[
        low_tau.index.get_level_values("uniprot_ids").isin(localization_ids)
        & low_tau.index.get_level_values("uniprot_ids").isin(fractionation_ids)
    ]
    scores = calculate_fractionation_scores(
        filtered_proteins.index.get_level_values("uniprot_ids"),
        tidy_data,
        high_fractions,
        low_fractions,
    )
    filtered_proteins = filtered_proteins.reset_index()
    filtered_proteins["ev_association_score"] = scores
    filtered_proteins["localization"] = localization_label
    return filtered_proteins


def extract_sample_npx(
    protein_list: list[str], tidy_data_sec: pd.DataFrame
) -> pd.DataFrame:
    """
    Processes a list of proteins to extract median NPX values for each sample.

    Parameters
    ----------
    protein_list: list
        List of protein UniProt IDs to process.
    tidy_data_sec: DataFrame
        Filtered dataset containing SEC data from healthy samples.

    Returns
    ----------
    DataFrame
        A DataFrame containing columns for protein ID, sample, and
        median NPX values.
    """

    proteins, samples, median_npx = [], [], []

    for column in protein_list:
        df = tidy_data_sec[column]
        for sample in df.index.get_level_values("Sample").unique():
            proteins.append(column)
            samples.append(sample)
            sample_df = df[df.index.get_level_values("Sample") == sample]
            median_npx.append(sample_df.median())

    return pd.DataFrame(
        {"protein": proteins, "sample": samples, "median_npx": median_npx}
    )
