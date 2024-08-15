import pandas as pd

def analyze_fractionation(
    tidy_dataframe,
    high_fractions,
    low_fractions,
    sample_health="all",
    mean_median_individual="individual_median",
) -> list:
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
    'mean_median' : {'mean', 'median'}
        How the groups of samples should be analyzed. Options:
            - 'mean': the mean of all high_fractions will be compared against the mean of all low_fractions
            - 'median': the median of all high_fractions will be compared against the median of all low_fractions
            - 'individual_mean': the mean of each high_fraction will be compared against the mean of each low_fraction
            - 'individual_median': the median of each high_fraction will be compared against the median of each low_fraction
        Default value: 'individual_median'
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
        if mean_median_individual == "individual_median":
            high_fractions_values = []
            for fraction in high_fractions:
                if not high_fractions_df[assay].isna().all():
                    high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].median())
            low_fractions_values = []
            for fraction in low_fractions: 
                if not low_fractions_df[assay].isna().all():
                    low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].median())
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
                    high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
            low_fractions_values = []
            for fraction in low_fractions: 
                if not low_fractions_df[assay].isna().all():
                    low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
            if high_fractions_values:
                high_fraction = min(high_fractions_values)
            if low_fractions_values:
                low_fraction = max(low_fractions_values)
            if high_fraction > low_fraction: 
                correct_fractionation.append(assay)
    return correct_fractionation