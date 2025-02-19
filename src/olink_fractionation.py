import pandas as pd

def filter_fractions(data: pd.DataFrame, fractions: list) -> pd.DataFrame:
    """
    Filters a DataFrame to include only the rows that match the specified fractions.

    Parameters
    ----------
    data : pandas.DataFrame
        The input DataFrame containing sample fraction data.
    fractions : list of str
        A list of fraction names to filter.

    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only the selected fractions.

    Notes
    -----
    - This function is used by analyze_fractionation() but can also be reused elsewhere.
    - Returns an empty DataFrame if no matching fractions are found.
    """
    return pd.concat([
        data[data.index.get_level_values("Sample").str.contains(frac)]
        for frac in fractions
    ], axis=0) if fractions else pd.DataFrame()


# def analyze_fractionation(
#     tidy_dataframe: pd.DataFrame,
#     high_fractions: list,
#     low_fractions: list,
#     sample_health: str = "all",
#     mean_median_individual: str = "individual_median",
# ) -> list:
#     """
#     Identifies proteins that show higher expression in high fractions than low fractions.

#     Parameters
#     ----------
#     tidy_dataframe : pandas.DataFrame
#         A DataFrame containing one column per assay, one row per sample, and linearized NPX values.
#         Must include the following index levels:
#             - 'SampleID'
#             - 'Health'
#             - 'Sample'
#             - 'CSF_sample'
#     high_fractions : list of str
#         List of fractions expected to have higher expression levels.
#     low_fractions : list of str
#         List of fractions expected to have lower expression levels.
#     sample_health : {'all', 'healthy', 'ad', 'mci', 'mci_spectrum'}, optional
#         Filter for sample health status. Options:
#         - "'healthy'" : Only samples from healthy individuals.
#         - "'all'" (default) : Includes all available health groups.
#         - "'ad'" : Includes samples from individuals diagnosed with Alzheimer's Disease (AD).
#         - "'mci'" : Includes samples from individuals diagnosed with mild cognitive impairment (MCI).
#         - "'mci_spectrum'" : Includes both MCI and AD samples.
#     mean_median_individual : {'mean', 'median', 'individual_mean', 'individual_median', 'individual_max'}, optional
#         Defines how expression levels should be compared:
#         - "'mean'" : Compare mean of all high fractions vs. mean of all low fractions.
#         - "'median'" : Compare median of all high fractions vs. median of all low fractions.
#         - "'individual_mean'" : Compare mean per fraction, taking the minimum high vs. maximum low.
#         - "'individual_median'" (default) : Compare median per fraction, taking the minimum high vs. maximum low.
#         - "'individual_max'" : Compare maximum expression value per fraction.

#     Returns
#     -------
#     list
#         A list of UniProt assay IDs that satisfy the fractionation criteria.

#     Notes
#     -----
#     - The function filters samples based on "sample_health", selects "high_fractions" and "low_fractions",
#       and applies different statistical comparisons based on "mean_median_individual".
#     - If no valid expression values exist in a fraction, that fraction is skipped in calculations.
#     """

#     non_ppa_data = tidy_dataframe[tidy_dataframe.index.get_level_values("Sample").str.contains("SEC")].copy()

#     health_filters = {
#         "healthy": "Healthy",
#         "ad": "AD",
#         "mci": "MCI",
#         "mci_spectrum": ["AD", "MCI"],
#         "all": None
#     }
#     selected_health = health_filters.get(sample_health, None)

#     if selected_health:
#         requested_health_data = non_ppa_data[non_ppa_data.index.get_level_values("Health").isin(
#             selected_health if isinstance(selected_health, list) else [selected_health]
#         )]
#     else:
#         requested_health_data = non_ppa_data

#     high_fractions_df = filter_fractions(requested_health_data, high_fractions)
#     low_fractions_df = filter_fractions(requested_health_data, low_fractions)

    correct_fractionation = []

    # for assay in non_ppa_data.columns:
    #     if high_fractions_df[assay].notna().any() and low_fractions_df[assay].notna().any():
            
    #         if mean_median_individual == "median":
    #             if high_fractions_df[assay].median() > low_fractions_df[assay].median():
    #                 correct_fractionation.append(assay)

    #         elif mean_median_individual == "mean":
    #             if high_fractions_df[assay].mean() > low_fractions_df[assay].mean():
    #                 correct_fractionation.append(assay)

    #         elif mean_median_individual in ["individual_median", "individual_mean", "individual_max"]:
    #             stat_func = high_fractions_df[assay].median if "median" in mean_median_individual else high_fractions_df[assay].mean
    #             high_values = [stat_func() for frac in high_fractions if high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(frac)][assay].notna().any()]
    #             low_values = [stat_func() for frac in low_fractions if low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(frac)][assay].notna().any()]

    #             if high_values and low_values:
    #                 high_stat = max(high_values) if mean_median_individual == "individual_max" else min(high_values)
    #                 low_stat = max(low_values)
    #                 if high_stat > low_stat:
    #                     correct_fractionation.append(assay)

    # return correct_fractionation

# def analyze_fractionation(
#     tidy_dataframe,
#     high_fractions,
#     low_fractions,
#     sample_health="all",
#     mean_median_individual="individual_median",
# ) -> list:
#     """
#     Parameters
#     ----------
#     'tidy_dataframe' : pandas.DataFrame
#         DataFrame with one column for each assay, one row for each sample, linearized NPX as the vlaues, and the following indices:
#             - 'SampleID'
#             - 'Health'
#             - 'Sample'
#             - 'CSF_sample'
#     'high_fractions' : list of strings
#         Fractions that should be higher than the fractions in the list of low fractions.
#     'low_fractions' : list of strings
#         Fractions that should be lower than the fractions in the list of high fractions.
#     'sample_health' : {'all', 'ad', 'mci', 'mci_spectrum'}
#         Health of the sample requested. Options:
#             - 'healthy': only samples from healthy individuals
#             - 'all': all different health groups
#             - 'ad': samples from individuals diagnosed with Alzheimer's Disease (AD)
#             - 'mci': samples from individuals diagnosed with mild cognitive imapirment that has not yet progressed to AD
#             - 'mci_spectrum': samples from individuals diagnosed with mild cognitive impairment and samples from individuals that have been diagnosed with AD
#     'mean_median_individual' : {'mean', 'median'}
#         How the groups of samples should be analyzed. Options:
#             - 'mean': the mean of all high_fractions will be compared against the mean of all low_fractions
#             - 'median': the median of all high_fractions will be compared against the median of all low_fractions
#             - 'individual_mean': the mean of each high_fraction will be compared against the mean of each low_fraction
#             - 'individual_median': the median of each high_fraction will be compared against the median of each low_fraction
#         Default value: 'individual_median'
#     """

#     non_ppa_data = tidy_dataframe[
#         tidy_dataframe.index.get_level_values("Sample").str.contains("SEC")
#     ]
#     if sample_health == "healthy":
#         requested_health_data = non_ppa_data[
#             non_ppa_data.index.get_level_values("Health").str.contains("Healthy")
#         ]
#     if sample_health == "all":
#         requested_health_data = non_ppa_data
#     if sample_health == "ad":
#         requested_health_data = non_ppa_data[
#             non_ppa_data.index.get_level_values("Health").str.contains("AD")
#         ]
#     if sample_health == "mci":
#         requested_health_data = non_ppa_data[
#             non_ppa_data.index.get_level_values("Health").str.contains("MCI")
#         ]
#     if sample_health == "mci_spectrum":
#         requested_health_data = non_ppa_data[
#             (non_ppa_data.index.get_level_values("Health").str.contains("AD"))
#             | (non_ppa_data.index.get_level_values("Health").str.contains("MCI"))
#         ]

#     high_fractions_dataframes = {}
#     for fraction in high_fractions:
#         high_fractions_dataframes[fraction] = requested_health_data[
#             requested_health_data.index.get_level_values("Sample").str.contains(
#                 fraction
#             )
#         ]
#     high_fractions_df = pd.concat(high_fractions_dataframes.values(), axis=0)
#     low_fractions_dataframes = {}
#     for fraction in low_fractions:
#         low_fractions_dataframes[fraction] = requested_health_data[
#             requested_health_data.index.get_level_values("Sample").str.contains(
#                 fraction
#             )
#         ]
#     low_fractions_df = pd.concat(low_fractions_dataframes.values(), axis=0)

#     correct_fractionation = []

#     for assay in list(non_ppa_data.columns):
#         if mean_median_individual == "median":
#             if not high_fractions_df[assay].isna().all():
#                 high_fractions = high_fractions_df[assay].median()

#             if not low_fractions_df[assay].isna().all():
#                 low_fractions = low_fractions_df[assay].median()
#             if high_fractions > low_fractions:
#                 correct_fractionation.append(assay)
#         if mean_median_individual == "mean":
#             if not high_fractions_df[assay].isna().all():
#                 high_fractions = high_fractions_df[assay].mean()
#             if not low_fractions_df[assay].isna().all():
#                 low_fractions = low_fractions_df[assay].mean()
#             if high_fractions > low_fractions:
#                 correct_fractionation.append(assay)
#         if mean_median_individual == "individual_median":
#             high_fractions_values = []
#             for fraction in high_fractions:
#                 if not high_fractions_df[assay].isna().all():
#                     high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].median())
#             low_fractions_values = []
#             for fraction in low_fractions: 
#                 if not low_fractions_df[assay].isna().all():
#                     low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].median())
#             if high_fractions_values:
#                 high_fraction = min(high_fractions_values)
#             if low_fractions_values:
#                 low_fraction = max(low_fractions_values)
#             if high_fraction > low_fraction: 
#                 correct_fractionation.append(assay)
#         if mean_median_individual == "individual_mean":
#             high_fractions_values = []
#             for fraction in high_fractions:
#                 if not high_fractions_df[assay].isna().all():
#                     high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
#             low_fractions_values = []
#             for fraction in low_fractions: 
#                 if not low_fractions_df[assay].isna().all():
#                     low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
#             if high_fractions_values:
#                 high_fraction = min(high_fractions_values)
#             if low_fractions_values:
#                 low_fraction = max(low_fractions_values)
#             if high_fraction > low_fraction: 
#                 correct_fractionation.append(assay)
#         if mean_median_individual == "individual_max":
#             high_fractions_values = []
#             for fraction in high_fractions:
#                 if not high_fractions_df[assay].isna().all():
#                     high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
#             low_fractions_values = []
#             for fraction in low_fractions: 
#                 if not low_fractions_df[assay].isna().all():
#                     low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
#             if high_fractions_values:
#                 high_fraction = max(high_fractions_values)
#             if low_fractions_values:
#                 low_fraction = max(low_fractions_values)
#             if high_fraction > low_fraction: 
#                 correct_fractionation.append(assay)
#     return correct_fractionation

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
    'mean_median_individual' : {'mean', 'median'}
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

    health_filters = {
        "healthy": "Healthy",
        "ad": "AD",
        "mci": "MCI",
        "mci_spectrum": ["AD", "MCI"],
        "all": None
    }
    selected_health = health_filters.get(sample_health, None)

    if selected_health:
        requested_health_data = non_ppa_data[non_ppa_data.index.get_level_values("Health").isin(
            selected_health if isinstance(selected_health, list) else [selected_health]
        )]
    else:
        requested_health_data = non_ppa_data

    high_fractions_df = filter_fractions(requested_health_data, high_fractions)
    low_fractions_df = filter_fractions(requested_health_data, low_fractions)

    correct_fractionation = []

    for assay in list(non_ppa_data.columns):
        if high_fractions_df[assay].notna().any() and low_fractions_df[assay].notna().any():

            if mean_median_individual == "median":
                high_fractions = high_fractions_df[assay].median()
                low_fractions = low_fractions_df[assay].median()
                if high_fractions > low_fractions:
                    correct_fractionation.append(assay)
            if mean_median_individual == "mean":
                high_fractions = high_fractions_df[assay].mean()
                low_fractions = low_fractions_df[assay].mean()
                if high_fractions > low_fractions:
                    correct_fractionation.append(assay)
            if mean_median_individual == "individual_median":
                high_fractions_values = []
                for fraction in high_fractions:
                    high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].median())
                low_fractions_values = []
                for fraction in low_fractions: 
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
                    high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
                low_fractions_values = []
                for fraction in low_fractions: 
                    low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
                if high_fractions_values:
                    high_fraction = min(high_fractions_values)
                if low_fractions_values:
                    low_fraction = max(low_fractions_values)
                if high_fraction > low_fraction: 
                    correct_fractionation.append(assay)
            if mean_median_individual == "individual_max":
                high_fractions_values = []
                for fraction in high_fractions:
                    high_fractions_values.append(high_fractions_df[high_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
                low_fractions_values = []
                for fraction in low_fractions: 
                    low_fractions_values.append(low_fractions_df[low_fractions_df.index.get_level_values("Sample").str.contains(fraction)][assay].mean())
                if high_fractions_values:
                    high_fraction = max(high_fractions_values)
                if low_fractions_values:
                    low_fraction = max(low_fractions_values)
                if high_fraction > low_fraction: 
                    correct_fractionation.append(assay)
    return correct_fractionation