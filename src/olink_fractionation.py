import pandas as pd


def filter_fractions(data: pd.DataFrame, fractions: list[str]) -> pd.DataFrame:
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
    - This function is used by analyze_fractionation() but can also be
      reused elsewhere.
    - Returns an empty DataFrame if no matching fractions are found.
    """
    if fractions:
        selected_frac_dfs = [
            data[data.index.get_level_values("Sample").str.contains(frac)]
            for frac in fractions
        ]
        return pd.concat(selected_frac_dfs, axis=0)
    else:
        return pd.DataFrame()


def analyze_fractionation(
    tidy_dataframe: pd.DataFrame,
    high_fractions: list,
    low_fractions: list,
    sample_health: str="all",
    mean_median_individual="individual_median",
) -> list:
    """
    Parameters
    ----------
    'tidy_dataframe' : pandas.DataFrame
        DataFrame with one column for each assay, one row for each
        sample, linearized NPX as the vlaues, and the following indices:
            - 'SampleID'
            - 'Health'
            - 'Sample'
            - 'CSF_sample'
    'high_fractions' : list of strings
        Fractions that should be higher than the fractions in the list
        of low fractions.
    'low_fractions' : list of strings
        Fractions that should be lower than the fractions in the list
        of high fractions.
    'sample_health' : {'all', 'ad', 'mci', 'mci_spectrum'}
        Health of the sample requested. Options:
            - 'healthy': only samples from healthy individuals
            - 'all': all different health groups
            - 'ad': samples from individuals diagnosed with Alzheimer's
              Disease (AD)
            - 'mci': samples from individuals diagnosed with mild
              cognitive imapirment that has not yet progressed to AD
            - 'mci_spectrum': samples from individuals diagnosed with
              mild cognitive impairment and samples from individuals
              that have been diagnosed with AD
            - "mci_spectrum"
    'mean_median_individual' : {'mean', 'median'}
        How the groups of samples should be analyzed. Options:
            - 'mean': the mean of all high_fractions will be compared
              against the mean of all low_fractions
            - 'median': the median of all high_fractions will be
              compared against the median of all low_fractions
            - 'individual_mean': the mean of each high_fraction will be
              compared against the mean of each low_fraction
            - 'individual_median': the median of each high_fraction will
              be compared against the median of each low_fraction
        Default value: 'individual_median'
    """

    non_ppa_data = tidy_dataframe[
        tidy_dataframe.index.get_level_values("Sample").str.contains("SEC")
    ]

    health_filters = {
        "healthy": ["Healthy"],
        "ad": ["AD"],
        "mci": ["MCI"],
        "mci_spectrum": ["AD", "MCI"],
        "all": None,
    }
    selected_health = health_filters.get(sample_health, None)

    if selected_health:
        requested_health_data = non_ppa_data[
            non_ppa_data.index.get_level_values("Health").isin(
                selected_health
                if isinstance(selected_health, list)
                else [selected_health]
            )
        ]
    else:
        requested_health_data = non_ppa_data

    high_fractions_df = filter_fractions(requested_health_data, high_fractions)
    low_fractions_df = filter_fractions(requested_health_data, low_fractions)

    correct_fractionation = []
    aggfunc = "median" if mean_median_individual.endswith("median") else "mean"
    high_comparator = max if mean_median_individual.endswith("max") else min
    for assay in list(non_ppa_data.columns):
        if (
            high_fractions_df[assay].notna().any()
            and low_fractions_df[assay].notna().any()
        ):
            if mean_median_individual in ["median", "mean"]:
                high_m = high_fractions_df[assay].agg(aggfunc)
                low_m = low_fractions_df[assay].agg(aggfunc)
                if high_m > low_m:
                    correct_fractionation.append(assay)
            elif mean_median_individual.startswith("individual"):
                high_fractions_values = []
                for fraction in high_fractions:
                    high_fractions_values.append(
                        high_fractions_df[
                            high_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].agg(aggfunc)
                    )
                low_fractions_values = []
                for fraction in low_fractions:
                    low_fractions_values.append(
                        low_fractions_df[
                            low_fractions_df.index.get_level_values(
                                "Sample"
                            ).str.contains(fraction)
                        ][assay].agg(aggfunc)
                    )
                if high_fractions_values and low_fractions_values:
                    high_fraction = high_comparator(high_fractions_values)
                    low_fraction = max(low_fractions_values)
                    if high_fraction > low_fraction:
                        correct_fractionation.append(assay)
    return correct_fractionation
