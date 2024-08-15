import matplotlib.pyplot as plt
import pandas as pd

CSF_SAMPLES = [
    "SEC Fract 6 ",
    "SEC Fract 7",
    "SEC Fract 8",
    "SEC Fract 9",
    "SEC Fract 10",
    "SEC Fract 11",
    "SEC Fract 12",
    "SEC Fract 13",
    "SEC Fract 14",
    "SEC Fract 15",
]


def clean_strings(strings) -> list:
    """
    Gets rid of possible sources of error among strings in a list.
    Parameters
    ----------
    strings : a list of strings.
    """
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


def clean_up_raw_data(raw_data, plate_layout_path) -> pd.DataFrame:
    """
    Normalizes all data points and converts to a tidy dataframe.
    Parameters
    ----------
    raw_data: .parquet file path
        Path to the raw parquet file provided by Olink containing sample data.
    plate_layout_path: .xlsx file path
        Path to an excel file used to map the SampleIDs to any information.
    """
    ht_data = pd.read_parquet(raw_data)
    plate_layout = pd.read_excel(plate_layout_path)
    data = pd.merge(ht_data, plate_layout, how="left", on="SampleID")
    # repeat assays: P32455, Q02750
    replicate_assays = data[["UniProt"]].value_counts()
    replicate_assays = replicate_assays[replicate_assays > 94]
    replicate_assay_list = clean_strings(replicate_assays.index.tolist())
    unique_data = data[~data["UniProt"].isin(replicate_assay_list)]
    unique_data.loc[:, "Linear NPX"] = unique_data["PCNormalizedNPX"].map(lambda x: 2**x)
    tidy_data = unique_data[unique_data["SampleType"] == "SAMPLE"].pivot(
        columns="UniProt",
        index=["SampleID", "Health", "Sample", "CSF_sample"],
        values="Linear NPX",
    )
    return tidy_data


def find_ratio(df, high_fractions, low_fractions) -> int:
    """
    Calculates the EV association score for each column of a tidy dataframe.
    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe containing the data for one protein collected by the Olink HT panel.
    high_fractions: list
        A list of fractions that would be expected to form the EV peak if the protein is associated with EVs.
    low_fractions: list
        A list of the fractions that are expected to be present in lower concentrations if the protein is associated with EVs.
    """
    high_fractions_or = '|'.join(high_fractions)
    peaking_fracts = df[
        (df.index.get_level_values("Sample").str.contains(high_fractions_or))
        & (df.index.get_level_values("Health") == "Healthy")
    ].median()

    low_fractions_or = '|'.join(low_fractions)
    low_fracts = df[
        (df.index.get_level_values("Sample").str.contains(low_fractions_or))
        & (df.index.get_level_values("Health") == "Healthy")
    ].median()
    return peaking_fracts / low_fracts

def ev_association_score_df(tidy_data, high_fractions, low_fractions) -> pd.DataFrame:
    """
    Returns a DataFrame containing the EV association score for each protein in the Olink HT panel.
    Parameters
    ----------
    tidy_data: pandas.DataFrame
        Tidy DataFrame containing the data collected by the Olink HT panel.
    high_fractions: list
        A list of fractions that would be expected to form the EV peak if the protein is associated with EVs.
    low_fractions: list
        A list of the fractions that are expected to be present in lower concentrations if the protein is associated with EVs.
    """
    ht_assay = []
    ht_ratio = []

    for assay in tidy_data.columns:
        df = tidy_data[
            (tidy_data.index.get_level_values("Health") == "Healthy")
            & ~(tidy_data.index.get_level_values("CSF_sample").str.contains("Internal"))
        ][assay]
        ratio = find_ratio(df, high_fractions, low_fractions)
        ht_assay.append(assay)
        ht_ratio.append(ratio)

    return(pd.DataFrame({"ht_assay": ht_assay, "ht_ratio": ht_ratio}))


def plot_protein_fractionation(tidy_data, uniprot_id) -> plt.boxplot:
    """
    Returns a box-and-whisker plot to depict the fractionation pattern for a specified protein using data collected by the Olink HT panel. Red lines represent the median, boxes represent the interquartile range, lines represent the range excluding outliers, and dots represent the outliers.
    Parameters
    ----------
    tidy_data: pandas.DataFrame
        Tidy DataFrame containing the data collected by the Olink HT panel.
    uniprot_id: string
        UniProt ID corresponding to the protein of interest.
    """
    df = tidy_data[uniprot_id]
    df = df[df.index.get_level_values("Health") == "Healthy"]
    df = df.reset_index(level=["SampleID", "Health", "Sample"])
    df["Sample"] = pd.Categorical(df["Sample"], categories=CSF_SAMPLES, ordered=True)
    df_sorted = df.sort_values("Sample")
    grouped_data = [
        group[uniprot_id].values for name, group in df_sorted.groupby("Sample")
    ]
    plt.boxplot(grouped_data, notch=None, vert=None, patch_artist=None, widths=None)
    plt.xlabel("Sample Description")
    plt.ylabel("Delta")
    plt.title(f"Healthy {uniprot_id} Fractionation Pattern, HT Panel")
    plt.xticks(range(1, len(CSF_SAMPLES) + 1), CSF_SAMPLES)
    plt.xticks(rotation=45, ha="right")
    plt.show()

