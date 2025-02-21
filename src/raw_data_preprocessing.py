import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes

from config import CSF_SAMPLES


def clean_strings(strings: list[str]) -> list[str]:
    """
    Gets rid of possible sources of error among strings in a list.

    Parameters
    ----------
    strings : a list of str
        A list containing strings that may contain unwanted characters.

    Returns
    ----------
    list of str
        A list with cleaned strings, where parenthases, single quotes,
        and commas are removed.

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
    Normalizes protein expression data and converts it to a tidy dataframe.

    This function loads raw expression data and a plate layout file,
    merges them to match sample IDs with metadata, removes duplicate
    assays, and normalizes NPX values.

    Parameters
    ----------
    raw_data: str or Path
        Path to the raw .parquet file provided by Olink containing
        sample data.
    plate_layout_path: str or Path
        Path to an excel (.xlsx) file used to map the SampleIDs to any
        information provided.

    Returns
    ----------
    pandas.DataFrame
        A tidy DataFrame where rows represent samples and columns
        represent proteins. Expression values are normalized and pivoted
        into a structured format.
    """
    ht_data = pd.read_parquet(raw_data)
    plate_layout = pd.read_excel(plate_layout_path)
    data = pd.merge(ht_data, plate_layout, how="left", on="SampleID")
    # repeat assays: P32455, Q02750
    replicate_assays = data[["UniProt"]].value_counts()
    replicate_assays = replicate_assays[replicate_assays > 94]
    replicate_assay_list = clean_strings(replicate_assays.index.tolist())
    unique_data = data[~data["UniProt"].isin(replicate_assay_list)].copy()
    unique_data.loc[:, "Linear NPX"] = unique_data["PCNormalizedNPX"].map(
        lambda x: 2**x
    )
    tidy_data = unique_data[unique_data["SampleType"] == "SAMPLE"].pivot(
        columns="UniProt",
        index=["SampleID", "Health", "Sample", "CSF_sample"],
        values="Linear NPX",
    )
    return tidy_data


def find_ratio(
    df: pd.DataFrame | pd.Series, high_fractions: list[str], low_fractions: list[str]
) -> int:
    """
    Calculates the extracellular vesicle (EV) association score for each protein.

    The EV association score is calculated as the ratio of median
    expression in the high fractions (expected EV peak) to median
    expression in low fractions (no signal from EVs expected).

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe containing the data for one protein collected by the
        Olink HT panel.
    high_fractions: list of str
        A list of fractions that would be expected to form the EV peak
        if the protein is associated with EVs.
    low_fractions: list of str
        A list of the fractions that are expected to be present in lower
        concentrations if the protein is associated with EVs.

    Returns
    ----------
    pandas.Series
        The calculated EV associated score for each protein.
    """
    healthy_ind = df.index.get_level_values("Health") == "Healthy"
    fracts = []
    for fractions in [high_fractions, low_fractions]:
        fractions_or = "|".join(fractions)
        fract_ind = df.index.get_level_values("Sample").str.contains(fractions_or)
        fracts.append(df[healthy_ind & fract_ind].median())
    return fracts[0] / fracts[1]


def calculate_fractionation_scores(proteins, tidy_data, high_fractions, low_fractions):
    return [
        find_ratio(tidy_data[protein], high_fractions, low_fractions)
        for protein in proteins
    ]


def ev_association_score_df(
    tidy_data: pd.DataFrame, high_fractions: list[str], low_fractions: list[str]
) -> pd.DataFrame:
    """
    Calculates the EV association score for each protein in the Olink HT panel.

    Parameters
    ----------
    tidy_data: pandas.DataFrame
        Tidy DataFrame containing the data collected by the Olink HT
        panel, where rows represent samples and columns represent
        proteins.
    high_fractions: list of str
        A list of fractions that would be expected to form the EV peak
        if the protein is associated with EVs.
    low_fractions: list of str
        A list of the fractions that are expected to be present in lower
        concentrations if the protein is associated with EVs.

    Returns
    ----------
    pandas.DataFrame
        A DataFrame with two columns:
        - 'ht_assay': UniProt IDs of proteins.
        - 'ht_ratio': Corresponding EV association scores.
    """

    internal_healthy = tidy_data[
        (tidy_data.index.get_level_values("Health") == "Healthy")
        & ~(tidy_data.index.get_level_values("CSF_sample").str.contains("Internal"))
    ]
    ht_ratio = tidy_data.columns.map(
        lambda assay: find_ratio(internal_healthy[assay], high_fractions, low_fractions)
    )

    return pd.DataFrame({"ht_assay": tidy_data.columns, "ht_ratio": ht_ratio})


def plot_protein_fractionation(tidy_data: pd.DataFrame, uniprot_id: str) -> Axes:
    """
    Box-and-whisker plot showing the fractionation pattern of a specified protein.

    The plot visualizes the distribution of protein expression across
    different SEC fractions.
    - Red lines represent the median.
    - Boxes represent the interquartile range.
    - Lines extend to the full range, excluding outliers.
    - Dots indicate outliers.

    Parameters
    ----------
    tidy_data: pandas.DataFrame
        Tidy DataFrame containing the data collected by the Olink HT
        panel, where rows represent samples and columns represent
        proteins.
    uniprot_id: str
        UniProt ID corresponding to the protein of interest.

    Returns
    ----------
    matplotlib.Axes
        A matplotlib Axes object containing the fractionation plot.
    """
    df = (
        tidy_data.loc[
            tidy_data.index.get_level_values("Health") == "Healthy", uniprot_id
        ]
        .reset_index(level=["SampleID", "Health", "Sample"])
        .copy()
    )
    df["Sample"] = pd.Categorical(df["Sample"], categories=CSF_SAMPLES, ordered=True)
    df_grouped = df.sort_values("Sample").groupby("Sample")
    grouped_data = [group[uniprot_id].values for __, group in df_grouped]

    fig, ax = plt.subplots()
    ax.boxplot(grouped_data, notch=None, vert=None, patch_artist=None, widths=None)
    ax.set_xlabel("Sample Description")
    ax.set_ylabel("Delta")
    ax.set_title(f"Healthy {uniprot_id} Fractionation Pattern, HT Panel")
    ax.set_xticks(range(1, len(CSF_SAMPLES) + 1))
    ax.set_xticklabels(CSF_SAMPLES, rotation=45, ha="right")
    plt.tight_layout()
    return ax
