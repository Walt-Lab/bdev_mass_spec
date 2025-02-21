import numpy as np
import pandas as pd
import scipy


def calculate_tau_score(array: np.ndarray) -> float:
    """
    Computes the Tau specificity score for gene expression specificity.

    The Tau score measures the degree of selective expression of a gene across
    multiple cell types, with higher values indicating greater specificity.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    float
        The Tau specificity score, ranging from 0 (ubiquitous
        expression) to 1 (high specificity).

    Notes
    -----
    - Normalizes expression by dividing by the maximum value.
    - Calculates the sum of (1 - normalized expression) across all cell
      types.
    - Divides the sum by the number of cell types minus one.
    """
    row_x = array / max(array)
    return np.sum(1 - row_x) / ((len(row_x)) - 1)


def calculate_tsi(array: np.ndarray) -> float:
    """
    Computes the Tissue Specificity Index (TSI) for a gene.

    The TSI measures how much a gene's expression is concentrated in a
    single cell type compared to others, with higher values indicating
    greater specificity.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    float
        The TSI value, where values closer to 1 indicate stronger
        specificity.

    Notes
    -----
    - The formula is: `TSI = max(expression) / sum(expression)`.
    - TSI values near 1 indicate that most expression is concentrated in
      a single cell type.
    - Lower values indicate broader expression across multiple cell
      types.
    """
    return max(array) / sum(array)


def calculate_gini(array: np.ndarray) -> float:
    """
    Computes the Gini coefficient for gene expression specificity.

    The Gini coefficient measures inequality in gene expression
    distribution, with higher values indicating more tissue-specific
    expression.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    float
        The Gini coefficient, where values closer to 1 indicate higher
        specificity.

    Notes
    -----
    - The Gini coefficient is derived from economics and adapted to gene
      expression.
    - A value of 0 represents equal expression across all cell types.
    - A value closer to 1 means that most expression is concentrated in
      a few cell types.
    - Uses the Lorenz curve and calculates areas under the curve for the
      computation.
    """
    sorted_types = np.sort(array)
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


def calculate_hg(array: np.ndarray) -> float:
    """
    Computes Shannon entropy (H(g)) for gene expression specificity.

    Shannon entropy quantifies the dispersion of gene expression across
    cell types, with lower values indicating more specific expression.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    float
        The entropy value, where lower values indicate higher
        specificity.

    Notes
    -----
    - The formula is: `H(g) = -Σ (p_i * log2(p_i))`, where p_i is the
      normalized expression.
    - Entropy values close to 0 indicate that expression is concentrated
      in a single cell type.
    - Higher entropy values indicate broader, less specific expression.
    """
    row_sum = np.sum(array)
    p_sub_i = array / row_sum
    return -1 * (sum(p_sub_i * np.log2(p_sub_i)))


def calculate_spm(array: np.ndarray) -> pd.Series:
    """
    Computes the Specificity Measure (SPM) for gene expression.

    The SPM method squares expression values to emphasize dominant cell
    types and normalizes them to derive specificity.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    pd.Series
        A pandas Series containing the normalized specificity values for
        each cell type.

    Notes
    -----
    - The formula is: `SPM_i = (expression_i)^2 / Σ (expression^2)`.
    - Higher values indicate stronger specificity in a given cell type.
    - Returns a Series with the same index as the input.
    """
    squared_array = array**2
    sum_squared_array = np.sum(squared_array)
    return pd.Series(squared_array / sum_squared_array)


def calculate_zscore(array: np.ndarray) -> pd.Series:
    """
    Computes the Z-score for gene expression across cell types.

    The Z-score standardizes gene expression values relative to the
    dataset's mean and standard deviation.

    Parameters
    ----------
    array : np.ndarray
        A 1D numpy array containing expression levels of a gene across
        multiple cell types.

    Returns
    -------
    pd.Series
        A pandas Series containing the Z-score for each cell type.

    Notes
    -----
    - The formula is: `Z = (X - μ) / σ`, where:
      - X is the expression value
      - μ is the mean expression
      - σ is the standard deviation
    - Higher absolute Z-scores indicate stronger deviation from the
      mean.
    """
    mean_array = np.mean(array)
    std_array = np.std(array)
    return pd.Series((array - mean_array) / std_array)


SPECIFICITY_FUNCTIONS = {
    "tau": calculate_tau_score,
    "tsi": calculate_tsi,
    "gini": calculate_gini,
    "hg": calculate_hg,
    "spm": lambda row: pd.Series(calculate_spm(row), index=row.index),
    "zscore": lambda row: pd.Series(calculate_zscore(row), index=row.index),
}


def calculate_enrichment(row, specificity_metric: str):
    """
    Computes gene expression specificity using a chosen specificity metric.

    Parameters
    ----------
    row : pandas.Series
        A Series containing expression data for a single gene across
        multiple cell types.
    specificity_metric : {'tau', 'tsi', 'gini', 'hg', 'spm', 'zscore'}
        The specificity metric used for calculation. Options:
        - 'tau': Tau specificity score
        - 'tsi': Tissue Specificity Index
        - 'gini': Gini coefficient
        - 'hg': Shannon entropy of gene expression distribution
        - 'spm': Specificity Measure
        - 'zscore': Z-score normalization

    Returns
    -------
    Union[pd.DataFrame, pd.Series, float]
        - If the specificity metric is 'spm' or 'zscore', returns a
          pandas.Series with computed values.
        - Otherwise, returns a float representing the specificity score.

    Raises
    ------
    ValueError
        If an invalid specificity metric is provided.

    Notes
    -----
    - The function applies different specificity calculations depending
      on the chosen metric.
    - The 'spm' and 'zscore' metrics return per-gene values as
      pandas.Series.

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
    row_array = np.array(row)

    if specificity_metric in SPECIFICITY_FUNCTIONS:
        return SPECIFICITY_FUNCTIONS[specificity_metric](row_array)

    raise ValueError(f"Invalid specificity_metric: {specificity_metric}")
