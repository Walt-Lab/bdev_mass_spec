import biolib
import gzip
import os
import pathlib
import pandas as pd

def parse_gz_file(file_path: str) -> dict:
    """
    Parses a .gz FASTA file to extract UniProt IDs and their corresponding protein sequences.

    Parameters
    ----------
    file_path : str
        Path to the .gz file containing UniProt IDs and FASTA sequences.

    Returns
    -------
    dict
        A dictionary where keys are UniProt IDs and values are the corresponding FASTA sequences.

    Notes
    -----
    - The function assumes that the FASTA header follows the UniProt format: `>sp|UniProtID|Protein Name`.
    - If a line in the FASTA file does not follow the expected format, it is skipped with a warning.
    - The function processes multi-line FASTA sequences by concatenating them.
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


def tmhmm_localization(assays: pd.DataFrame, output_directory: str):
    """
    Runs DeepTMHMM to predict transmembrane regions and protein localization.

    This function creates a FASTA file from the input DataFrame, runs DeepTMHMM locally, 
    and saves the results in the specified output directory.

    Parameters
    ----------
    assays : pandas.DataFrame
        A DataFrame containing at least two columns:
        - 'UniProt ID': The UniProt identifier of the protein.
        - 'Sequence': The amino acid sequence of the protein.
    output_directory : str
        Path to the directory where DeepTMHMM results will be saved.

    Notes
    -----
    - This function requires the DeepTMHMM model (`DTU/DeepTMHMM:1.0.24`) to be available.
    - The FASTA file is temporarily created as `query.fasta`.
    - The function runs DeepTMHMM locally using `biolib`.
    """
    deeptmhmm = biolib.load('DTU/DeepTMHMM:1.0.24')
    assay_list = []
    with open("query.fasta", "w") as fasta_file:
        for _, row in assays.iterrows():
            fasta_line = f">{row['UniProt ID']}\n{row['Sequence']}\n"
            fasta_file.write(fasta_line)
            assay_list.append(fasta_line)
    deeptmhmm_job = deeptmhmm.cli(args="--fasta query.fasta", machine="local")
    deeptmhmm_job.save_files(output_directory)

def get_regional_uniprots(localization_df: pd.DataFrame, region_name: str) -> set:
    """
    Extracts a set of UniProt IDs associated with a given localization region.

    Parameters
    ----------
    localization_df : pandas.DataFrame
        A DataFrame containing protein localization data from DeepTMHMM.
    region_name : str
        The localization region of interest (e.g., 'TMhelix', 'inside', 'outside').

    Returns
    -------
    set
        A set of UniProt IDs corresponding to proteins localized in the specified region.
    """
    return set(localization_df[localization_df["region_location"] == region_name]["uniprot_id"])


def identify_localization(
    assays: pd.DataFrame, region: str, output_directory: str = "ht_output"
) -> set:
    """
    Identifies proteins localized in specific subcellular regions using DeepTMHMM predictions.

    This function reads DeepTMHMM output and returns a set of UniProt IDs for proteins 
    found in a specified localization category.

    Parameters
    ----------
    assays : pandas.DataFrame
        A DataFrame containing at least two columns:
        - 'UniProt ID': The UniProt identifier of the protein.
        - 'Sequence': The amino acid sequence of the protein.
    region : {'TMhelix', 'inside', 'outside', 'internal', 'external'}
        The subcellular region of interest. Options:
        - 'TMhelix' : Transmembrane proteins
        - 'inside' : At least part of the protein is inside the cell/EV
        - 'outside' : At least part of the protein is outside the cell/EV
        - 'internal' : The protein is only found inside the cell, with no transmembrane or outside domains
        - 'external' : The protein is only found outside the cell, with no transmembrane or inside domains
    output_directory : str, optional
        Path to the directory where DeepTMHMM output files are stored. Default is "ht_output".

    Returns
    -------
    set
        A set of UniProt IDs corresponding to proteins localized in the specified region.

    Notes
    -----
    - If DeepTMHMM output does not exist in the specified directory, it will be generated.
    - The function reads localization data from the "TMRs.gff3" output file.
    - For 'internal' and 'external' localizations, proteins are filtered by excluding transmembrane and opposing-region proteins.
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

    if region == "internal":
        tm_uniprots = get_regional_uniprots(localization_df, "TMhelix")
        outside_uniprots = get_regional_uniprots(localization_df, "outside")
        inside_uniprots = get_regional_uniprots(localization_df, "inside")
        return inside_uniprots - tm_uniprots - outside_uniprots

    elif region == "external":
        tm_uniprots = get_regional_uniprots(localization_df, "TMhelix")
        outside_uniprots = get_regional_uniprots(localization_df, "outside")
        inside_uniprots = get_regional_uniprots(localization_df, "inside")
        return outside_uniprots - inside_uniprots - tm_uniprots
    
    return get_regional_uniprots(localization_df, region)

def get_localization_data(
        assays: pd.DataFrame, 
        fasta_sequences: dict, 
        localization_types: list, 
        output_dir: str
    ) -> dict:
    """
    Identifies protein localizations based on UniProt sequences and predefined localization types.
    
    Parameters
    ----------
    assays: DataFrame
        DataFrame containing assay information, including UniProt IDs.
    fasta_sequences: dict
        Dictionary mapping UniProt IDs to their respective protein sequences.
    localization_types: {"TMhelix", "internal", "external", "inside", "outside"}
        The subcellular region of interest. Options:
        - 'TMhelix' : Transmembrane proteins
        - 'inside' : At least part of the protein is inside the cell/EV
        - 'outside' : At least part of the protein is outside the cell/EV
        - 'internal' : The protein is only found inside the cell, with no transmembrane or outside domains
        - 'external' : The protein is only found outside the cell, with no transmembrane or inside domains
    output_dir: str
        Path to directory to store output files.
    
    Returns
    ----------
    dict
        A dictionary where keys are localization types and values are lists of UniProt IDs corresponding to each localization.
    """
    assays["Sequence"] = assays["UniProt ID"].map(lambda x: fasta_sequences.get(x, "N/A"))
    return {loc: identify_localization(assays, loc, output_dir) for loc in localization_types}