import pandas as pd

import requests
import time

from raw_data_preprocessing import clean_up_raw_data
from olink_fractionation import analyze_fractionation
from brainrnaseq_specificity import map_hgnc_ids, cell_type_enrichment, create_enrichment_dataframe
from deeptmhmm_localization import parse_gz_file, identify_localization

def uniprot_to_protein_name(uniprot_ids):
    from_db = "UniProtKB_AC-ID"
    to_db = "UniProtKB"

    mapped_ids = uniprot_id_mapping(from_db, to_db, list(uniprot_ids))
    result_dict = mapped_ids if mapped_ids else {}

    return set(result_dict.values())

def protein_name_to_uniprot(uniprot_ids):
    from_db = "UniProtKB" 
    to_db = "UniProtKB_AC-ID"

    mapped_ids = uniprot_id_mapping(from_db, to_db, list(uniprot_ids))
    result_dict = mapped_ids if mapped_ids else {}

    return set(result_dict.values())

def uniprot_id_mapping(from_db, to_db, ids):
    url = "https://rest.uniprot.org/idmapping/run"
    data = {
        "from": from_db,
        "to": to_db,
        "ids": " ".join(ids)
    }
    response = requests.post(url, data=data)
    response.raise_for_status()
    job_id = response.json()["jobId"]
    
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    
    max_attempts = 5
    attempts = 0
    while attempts < max_attempts:
        try:
            status_response = requests.get(status_url)
            status_response.raise_for_status()
            status_data = status_response.json()
            
            if "jobStatus" in status_data:
                if status_data["jobStatus"] == "RUNNING":
                    time.sleep(5)
                elif status_data["jobStatus"] == "FINISHED":
                    break
                else:
                    print(f"Unexpected status: {status_data['jobStatus']}")
                    return None
            elif "results" in status_data:
                return process_results(status_data)
            else:
                print("Unexpected response format")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Error checking job status: {e}")
            attempts += 1
            time.sleep(5)
    
    if attempts == max_attempts:
        print("Max attempts reached. Could not get job status.")
        return None
    
    try:
        results_response = requests.get(results_url)
        results_response.raise_for_status()
        return process_results(results_response.json())
    except requests.exceptions.RequestException as e:
        print(f"Error getting results: {e}")
        return None

def process_results(data):
    accession_to_id = {}
    if "results" in data:
        for result in data["results"]:
            to_data = result.get("to", {})
            uniprot_kb_id = to_data.get("uniProtkbId")
            if uniprot_kb_id:
                primary_accession = to_data.get("primaryAccession")
                if primary_accession:
                    accession_to_id[primary_accession] = uniprot_kb_id
                secondary_accessions = to_data.get("secondaryAccessions", [])
                for sec_acc in secondary_accessions:
                    accession_to_id[sec_acc] = uniprot_kb_id
    return accession_to_id


def identify_targets(
    assay_list_path,
    uniprot_fasta_database,
    brain_rna_seq_raw_path,
    region,
    cell_type,
    specificity_metric,
    specificity_cutoff,
    high_fractions,
    low_fractions,
    sample_health,
    mean_median_individual="median",
    raw_olink_data_file="none",
    plate_layout_dataframe="none",
    tidy_dataframe="none",
    output_directory="ht_output",
) -> set:
    """
    Identifies targets that meet specified fractionation, cell-type specificity, and localization criteria. Returns a set of UniProt IDs.
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
    'brain_rna_seq_raw_path': csv file path
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
    assays = pd.read_excel(assay_list_path)
    assays["Sequence"] = assays["UniProt ID"].map(
        lambda x: fasta_sequences.get(x, "N/A")
    )
    localization_uniprot_ids = identify_localization(assays, region, output_directory)

    brain_rna_seq = map_hgnc_ids(brain_rna_seq_raw_path)
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

    return (set(correct_fractionation_uniprot_ids) & set(cell_type_uniprot_ids) & set(localization_uniprot_ids))