from flask import Flask, request, jsonify
import ht_analysis
from ht_analysis import assay_list_path, uniprot_fasta_database, brain_rna_seq_path

app = Flask(__name__)

@app.route('/filter', methods=['GET'])
def filter():
    # Get arguments from the request
    region = request.args.get('Target Localization')
    cell_type = request.args.get('Cell Type of Interest')
    specificity_metric = request.args.get('Specificity Metric')
    specificity_cutoff = request.args.get('Specificity Cutoff')
    high_fractions = request.args.get('High Concentration Fractions')
    low_fractions = request.args.get('Low Concentration Fractions')
    sample_health = request.args.get('Sample Health')
    mean_median_individual = request.args.get('Aggregation Metric')

    # Call your function with the arguments
    result = ht_analysis.get_targets(
    assays_path = assay_list_path,
    uniprot_fasta_database = uniprot_fasta_database,
    region = region,
    brain_rna_seq_path = brain_rna_seq_path,
    cell_type = cell_type,
    specificity_metric = specificity_metric,
    specificity_cutoff = specificity_cutoff,
    high_fractions = high_fractions,
    low_fractions = low_fractions,
    sample_health = sample_health,
    mean_median_individual = "median",
    raw_olink_data_file="none",
    plate_layout_dataframe="none",
    tidy_dataframe="none",
    output_directory="ht_output",)

    # Return the result as JSON
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True)
