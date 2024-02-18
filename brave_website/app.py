from flask import Flask, render_template, request

import sys
sys.path.append("C:\\Users\\Wyss User\\Documents\\EVs\\OLINK\\brave_website\\olink_analysis_final.ipynb")

import olink_analysis_final.ipynb
app = Flask(__name__)

# Define routes
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    uniprot_id = request.form['uniprot_id']
    # Call a function or method from your Python module
    results = olink_analysis_final.graph_medians(uniprot_id)
    return render_template('search_results.html', results=results)

# Define other routes as needed

if __name__ == '__main__':
    app.run(debug=True)
