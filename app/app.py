from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from FluGibson.primer_designer import PrimerDesigner
from FluGibson.plasmids import plasmids
from IPython.display import HTML

import pandas as pd

app = Flask(__name__)


@app.route('/')
def main():
    # Get a list of all of the available backbones.
    backbones = sorted([i for i in plasmids.keys()])
    print(backbones)
    return render_template('index_new.html', plasmids=backbones)


def validate_sequence_characters(sequence):
    """
    Validates that the sequence passed into it consists only of A, T, G and Cs.
    """
    sequence = set(sequence.upper())
    allowed = ['A', 'T', 'G', 'C']
    for s in sequence:
        try:
            assert s in allowed
            return True
        except AssertionError:
            return False


def validate_sequence_length(sequence):
    """
    Validates that the sequence passed into it has a minimum length of 100 n.t.
    """
    try:
        assert len(sequence) >= 100
        return True
    except AssertionError:
        return False


# @app.errorhandler(600)
def error(message):
    """
    Displays the error message.
    """
    return render_template('error.html', message=message)


@app.route('/compute_primers', methods=['POST'])
def compute_primers():
    # Get the part names.
    sequence_name = request.form['sequence_name']
    sequence = request.form['sequence'].replace('\n', '')
    backbone = request.form['backbone']

    # Validation checks
    if not validate_sequence_characters(sequence):
        return error('Sequence contains non-ATGC letters')
    if not validate_sequence_length(sequence):
        return error('Sequence length has to be at least 100 n.t.')
    else:
        # Create a SeqRecord object.
        parts = [SeqRecord(Seq(sequence.upper(), alphabet=generic_dna),
                           name=sequence_name, id=sequence_name)]
        parts.extend(plasmids[backbone])

        p = PrimerDesigner()
        p.set_sequences(parts)
        p.construct_graph()

        nodes = p.nodes(data=True)
        edges = p.edges(data=True)
        protocol = p.pcr_protocol()

        return render_template('primers.html',
                               nodes=nodes,
                               edges=edges,
                               protocol=protocol)

if __name__ == '__main__':
    app.run(debug=True, port=5005)
