from FluGibson.primer_designer import PrimerDesigner
import os
from Bio.SeqRecord import SeqRecord
import pytest

# Get the directory of this test file.
package_directory = os.path.dirname(os.path.abspath(__file__))

os.chdir(os.path.join(package_directory, 'test_data'))

p = PrimerDesigner()
p.read_sequences('victoria_np.fasta')
p.construct_graph()


def test_read_sequences():
    assert len(p.sequences) == 2


def test_construct_graph():
    assert len(p.nodes()) == 2
    assert len(p.edges()) == 2


def test_has_part():
    assert p.has_part('Vic_NP')
    assert p.has_part('pCI')


def test_compute_pcr_protocol():
    pass


def test_get_all_assembly_primers():

    primers = p.get_all_assembly_primers()
    assert len(primers) == 4


def test_design_fragment_sequencing_primers():
    for n, d in p.nodes(data=True):
        frag_seq_primers = d['fragment_sequencing_primers']

        # Check that the forward primers are designed correctly.
        for i, primer in enumerate(frag_seq_primers['fw']):
            if i >= 1:
                sequence = p.node[n]['object'].seq
                assert str(primer) in str(sequence)
            else:
                # assert str(primer) in str(p.predecessors(n)[0].seq)
                predecessor_sequence = \
                    p.node[p.predecessors(n)[0]]['object'].seq
                assert str(primer) in predecessor_sequence
        # Check that the reverse primers are designed correctly.
        for i, primer in enumerate(frag_seq_primers['re']):
            if i >= 1:
                sequence = p.node[n]['object'].seq.reverse_complement()
                assert str(primer) in str(sequence)

            else:
                successor_sequence = \
                    p.node[p.successors(n)[0]]['object']\
                    .seq.reverse_complement()
                assert str(primer) in str(successor_sequence)


def test_node_attributes_are_correct():
    """
    Ensures that each node in the graph conforms to the data model specified.
    """
    for n, d in p.nodes(data=True):
        assert isinstance(d['object'], SeqRecord)
        assert isinstance(n, str)
        assert 'fw_cloning_primer' in d.keys()
        assert 're_cloning_primer' in d.keys()

    for u, v, d in p.edges(data=True):
        assert 'fw_sequencing_primer' in d.keys()
        assert 're_sequencing_primer' in d.keys()


def test_get_fragment_sequencing_primers():
    primers = p.get_fragment_sequencing_primers("Vic_NP")
    for _, pr in primers.items():  # _ = direction, pr = primers
        assert len(pr) == 4, print(pr)

    with pytest.raises(ValueError):
        p.get_fragment_sequencing_primers("Hello")
