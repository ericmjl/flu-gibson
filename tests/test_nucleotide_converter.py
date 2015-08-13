from FluGibson.nucleotide_converter import NucleotideConverter
from collections import defaultdict
from Levenshtein import distance
import os

# Get the directory of this test file.
package_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(package_directory, 'test_data'))

nc = NucleotideConverter()
nc.read_sequences('vic_np_mut1.fasta', 'vic_np_mut2.fasta')
nc.set_backbone('pCI')
nc.compute_diff_codon_positions()
nc.compute_distance_graph()
nc.compute_assembly_steps()
nc.compute_intermediate_sequences()
nc.compute_pcr_fragments()
nc.compute_mutagenesis_primers()


def test_compute_distance_graph():
    assert len(nc.distance_graph.edges()) == 5


def test_compute_assembly_steps():
    true_protocol = defaultdict(set)
    true_protocol[1] = set([26, 362, 223])
    true_protocol[2] = set([224])
    assert nc.protocol == true_protocol


def test_compute_diff_codon_positions():
    assert nc.codon_positions == set([224, 26, 362, 223])


def test_compute_intermediate_sequences():
    for k, intermediate in nc.intermediates.items():
        if k == 1:
            assert distance(str(intermediate.seq), str(nc.src.seq)) == 6
        if k == 2:
            assert distance(str(intermediate.seq), str(nc.src.seq)) == 8


def test_compute_pcr_fragments():
    frag_len = dict()
    frag_len[1] = set([591, 417, 4490])
    frag_len[2] = set([5498])

    for step, frags in nc.fragments.items():
        lengths = set([len(frag) for frag in frags])
        assert lengths == frag_len[step]
