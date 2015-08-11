from FluGibson.nucleotide_converter import NucleotideConverter
from collections import defaultdict
import os

# Get the directory of this test file.
package_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(package_directory, 'test_data'))

nc = NucleotideConverter()
nc.read_sequences('vic_np_mut1.fasta', 'vic_np_mut2.fasta')
nc.compute_diff_codon_positions()
nc.compute_distance_graph()
nc.compute_assembly_steps()

for i, _ in enumerate(nc.source.seq[::3]):
    codon1 = nc.source.seq[(i * 3):(i * 3 + 3)]
    codon2 = nc.destination.seq[(i * 3):(i * 3 + 3)]
    if str(codon1) != str(codon2):
        print(i, codon1, codon2)

print(nc.codon_positions)
print(nc.distance_graph.edges())
print(nc.protocol)


def test_compute_distance_graph():
    assert len(nc.distance_graph.edges()) == 5

def test_compute_assembly_steps():
    true_protocol = defaultdict(set)
    true_protocol[1] = set([26, 362, 223])
    true_protocol[2] = set([224])
    assert nc.protocol == true_protocol

def test_compute_diff_codon_positions():
    assert nc.codon_positions == set([224, 26, 362, 223])