from FluGibson.nucleotide_protein_converter import NucleotideProteinConverter
from Levenshtein import distance
import os

package_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(package_directory, 'test_data'))

np = NucleotideProteinConverter()
np.read_sequences('vic_np_mut2.fasta', 'vic_np_mut1_aa.fasta')
np.convert()

print(distance(str(np.src_nt.seq), str(np.des_nt.seq)))

def test_convert():
    # There are 4 amino acid changes, but because codons are chosen randomly,
    # there could be anywhere between 10 and 12 changes inclusive.
    
    d = distance(str(np.src_nt.seq), str(np.des_nt.seq))
    assert d <= 12
    assert d > 9