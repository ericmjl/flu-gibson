from FluGibson.nucleotide_protein_converter import NucleotideProteinConverter
from Levenshtein import distance
from Bio import SeqIO
import os
import pytest

package_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.join(package_directory, 'test_data'))

np = NucleotideProteinConverter()
np.read_sequences('vic_np_mut2.fasta', 'vic_np_mut1_aa.fasta')
np.convert()

# print(distance(str(np.src_nt.seq), str(np.des_nt.seq)))
print(np.src_nt.seq.alphabet)


def test_convert():
    # There are 4 amino acid changes, but because codons are chosen randomly,
    # based on experimental tests, there could be anywhere between 8 and 12
    # changes, inclusive.

    d = distance(str(np.src_nt.seq), str(np.des_nt.seq))
    assert d <= 12
    assert d >= 8


def test_set_sequences():
    src = SeqIO.read('vic_np_mut2.fasta', 'fasta')
    des = SeqIO.read('vic_np_mut1_aa.fasta', 'fasta')

    # The following should raise an error.
    with pytest.raises(AssertionError):
        np.set_sequences(des, src)
