from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data import CodonTable
from FluGibson import utils

class NucleotideProteinConverter(object):
    """
    A class that performs converts one nucleotide sequence into a destination
    protein sequence, by performing the necessary mutations required.
    """
    def __init__(self):
        super(NucleotideProteinConverter, self).__init__()
        # Source sequence. Should be nucleotides.
        self.src_nt = None

        # Source sequence, in amino acids.
        self.src_aa = None

        # Destination sequence. Should be amino acids.
        self.des_aa = None

        # Destination sequence, in nucleotides.
        self.des_nt = None

    def read_sequences(self, source, destination):
        """
        Reads in the source and destination FASTA files.
        """

        self.src_nt = SeqIO.read(source, 'fasta', alphabet=generic_dna)
        self.des_aa = SeqIO.read(destination, 
                                 'fasta', 
                                 alphabet=generic_protein)

    def set_sequences(self, source, destination):
        """
        As an alternative to reading the sequences, directly set the
        attributes while providing type checking.

        Assume that the source is a nucleotide SeqRecord, and the destination
        is a protein SeqRecord.
        """

        assert isinstance(source, SeqRecord)
        assert isinstance(source.seq.alphabet, generic_dna)
        assert isinstance(destination, SeqRecord)
        assert isinstance(destination.seq.alphabet, generic_protein)

        self.src = source
        self.des = destination
    
    def convert(self):
        """
        Performs the conversion from source nucleotide sequence to destination
        nucleotide sequence.

        Because we are seeking to mutate specific codons, and do not want to
        change any other 3rd codons, we will have to build the new sequence
        iteratively.
        """
        new_sequence = ''
        for codon_pos, _ in enumerate(self.src_nt.seq[::3]):
            # Get the codon at that codon position
            codon = self.src_nt.seq[codon_pos*3:codon_pos*3 + 3]
            # print(codon_pos, codon)
            # Get the translated codon
            codon_tr = codon.translate()

            # Check if translated codon is the same as the amino acid sequence.
            if str(codon_tr) != self.des_aa.seq[codon_pos]:
                print('Two codons not the same!')
                print(codon_tr, self.des_aa.seq[codon_pos])

                # If they aren't the same, add a reverse translated codon. 
                # The get_codon function is called from utils.py
                new_codon = utils.get_codon(self.des_aa.seq[codon_pos])

                new_sequence += str(new_codon)

            else:
                new_sequence += str(codon)

        self.des_nt = SeqRecord(new_sequence, id='{0}_nt'.format(
                self.des_aa.id))

        print(self.des_nt)

    def save_mutated_sequence(self, filename):
        """
        Saves the mutated sequence in a FASTA file format.
        """

        SeqIO.write(self.des_nt, filename, 'fasta')


