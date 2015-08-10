from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import matplotlib.pyplot as plt
import sys
import pandas as pd
import networkx as nx

"""
A small utility that produces Gibson/CPEC assembly primers to be assembled
into a circular plasmid.

Input: A single FASTA file containing nucleotide sequences in clockwise order
of assembly. The nucleotide sequence should all be the same strand.

Output: A set of named primers required for the assembly, along with their
sequences, and the predicted PCR product size.

Assumptions:
- The region of annealing is not repeated anywhere, i.e. it is unique.
- The plasmid sequence is predicted to be stable (no homologous recombination
  possible).

Dependencies:
- biopython
- pandas
- networkx
"""


class PrimerDesigner(object):

    """docstring for PrimerDesigner"""

    def __init__(self, filename):
        super(PrimerDesigner, self).__init__()
        self.filename = filename
        self.sequences = [s for s in SeqIO.parse(filename, 'fasta')]
        self.graph = nx.DiGraph()

    def construct_graph(self):
        for i, s in enumerate(self.sequences):
            self.graph.add_node(s)
            if i > 0:
                prev_seq = self.sequences[i - 1]
                self.graph.add_edge(prev_seq, s)
            if i == len(self.sequences) - 1:
                zeroth_seq = self.sequences[0]
                self.graph.add_edge(s, zeroth_seq)

    def design_assembly_primers(self):
        '''
        Given the sequences present in the graph, design primers that are
        15 n.t. overhang and 25 n.t. annealing.
        '''
        for n, d in self.graph.nodes(data=True):
            current = n
            predecessor = self.graph.predecessors(n)[0]
            successor = self.graph.successors(n)[0]

            fw_primer = SeqRecord(predecessor.seq[-15:] + current.seq[0:25])
            re_primer = SeqRecord(current.seq[-25:] +
                                  successor.seq[0:15]).reverse_complement()

            self.graph.node[n]['fw_sequence'] = fw_primer
            self.graph.node[n]['re_sequence'] = re_primer
            self.graph.node[n]['fw_primer'] = self.filename.split(
                '.')[0] + '_' + n.id + '_fw'
            self.graph.node[n]['re_primer'] = self.filename.split(
                '.')[0] + '_' + n.id + '_re'

    def design_sequencing_primers(self):
    	"""
    	For each junction (i.e. edge) in the graph, design primers that are 
    	positioned at -100 on the upstream part relative to the junction, and 
    	+100 on the downstream part relative to the junction.
    	"""

    	for upstream, downstream, d in self.graph.edges(data=True):
    		fw_primer = SeqRecord(upstream.seq[-120:-100])
    		re_primer = SeqRecord(downstream.seq[100:120]).reverse_complement()

    		self.graph.node[upstream]['fw_sequencing_primer'] = fw_primer
    		self.graph.node[downstream]['re_sequencing_primer'] = re_primer

    def compute_pcr_protocol(self):
        """
        Returns a list of dictionaries.
        This format is really flexible, can be converted into a pandas
        dataframe later on.
        """
        pcr_protocol = list()
        for n, d in self.graph.nodes(data=True):
            primers = dict()

            f_id = self.filename.split('.')[0] + '_' + n.id + '_fw'
            r_id = self.filename.split('.')[0] + '_' + n.id + '_re'

            primers['fw_sequence'] = d['fw_sequence'].seq
            primers['re_sequence'] = d['re_sequence'].seq
            primers['fw_primer'] = d['fw_primer']
            primers['re_primer'] = d['re_primer']
            primers['fw_len'] = len(d['fw_sequence'])
            primers['re_len'] = len(d['re_sequence'])
            primers['template'] = n.id
            primers['pcr_length'] = len(n) + 30
            primers['phusion_extension_minutes'] = (len(n) + 30) / 1000 * 0.5
            primers['fw_sequencing_primer'] = d['fw_sequencing_primer'].seq
            primers['re_sequencing_primer'] = d['re_sequencing_primer'].seq
            pcr_protocol.append(primers)

        self.pcr_protocol = pcr_protocol

    def save_pcr_protocol(self):
        """
        Saves the PCR protocol as a CSV file.
        """

        pd.DataFrame(self.pcr_protocol).to_csv(
            "{filename}.csv".format(filename=self.filename.split('.')[0]))


# if __name__ == '__main__':

#     filename = sys.argv[1]
#     pdesigner = PrimerDesigner(filename)
#     pdesigner.design_primers()
#     pdesigner.compute_pcr_protocol()
#     pdesigner.save_pcr_protocol()
