"""
Author: Eric J. Ma

Purpose: A small utility that produces Gibson/CPEC assembly primers to be
assembled into a circular plasmid.

Input: One of two ways to feed sequences in:
- A single FASTA file containing nucleotide sequences in clockwise order of
  assembly. The nucleotide sequence should all be the same strand.
- Pass in a list of parts using the set_sequences() function, ensuring that
  they are:
    1. In the correct order of assembly, and
    2. All BioPython SeqRecord objects.

Output: A set of named 40-mer primers required for the assembly, along with
their sequences, and the predicted PCR product size.

Assumptions:
- The region of annealing is not repeated anywhere, i.e. it is unique.
- The plasmid sequence is predicted to be stable (no homologous recombination
  possible).

Dependencies:
- biopython
- pandas
- networkx
- matplotlib
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from collections import defaultdict

import networkx as nx


class PrimerDesigner(nx.DiGraph):
    """
    Sub-classes nx.DiGraph.
    """

    def __init__(self):
        super(PrimerDesigner, self).__init__()
        # The sequences to assemble.
        self.sequences = []

    def read_sequences(self, filename):
        """
        Reads in the PCR products to assemble, which have been formatted as a
        FASTA file. This method does not construct the graph; the separate
        construct_graph() method must be called explicitly.

        Parameters:
        ===========
        - filename: (str) the name of the FASTA file that contains the DNA
                    parts to be assembled.
        """
        for s in SeqIO.parse(filename, 'fasta'):
            s.seq.alphabet = DNAAlphabet()
            self.sequences.append(s)

    def add_sequences(self, seq_records):
        """
        Takes in a list of BioPython SeqRecord objects, in the correct order
        to be assembled, and appends them to the self.sequences attribute.
        """
        self.check_seqrecords(seq_records)
        self.sequences.extend(seq_records)

    def check_seqrecords(self, seq_records):
        """
        A utility function to check that the seqrecords are valid ones.
        """
        # Make sure that seq_records is a list.
        assert isinstance(seq_records, list), "A list of SeqRecords must be\
            passed in."

        # checks on each SeqRecord
        for s in seq_records:
            # Ensure that it is a SeqRecord
            assert isinstance(s, SeqRecord), "A list of SeqRecords must be\
                passed in"

            # Ensure that it has an id that isn't some variant of 'None'
            assert s.id != ''
            assert s.id is not None
            assert s.id != 'None'

            # Automatically change the alphabet of the SeqRecord's Seq object
            # if it is not DNAAlphabet.
            if not isinstance(s.seq.alphabet, DNAAlphabet):
                s.seq.alphabet = DNAAlphabet()

    def set_sequences(self, seq_records):
        """
        Takes in a list of BioPython SeqRecord objects, in the order that they
        are to be assembled, and sets the self.sequences attribute to that.

        This method exists rather than setting the attribute directly, so as
        to perform some basic checks on the SeqRecords that are passed in.
        """
        self.check_seqrecords(seq_records)
        self.sequences = seq_records

    def construct_graph(self):
        """
        Constructs the graph from the list of SeqRecords. Automatically
        computes the necessary primers as well.

        The graph definition is as such:
        - nodes: SeqRecord object
            - node attributes:
                - fw_cloning_primer
                - re_cloning_primer
        - edges: delineating the junction between SeqRecord objects.
            - edge attributes:
                - fw_sequencing_primer: given two adjacent nodes ([a], [b]),
                  selects a sequence on [a] 100 b.p. upstream of the [a], [b]
                  junction.
                - re_sequencing_primer: given two adjacent nodes ([a], [b]),
                  selects a sequence on [b] 100 b.p. downstream of the [a],
                  [b] junction, and takes the reverse complement.
        """
        for i, s in enumerate(self.sequences):
            self.add_node(s)
            if i > 0:
                prev_seq = self.sequences[i - 1]
                self.add_edge(prev_seq, s)
            if i == len(self.sequences) - 1:
                zeroth_seq = self.sequences[0]
                self.add_edge(s, zeroth_seq)

        self.compute_assembly_primers()
        self.compute_junction_sequencing_primers()
        self.compute_fragment_sequencing_primers()

    def compute_assembly_primers(self):
        """
        Given the sequences present in the graph, design primers that are
        15 n.t. overhang and 25 n.t. annealing.
        """
        for n, d in self.nodes(data=True):
            current = n
            predecessor = self.predecessors(n)[0]
            successor = self.successors(n)[0]

            fw_primer = SeqRecord(predecessor.seq[-15:] + current.seq[0:25])
            re_primer = SeqRecord(current.seq[-25:] +
                                  successor.seq[0:15]).reverse_complement()

            self.node[n]['fw_cloning_primer'] = fw_primer
            self.node[n]['re_cloning_primer'] = re_primer

    def compute_junction_sequencing_primers(self):
        """
        For each junction (i.e. edge) in the graph, design primers that are
        positioned at -100 on the upstream part relative to the junction, and
        +100 on the downstream part relative to the junction.

        The junction sequencing primers are stored as edge attributes.
        """
        for upstream, downstream, d in self.edges(data=True):
            fw_primer = SeqRecord(upstream.seq[-125:-100])
            re_primer = SeqRecord(downstream.seq[100:125]).reverse_complement()

            self.edge[upstream][downstream]['fw_sequencing_primer'] = fw_primer
            self.edge[upstream][downstream]['re_sequencing_primer'] = re_primer

    def compute_fragment_sequencing_primers(self):
        """
        For each node in the graph, design primers for sequencing that node.

        The sequencing primers are designed to be 500 n.t. apart from one
        another, which is what is suitable for Sanger sequencing.
        """

        for n, d in self.nodes(data=True):
            # Identify the upstream and downstream parts.
            upstream = self.predecessors(n)[0]
            downstream = self.successors(n)[0]

            # Initialize a list of sequencing primers.
            sequencing_primers = defaultdict(list)

            # Add in fw sequencing primer from the upstream part.
            sequencing_primers['fw'].append(upstream.seq[-125:-100])
            # Add in fw sequencing primers from the current part.
            for pos in range(400, len(n.seq), 500):
                sequencing_primers['fw'].append(n.seq[pos-25:pos])

            # Add in re sequencing primers from the downstream part.
            sequencing_primers['re'].append(
                downstream.seq[100:125].reverse_complement())
            # Add in re sequencing primers from the current part.
            for pos in range(400, len(n.seq), 500):
                sequencing_primers['re'].append(n.seq.reverse_complement()[
                    pos-25:pos])

            # Assign the sequencing primers to the node metadata
            self.node[n]['fragment_sequencing_primers'] = sequencing_primers

    def get_part_ids(self):
        """
        Returns a list of Part IDs.
        """
        return [n.id for n in self.nodes()]

    def has_part(self, id):
        """
        Checks to make sure that the DNA part exists.

        Parameters:
        ===========
        - id: (str) the string that should match one SeqRecord ID.

        Returns:
        ========
        - has_part: (bool) True = the part exists.
        """

        has_part = id in self.get_part_ids()

        return has_part

    def get_part(self, part_name):
        """
        Returns the SeqRecord object whose id is the part_name.

        Paramters:
        ==========
        - part_name:   (str) the `id` of the SeqRecord

        Returns:
        ========
        - part:        (BioPython SeqRecord) the node in the graph that
                       contains the part of interest
        """
        if self.has_part(part_name):
            for n, d in self.nodes(data=True):
                if n.id == part_name:
                    return n
        else:
            raise ValueError('Part {0} not present.'.format(part_name))

    def get_fragment_sequencing_primers(self, part_name):
        """
        Returns the fragment sequencing primers dictionary.
        """
        part = self.get_part(part_name)

        return self.node[part]['fragment_sequencing_primers']

    def pcr_protocol(self):
        """
        Returns a list of dictionaries. The PCR protocol involves the
        following:

        - template: self.node.id
        - fw_primer: self.node[n]['fw_cloning_primer']
        - re_primer: self.node[n]['re_cloning_primer']
        - product_length: len(n.seq) + 30
        - phusion_extension_time: in minutes

        Design note: This format is really flexible, can be converted into a
        pandas dataframe later on. I will not opt to provide a
        save_pcr_protocol function later.
        """

        pcr_protocol = list()
        for n, d in self.nodes(data=True):
            primers = dict()
            primers['template'] = n.id
            primers['fw_cloning_primer'] = d['fw_cloning_primer']
            primers['fw_cloning_primer_name'] = '{0}_fw_cloning_primer'.format(
                n.id)
            primers['re_cloning_primer'] = d['re_cloning_primer']
            primers['re_cloning_primer_name'] = '{0}_re_cloning_primer'.format(
                n.id)
            primers['product_length'] = len(n) + 30
            primers['phusion_extension_time'] = (len(n) + 30) / 1000 * 0.5
            pcr_protocol.append(primers)

        return pcr_protocol

    def get_all_assembly_primers(self):
        """
        Gets all of the assembly primers from each of the nodes.
        """

        primers = []
        for n, d in self.nodes(data=True):
            primers.append(d['fw_cloning_primer'])
            primers.append(d['re_cloning_primer'])

        return primers

    def assembled_plasmid(self):
        """
        Returns a SeqRecord object containing the sequence of the final
        assembled plasmid.
        """

        pass
