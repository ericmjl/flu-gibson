from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import defaultdict

import matplotlib.pyplot as plt
import sys
import pandas as pd
import networkx as nx


class PrimerDesigner(object):
    """
    A small utility that produces Gibson/CPEC assembly primers to be assembled
    into a circular plasmid.

    Input: One of two ways to feed sequences in:
    - A single FASTA file containing nucleotide sequences in clockwise order
      of assembly. The nucleotide sequence should all be the same strand.
    - Pass in a list of parts using the set_sequences() function, ensuring
      that they are:
        1. In the correct order of assembly, and
        2. All BioPython SeqRecord objects.

    Output: A set of named 40-mer primers required for the assembly, along with
    their sequences, and the predicted PCR product size.

    Assumptions:
    - The region of annealing is not repeated anywhere, i.e. it is unique.
    - The plasmid sequence is predicted to be stable (no homologous
      recombination possible).

    Dependencies:
    - biopython
    - pandas
    - networkx
    - matplotlib
    """

    def __init__(self):
        super(PrimerDesigner, self).__init__()
        # The FASTA file to assemble.
        self.filename = None

        # The sequences to assemble.
        self.sequences = None

        # A NetworkX graph that stores the assembly graph.
        self.graph = nx.DiGraph()

    def read_sequences(self, filename):
        """
        Reads in the PCR products to assemble, which have been formatted as a
        FASTA file.
        """
        self.filename = filename
        self.sequences = [s for s in SeqIO.parse(filename, 'fasta')]

    def set_sequences(self, seq_records):
        """
        Takes in a list of BioPython SeqRecord objects, in the order that they
        are to be assembled, and sets the self.sequences attribute to that.

        This method exists rather than setting the attribute directly, so as
        to perform some basic checks on the SeqRecords that are passed in.
        """
        # Make sure that seq_records is a list.
        assert isinstance(seq_records, list)

        # checks on each SeqRecord
        for s in seq_records:
            # Ensure that it is a SeqRecord
            assert isinstance(s, SeqRecord)

            # Ensure that it has an id that isn't some variant of 'None'
            assert s.id != ''
            assert s.id is not None
            assert s.id != 'None'

        self.sequences = seq_records

    def construct_graph(self):
        """
        Constructs the graph from the list of SeqRecords.

        The graph definition is as such:
        - nodes: SeqRecord object
        - edges: delineating the order of SeqRecord objects.
        """
        for i, s in enumerate(self.sequences):
            self.graph.add_node(s)
            if i > 0:
                prev_seq = self.sequences[i - 1]
                self.graph.add_edge(prev_seq, s)
            if i == len(self.sequences) - 1:
                zeroth_seq = self.sequences[0]
                self.graph.add_edge(s, zeroth_seq)

    def list_part_ids(self):
        """
        Returns a list of Part IDs.
        """
        return [n.id for n in self.graph.nodes()]

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

        part = None
        for n, d in self.graph.nodes(data=True):
            if n.id == part_name:
                part = n

        if part is None:
            raise ValueError('Part not present.')
        else:
            return part

    def _fragment_sequencing_primers(self, part_name):
        """
        Returns the fragment sequencing primers dictionary.
        """
        part = self.get_part(part_name)

        return self.graph.node[part]['fragment_sequencing_primers']

    def get_fragment_sequencing_primers(self, part_name):
        """
        Returns a pandas DataFrame of the sequencing primers that need to be
        ordered.
        """
        primer_dict = self._fragment_sequencing_primers(part_name)

        primer_list = []
        for direction, primers in primer_dict.items():
            for i, primer in enumerate(primers):
                primer_entry = dict()
                primer_entry['part'] = part_name
                primer_entry['direction'] = direction
                primer_entry['sequence'] = str(primer)
                primer_entry['idx'] = i

                primer_list.append(primer_entry)

        return pd.DataFrame(primer_list)

    def print_fragment_sequencing_primers(self, part_name):
        """
        Uses the _fragment_sequencing_primers() function to get a list of
        sequencing primers for a given part. It then pretty prints those
        primers.

        Paramters:
        ==========
        - part_name:    (str) the name of the DNA part of interest.
        """

        primer_dict = self._fragment_sequencing_primers(part_name)

        for direction, primers in primer_dict.items():
            print(direction)
            for primer in primers:
                print(primer)

    def design_assembly_primers(self):
        """
        Given the sequences present in the graph, design primers that are
        15 n.t. overhang and 25 n.t. annealing.
        """
        for n, d in self.graph.nodes(data=True):
            current = n
            predecessor = self.graph.predecessors(n)[0]
            successor = self.graph.successors(n)[0]

            fw_primer = SeqRecord(predecessor.seq[-15:] + current.seq[0:25])
            re_primer = SeqRecord(current.seq[-25:] +
                                  successor.seq[0:15]).reverse_complement()

            self.graph.node[n]['fw_sequence'] = fw_primer
            self.graph.node[n]['re_sequence'] = re_primer
            self.graph.node[n]['fw_primer_name'] = self.filename.split('.')[0] + '_' + n.id + '_fw'
            self.graph.node[n]['re_primer_name'] = self.filename.split('.')[0] + '_' + n.id + '_re'

    def design_junction_sequencing_primers(self):
        """
        For each junction (i.e. edge) in the graph, design primers that are
        positioned at -100 on the upstream part relative to the junction, and
        +100 on the downstream part relative to the junction.
        """

        for upstream, downstream, d in self.graph.edges(data=True):
            fw_primer = SeqRecord(upstream.seq[-125:-100])
            re_primer = SeqRecord(downstream.seq[100:125]).reverse_complement()

            self.graph.node[upstream]['fw_sequencing_primer'] = fw_primer
            self.graph.node[downstream]['re_sequencing_primer'] = re_primer

    def design_fragment_sequencing_primers(self):
        """
        For each node in the graph, design primers for sequencing that node.
        """

        for n, d in self.graph.nodes(data=True):
            # Identify the upstream and downstream parts.
            upstream = self.graph.predecessors(n)[0]
            downstream = self.graph.successors(n)[0]

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
            self.graph.node[n]['fragment_sequencing_primers'] = \
                sequencing_primers

    def compute_pcr_protocol(self):
        """
        Returns a list of dictionaries.

        Design note: This format is really flexible, can be converted into a
        pandas dataframe later on.
        """
        assert self.filename is not None
        assert self.filename != ''

        pcr_protocol = list()
        for n, d in self.graph.nodes(data=True):
            primers = dict()

            f_id = self.filename.split('.')[0] + '_' + n.id + '_fw'
            r_id = self.filename.split('.')[0] + '_' + n.id + '_re'

            primers['fw_sequence'] = d['fw_sequence'].seq
            primers['re_sequence'] = d['re_sequence'].seq
            primers['fw_primer_name'] = d['fw_primer_name']
            primers['re_primer_name'] = d['re_primer_name']
            primers['fw_len'] = len(d['fw_sequence'])
            primers['re_len'] = len(d['re_sequence'])
            primers['template'] = n.id
            primers['pcr_length'] = len(n) + 30
            primers['phusion_extension_minutes'] = (len(n) + 30) / 1000 * 0.5
            primers['fw_sequencing_primer'] = d['fw_sequencing_primer'].seq
            primers['re_sequencing_primer'] = d['re_sequencing_primer'].seq
            primers['pcr_template'] = n.id
            pcr_protocol.append(primers)

        self.pcr_protocol = pcr_protocol

    def save_pcr_protocol(self):
        """
        Saves the PCR protocol as a CSV file.
        """

        pd.DataFrame(self.pcr_protocol).to_csv(
            "{filename}.csv".format(filename=self.filename.split('.')[0]))
