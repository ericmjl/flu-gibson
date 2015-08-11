"""
Author: Eric J. Ma, MIT
Date: 11 August 2015
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import defaultdict
from itertools import combinations
import networkx as nx


class NucleotideConverter(object):

    """
    A class that takes in two nucleotide FASTA files and returns the primers
    needed to convert the first one to the other.
    """

    def __init__(self):
        super(NucleotideConverter, self).__init__()

        # The source sequence.
        self.source = None

        # The destination sequence.
        self.destination = None

        # A list of codon positions where source and destination differ
        self.codon_positions = None

        # A NetworkX graph that stores the distance between each codon that
        # needs to be mutated.
        self.distance_graph = None

    def read_sequences(self, source, destination):
        """
        Reads in the source and destination sequences.
        """
        src = SeqIO.read(source, 'fasta', alphabet=generic_dna)
        des = SeqIO.read(destination, 'fasta', alphabet=generic_dna)

        assert len(src.seq) == len(des.seq)

        self.source = src
        self.destination = des

    def compute_diff_codon_positions(self):
        """
        Sets the self.codon_positions attribute to a set of codon_positions
        where source sequence and destination sequence differ.
        """

        # Initialize an empty set of codon positions
        codon_positions = set()
        for i, _ in enumerate(self.source.seq[::3]):
            codon1 = self.source.seq[(i * 3):(i * 3 + 3)]
            codon2 = self.destination.seq[(i * 3):(i * 3 + 3)]
            if str(codon1) != str(codon2):
                codon_positions.add(i)

        self.codon_positions = codon_positions

    def compute_distance_graph(self):
        """
        Computes the distance between each pair of codons. They should be at
        least 50 codons (i.e. 150 n.t.) apart for them to be joined together.
        """

        assert isinstance(self.codon_positions, set)

        # Initialize a networkx directed graph.
        G = nx.DiGraph()

        # Iterate over every pair of positions. If they are greater than 50
        # codons apart, draw an edge between them.
        for p1, p2 in combinations(self.codon_positions, 2):
            dist = abs(p2 - p1)
            if dist > 50:
                G.add_edge(min(p1, p2), max(p1, p2), dist=dist)

        self.distance_graph = G

    def compute_assembly_steps(self):
        """
        Computes the mutation assembly steps based on the distance graph.

        The algorithm is as such:

        - while number of nodes > 0:
            - find smallest node
            - greedily add downstream until no more destination nodes
        """
        # Check first to make sure that the distance_graph has been made.
        assert isinstance(self.distance_graph, nx.DiGraph)

        # Copy the distance graph.
        G = self.distance_graph.copy()

        def _assembly_step(G, current):
            """
            Helper function that computes the current iteration assembly step.
            """
            step = set([current])

            def next_nearest_node(G, n):
                if len(G.successors(n)) == 0:
                    return None
                else:
                    return min(G.successors(n))

            while next_nearest_node(G, current):
                step.add(next_nearest_node(G, current))
                current = next_nearest_node(G, current)

            return step

        # Initialize the assembly step dictionary
        protocol = defaultdict(set)

        # Iterate over nodes
        i = 0
        while len(G.nodes()) > 0:
            i += 1
            # Find node with smallest number.
            start = min(G.nodes())
            protocol[i] = _assembly_step(G, start)

            G.remove_nodes_from(protocol[i])

        self.protocol = protocol
