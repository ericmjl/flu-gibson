"""
Author: Eric J. Ma, MIT
Date: 11 August 2015
"""

from Bio import SeqIO
from itertools import combinations
import networkx as nx

class NucleotideConverter(object):

    """
    A class that takes in two nucleotide FASTA files and returns the primers
    needed to convert the first one to the other.
    """

    def __init__(self, source, destination):
        super(NucleotideConverter, self).__init__()
        # Read in the source and destination FASTA files.
        # They should only contain one sequence.
        self.source = SeqIO.read(source)
        self.destination = SeqIO.read(destination)

        # Check to make sure that source and destination are of the same 
        # length.
        assert len(self.source.seq) == len(self.destination.seq),
            "Source and destination sequences must be of the same length!"

        # Initialize the positions to mutate.
        self.positions = None

    def diff_positions(self):
        """
        Returns a list of positions where source sequence and destination
        sequence differ.
        """
        positions = []
        for i, l1 in enumerate(self.source.seq):
            l2 = self.destination.seq[i]
            if l1 != l2:
                positions.append(i)

        self.positions = positions

    def distance_graph(self):
        """
        Computes the distance between each set of .
        """
        # Initialize a distance graph that 
        cg = nx.DiGraph()
        for p1, p2 in combinations(self.positions, 2):
            dist = abs(p2 - p1)
            if dist > 200:
                cg.add_edge(p1, p2, dist=dist)




