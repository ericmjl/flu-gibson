"""
Author: Eric J. Ma

Purpose:
This Python module provides a dictionary for retrieving plasmid backbones as a
BioPython SeqRecord object.
"""

import os
from Bio import SeqIO


plasmid_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'plasmid_backbones')

# print(os.listdir(plasmid_dir))

plasmids = {f.strip('.fasta'): SeqIO.read(os.path.join(plasmid_dir, f),
                                         'fasta') for f in
            os.listdir(plasmid_dir)}
print(plasmids)