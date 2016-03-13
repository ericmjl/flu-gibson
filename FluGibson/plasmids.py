"""
Author: Eric J. Ma

Purpose:
This Python module provides a dictionary for retrieving plasmid backbones as a
BioPython SeqRecord object.
"""

import os
from Bio import SeqIO

plasmid_dir = 'plasmid_backbones'
plasmids = {f.strip('fasta'): SeqIO.read(f) for f in os.listidr(plasmid_dir)}
