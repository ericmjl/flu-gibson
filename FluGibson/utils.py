"""
A series of utility functions that provide useful information.
"""
from random import choice
from Bio.Data import IUPACData

codons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']
    }


def get_codon(amino_acid):
    if len(amino_acid) == 1:
        amino_acid = IUPACData.protein_letters_1to3[amino_acid].upper()
    return choice(codons[amino_acid])
