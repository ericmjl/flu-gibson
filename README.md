# FluGibson

![Travis Status](https://travis-ci.org/ericmjl/flu-gibson.svg) [![Code Climate](https://codeclimate.com/github/ericmjl/flu-gibson/badges/gpa.svg)](https://codeclimate.com/github/ericmjl/flu-gibson)

A tool for designing primers to clone influenza polymerase segments from viral cDNA.

# Installation

The installation requires the following packages:

1. `networkx`
2. `biopython`
3. `pandas`
4. `matplotlib`
5. `flask`
6. `ipython`

I have provided an `environment.yml` file for `conda` users.

*From Github (recommended):*

1. Download this repository as a Zip file.
2. Unzip the file.
3. In your terminal, navigate to the FluGibson directory.
4. Run command: `python setup.py install`

*From PyPI:*

1. (if applicable) Switch to your proper Python environment.
2. Run command: `pip install FluGibson`

*Using Conda:*

1. (if applicable) Switch to your proper Python environment.
2. Run command: `conda install FluGibson`

# Usage

## Web App

The web app in the `app` directory provides a graphical user interface to operate FluGibson. To run it:

    $ cd app
    $ python app.py

## Scripted

One way to use FluGibson is to use the provided script in the `/examples` directory. Copy the script to your working directory. 

Create the FASTA formatted files containing the DNA parts that you want to stitch together. For example, you would use the following FASTA definition to stitch the following 3 parts together:

    >PART_1
    >CATCTATCTCTCTACTGCGAGGCTATTCGACTGGCCGTTACTCGCCGGTACGTAGCTCGGTCTCGATCATCAGTACGTCTACGTGTCGTCGTACTTACACGGTCGCTCGGACTGACGTACGTCTACGTCGTCTGACTGA
    
    >PART_2
    >CTACTGTCTGCTGATGGTACGTACGTGAGTACGCGCAGCACAGACACTACTTACTCTCGCGCGAGAGCTATCTACGACTACGTACTCGTCGTACGAGCTGACTGATCGACGTAGCTTGACGTACGTATCACGTACGTATCG
    
    >PART_3
    >CAGCTTCGGCGCGATTACTCTACGAGCACGACGCAGCTGTCGCTGTCTGGTCTACGCTAGCGCTACGACTATCGATCAGCGTCGTACTGACGTGACGCGCATCGACGTTCGGACGTCGTCGTCGTACGACGTCTACGATGC

The parts will be joined in the order `PART_1-->PART_2-->PART_3`.

To produce the CSV file that has all of the primers listed, from the command line, run `python compute_primers.py`. You will get a CSV file, named `all_primers.csv`, that will house the primers that you will need to order.


# Changelog

## Version 1.3

- API changes that allow for interfacing with a [browser-based GUI][1] that is being developed in parallel.
- Added standard influenza cloning plasmid backbones.

## Version 1.2

- Added a class that converts one nucleotide sequence into another, using Gibson assembly primers.

[1]: http://github.com/ericmjl/flu-gibson-webui/