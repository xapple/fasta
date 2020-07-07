[![PyPI version](https://badge.fury.io/py/fasta.svg)](https://badge.fury.io/py/fasta)

# `fasta` version 2.0.1

This python package enables you to deal with biological sequence files easily.

The FASTA file format is a standard for storing several short or long DNA sequences inside a text file, use this package to manipulate and transform these files quickly and with short instructions.

## Prerequisites

Since `fasta` is written in python it is compatible with all operating systems: Linux, macOS and Windows. The only prerequisite is `python3` which is often installed by default and the `pip3` package manager.

To check you have `python3` installed, type the following on your terminal:

    $ python3 -V

If you do not have `python3` installed, please refer to the section [getting python](docs/installing_tips.md#obtaining-python3).

To check you have `pip3` installed, type the following on your terminal:

    $ pip3 -V

If you do not have `pip3` installed, please refer to the section [getting pip](docs/installing_tips.md#obtaining-pip3).

## Installing

To install the `fasta` package, simply type the following commands on your terminal:

    $ pip3 install --user waste_flow

Alternatively, if you want to install it for all users of the system:

    $ sudo pip3 install waste_flow

## Usage

Bellow are some examples to illustrate the various ways there are to use this package.

To xyz you can do the following:

    from fasta import FASTQ
    print(FASTQ)
    
# Extra documentation 

More documentation is available at:

<http://xapple.github.io/fasta/fasta>

This documentation is simply generated from the source code with:

    $ pdoc --html --output-dir docs --force fasta