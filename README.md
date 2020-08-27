[![PyPI version](https://badge.fury.io/py/fasta.svg)](https://badge.fury.io/py/fasta)

# `fasta` version 2.1.8

This python package enables you to deal with biological sequence files easily.

The FASTA file format is a standard for storing several short or long DNA sequences inside a text file, use this package to manipulate and transform these files quickly and with short instructions.

## Prerequisites

Since `fasta` is written in python, it is compatible with all operating systems: Linux, macOS and Windows. The only prerequisite is `python3` (which is often installed by default) along with the `pip3` package manager.

To check if you have `python3` installed, type the following on your terminal:

    $ python3 -V

If you do not have `python3` installed, please refer to the section [obtaining python3](docs/installing_tips.md#obtaining-python3).

To check you have `pip3` installed, type the following on your terminal:

    $ pip3 -V

If you do not have `pip3` installed, please refer to the section [obtaining pip3](docs/installing_tips.md#obtaining-pip3).

## Installing

To install the `fasta` package, simply type the following commands on your terminal:

    $ pip3 install --user fasta

Alternatively, if you want to install it for all users of the system:

    $ sudo pip3 install fasta

## Usage

Bellow are some examples to illustrate the various ways there are to use this package.

Let's say you have a FASTQ file somewhere inside your home directory and you want to analyze it. To validate it, you can start by doing the following:

    >>> from fasta import FASTQ
    >>> fastq = FASTQ("~/repos/fasta/test/data/seqs.fastq")
    >>> print(fastq.validator())
    True
 
To check the number of reads inside the file, do the following:

    >>> from fasta import FASTQ
    >>> fastq = FASTQ("~/repos/fasta/test/data/seqs.fastq")
    >>> print(len(fastq))
    1401

Then, to run the FastQC software on that file automatically, do the following:

    >>> from fasta import FASTQ
    >>> fastq = FASTQ("~/repos/fasta/test/data/seqs.fastq")
    >>> print(fastq.fastqc())
    ~/repos/fasta/test/data/seqs.fastqc/

### Subsampling sequences

Next, to randomly pick a hundred sequences from the FASTQ file and put them in a new FASTQ file, use these commands:

    from fasta import FASTQ
    #TODO this example is not completed yet

### Working with forward and reverse reads

The `fasta` package also offers convenient ways of dealing with paired sequence files, where one has two FASTQ files with the same number of sequences in each file. Here is an example:

    from fasta import FASTQ
    #TODO this example is not completed yet

### Splitting FASTA files into sub-files

The `fasta` package also offers convenient ways of dealing with large number of sequences by automatically splitting them into an arbitrary number of smaller FASTA files. This is useful for the parallelization of certain operations. Here is an example:

    from fasta import FASTQ
    #TODO this example is not completed yet

### Parsing FASTA files with primers

The `fasta` package also offers functionality to parse reads from a FASTA file while automatically detecting the position of any forward and reverse primers, as well as the lack thereof. This is useful for filtering sequences and controlling quality. Here is an example:

    from fasta import FASTQ
    #TODO this example is not completed yet

### Producing visualizations

The `fasta` package is capable of producing certain types of graphs, such as a histogram of the sequence length distribution within a FASTA file:

    from fasta import FASTQ
    #TODO this example is not completed yet

### Others

The `fasta` package offers many other functions which have not been documented here yet. They can be discovered by looking at the source code or exploring the extra documentation below.

## Extra documentation 

More documentation is available at:

<http://xapple.github.io/fasta/fasta>

This documentation is simply generated from the source code with:

    $ pdoc --html --output-dir docs --force fasta