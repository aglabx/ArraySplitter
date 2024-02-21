# ArraySplitter
ArraySplitter: De Novo Decomposition of Satellite DNA Arrays into Monomers within Telomere-to-Telomere Assemblies

## Introduction

ArraySplitter is a tool for the de novo decomposition of satellite DNA arrays into monomers within telomere-to-telomere assemblies. It is designed to work with a very long satDNA arrays from T2T assemblies, such as those found in centromeric and pericentromeric regions. ArraySplitter is implemented in Python and is available as a standalone tool.

## Installation

ArraySplitter is implemented in Python and requires Python 3.6 or later. 

To install ArraySplitter, clone the repository from GitHub and install the required dependencies using pip:

```bash
pip install -r requirements.txt

python setup.py install
```

## Usage

ArraySplitter is a command-line tool. To see the available options, run:

```bash
time python ./ArraySplitter/main.py --help
```

The main input to ArraySplitter is a FASTA file containing the telomere-to-telomere assembly. The output is a FASTA file containing the monomers of the satellite DNA arrays seperated by spaces.

## Example

To run ArraySplitter on the provided test data, run:

```bash
time python ./ArraySplitter/main.py --input ./test_data/chr1.fa --output ./test_data/chr1_monomers
```

