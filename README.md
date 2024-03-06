# ArraySplitter: De Novo Decomposition of Satellite DNA Arrays

Decomposes satellite DNA arrays into monomers within telomere-to-telomere (T2T) assemblies. Ideal for analyzing centromeric and pericentromeric regions on monomeric level.

**Status:** In development. Optimized for 100Kb scale arrays; longer arrays will work but may take longer to process. Signigicanlty longer time.

**Update:** From 1.1.6, ArraySplitter now successfully decomposes arrays on megabase scale. Largest arrays takes around 5 minutes to process. Fortunatelly, there are only 41 arrays large 1 Mb in CHM13v20 assembly. And I'm going to add parallel processing to speed up singificantly the process. Currently, it is single-threaded.

**Update:** Monomers are required some polising of borders, I am working on it.

**Update:** To test ArraySplitter, I used CHM13v20 assembly, it takes around 3 hours, to decompose all arrays longer than 1 Kb (13K arrays).

## Installation

**Prerequisites**

* Python 3.6 or later

**Installation with pip:**

```bash
pip install arraysplitter
```

## Usage

**Basic Example**

```bash
time arraysplitter -i chr1.arrays.fa -o chr1.arrays
```

It will create a FASTA file with monomers separated by spaces.

**Explanation**

* **`-i chr1.arrays.fa`:**  FASTA file of satDNA arrays.
* **`-o chr1.arrays`:** Prefix for the output FASTA containing decomposed monomers (separated by spaces).

**All Options** 

```bash
arraysplitter --help 
```

## Rotating monomers to start with the same sequence

We found that different arrays of the same repeat family can be decomposed sligtly differently. To make them comparable, ArraySplitter can rotate monomers to start with the same sequence. 

```bash
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa
```

And you can give the sequence to start with:

```bash
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa -s TTTC
```

**Explanation**

* **`-i arrays.fa`:**  FASTA file of monomers.
* **`-o arrays.norm.fa`:** Output FASTA file with rotated monomers.

## Extracting and counting monomers

And finally, you can extract and count monomers from the arrays:

```bash
arraysplitter_extract -i arrays.norm.fa -o arrays.norm
```

It will create a file with monomer length, monomer frequency, and monomer sequence (ordered by frequency). For example, for the arrays.norm.fa file above, the output will be like this:

```bash
514     10      ATCCCATTCC
514     10      GATTGGAGTG
514     6       TCCTTT
514     5       TGCTG
514     10      ATTGAATGGA
514     10      ATGCAATGGA
514     5       TCCTA
```


## Contact

For questions or support: ad3002@gmail.com
