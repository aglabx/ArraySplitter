# ArraySplitter: De Novo Decomposition of Satellite DNA Arrays

Decomposes satellite DNA arrays into monomers within telomere-to-telomere (T2T) assemblies. Ideal for analyzing centromeric and pericentromeric regions on monomeric level.

**Status:** In development. Optimized for 100Kb scale arrays; longer arrays will work but may take longer to process. Signigicanlty longer time.

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
time arraysplitter --input chr1.arrays.fa --output chr1.arrays
```

**Explanation**

* **`--input chr1.arrays.fa`:**  FASTA file of satDNA arrays.
* **`--output chr1.arrays`:** Output FASTA containing decomposed monomers (separated by spaces).

**All Options** 

```bash
arraysplitter --help 
```

## Contact

For questions or support: ad3002@gmail.com
