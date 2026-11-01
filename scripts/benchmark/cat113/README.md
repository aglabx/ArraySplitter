# Domestic cat 113-bp benchmark (paper §3.5, Tables 4-5, Figure 2)

Heterogeneous centromeric HORs of the diploid T2T domestic cat assembly F.catus_OcraZ1:
pat1.0 `GCA_051911875.1` and mat1.0 `GCA_051911995.1` (NCBI GenBank).

ArraySplitter v1.9.0 (regime router + phased 7-mer constellation) is compared with the single-anchor
baseline v1.8.2 (commit `433722a`), FasTAN (commit `963ac7f`) and centroAnno (commit `0b1f7dc`).

## Reproduce

```bash
# 1. inputs: downloads both assemblies from NCBI (~1.6 GB), checks NCBI MD5s, cuts the arrays,
#    and verifies the three FASTAs byte for byte against data/inputs.md5
bash fetch_inputs.sh inputs/

# 2. all runs, the Table 4-5 evaluator and the Figure 2 script
bash run_cat113.sh work/ inputs/ \
  /path/to/arraysplitter-1.9.0 /path/to/arraysplitter-1.8.2 \
  /path/to/FASTGA /path/to/FasTAN /path/to/centroAnno

# 3. the TRF columns of Table 4 (TRF 4.10.0)
bash run_trf.sh work/ inputs/ /path/to/trf
```

`work/eval_summary.txt` holds Tables 4 and 5, `work/figures/fig_cat113_heterogeneous_hor.pdf` is Figure 2, and
`work/MANIFEST.md5` pins every deterministic output. Requirements: `samtools`, `curl`, GNU `time`, `cargo`,
Python 3 with numpy and matplotlib, and `mafft` on `PATH` (Figure 2 aligns the most frequent monomers).

## Contents

| Path | What |
|---|---|
| `data/cat7_arrays.bed` | the 7-array panel on pat1.0: six 113-bp HOR arrays and the 483-bp `A3_483_control` (0-based BED, label = FASTA header) |
| `data/fields_pat.tsv`, `data/fields_mat.tsv` | the 41 pat1.0 and 51 mat1.0 113-bp satellite fields (1-based, inclusive) with the number of nhmmer monomer calls in each |
| `data/labels7.txt`, `data/inputs.md5` | panel order; MD5s of `cat7.fa`, `pat.fa`, `mat.fa` |
| `profiles/pass1_pat/`, `profiles/pass1_mat/` | first-pass nhmmer profiles (HMMER 3.3.2), built from ArraySplitter monomers of the E1 and E2 arrays of each haplotype |
| `profiles/pass2_pat/` | second-pass per-chromosome profiles, pat1.0 only |
| `reference/calls_{pat,mat}.bed`, `reference/groups_{pat,mat}.tsv` | the nhmmer monomer calls and the arrays they define. pat1.0 calls are two-pass, mat1.0 calls are first-pass only |
| `pipeline/` | the scripts that built the profiles, the calls and the fields (`cen113_profile_search.sh` → `cen113_monomer_table.sh` → `cen113_second_pass.sh` (pat) → `cen113_array_coverage.sh` → `cen113_fields_fasta.sh`; `trf_on_arrays.sh` cut the 7-array panel). They read a 2bit of each assembly; their headers give the arguments |
| `fetch_inputs.sh`, `run_cat113.sh` | inputs and runs |
| `run_trf.sh`, `trf_period_summary.py` | TRF at satellome's parameters (`2 5 7 80 10 50 2000`) and permissive ones (`2 3 3 75 20 30 2000`), and the per-period coverage summary behind Table 4 |
| `eval/` | Rust evaluator for Tables 4-5 (`cat113_constellation_eval WORKDIR ONEview`); writes `table5.tsv` with one definition per column for every tool (longest period >= 300 bp, a 108-116-bp unit, both nested) |
| `plot_cat113_fields.py` | Figure 2 (`RUN_DIR TABLE5_TSV OUT_DIR`); every number in a panel title is computed, panel D is `table5.tsv` |

The nhmmer reference is an annotation with its own model, not ground truth: Table 5 reports monomer counts as
a share of it, and a share above 100% means ArraySplitter emits more 108-116-bp units than nhmmer calls.
