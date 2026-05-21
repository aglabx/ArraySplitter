# Regression baseline

Locks the byte-exact output of `arraysplitter` on a fixed input panel so
intentional schema changes can be distinguished from accidental drift.

## What's here

| File | Role |
|---|---|
| `zfinch_iter2.manifest` | Locked md5/size/line-count for the 5 outputs on `test_data/zebra_finch_satdna.fasta`. Header comments record the binary SHA, command, and lock date. |
| `run.sh` | Builds the release binary, runs it on the input, and diffs the 5 outputs against the manifest. Exits 0 on full match, 1 otherwise. |

## How to run

```bash
bash tests/regression/run.sh             # build + run + check
bash tests/regression/run.sh --no-build  # skip cargo build (binary already present)
bash tests/regression/run.sh --keep      # don't delete the run output dir
```

Runs land under `results/regression_runs/zfinch_<timestamp>_<pid>/` (deleted
on exit unless `--keep`). `cargo build --release` is run from
`src/rust/arraysplitter` so `target/` stays where the binary expects it.

The script tolerates both macOS (`md5`, `stat -f%z`) and Linux (`md5sum`,
`stat -c%s`).

## When the manifest needs to change

The manifest is **load-bearing**: a failing `run.sh` after a commit that
wasn't supposed to change output is a regression and must be investigated.

When a commit **intentionally** changes output schema (new column, new row
type, semantics change), update the manifest **in the same commit** that
introduces the schema change:

1. Build the new binary.
2. Run it once: `./target/release/arraysplitter -i test_data/zebra_finch_satdna.fasta -o /tmp/new -t 4 --method autocorr`.
3. Compute new md5/size/lines for each of the 5 outputs.
4. Replace the rows in `zfinch_iter2.manifest` (and bump the comment header to mention the new commit SHA and what changed).
5. Re-run `bash tests/regression/run.sh` to confirm green.
6. Commit manifest + code together. Reviewers should be able to read the
   commit message and understand why the schema change was deliberate.

If only some files change schema, only those rows in the manifest need
updating; the others stay locked.

## Why locked outputs and not unit tests

These five files are the user-facing contract. A unit test on
`fast_edit_distance_bounded` doesn't catch a writer that silently reorders
rows; a unit test on the writer doesn't catch a sort that breaks when a new
type is added. Byte-comparison on a realistic 60 Mbp panel does both.

## Cost

On the locked iter-2 binary the run takes ~4-7 min depending on host (laptop
4-core: ~6:30; AMD Ryzen 5 3600 12-core: ~3:30). Build is ~30 s incremental,
~30 s clean release.
