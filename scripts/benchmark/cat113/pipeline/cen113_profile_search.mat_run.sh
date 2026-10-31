#!/usr/bin/env bash
# Profile search for the cat 113-bp centromeric satellite on whole chromosomes (2026-10-30; PI: it must sit on every
# centromere that is assembled, «Да должен везде если собран»; its exact k-mers and alphabet's anchors found it on 6
# of 19 chromosomes only, and its copies are too diverged for TRF at satellome's setting).
#
#   cen113_profile_search.sh MONOMERS_TSV TWOBIT CHROM_LIST WORKDIR [THREADS]
#
# 1. Profiles: from MONOMERS_TSV (ArraySplitter monomers.tsv of the E1 and E2 arrays, 113 bp HOR / 107-110 bp
#    monomer), the base monomers of 95-130 bp of each array, every n-th to keep 200 per array (deterministic),
#    aligned with MAFFT (--auto) and built with hmmbuild --dna, one model per array (their strands are not
#    reconciled; nhmmer searches both strands). The sample size and length window are choices of this check.
# 2. nhmmer (E <= 1e-5, both strands) of both models against each chromosome in CHROM_LIST (one per line), cut
#    from TWOBIT.
# 3. WORKDIR/clusters.tsv: hits of either model merged when <= 5 kb apart; clusters >= 10 kb with their hit count
#    and the model with more hits. WORKDIR/per_chrom.tsv: per chromosome, hits, bp covered by hits, largest cluster.
# A chromosome whose nhmmer fails stops the run; a chromosome without hits gets a row with 0.
set -euo pipefail
M=${1:?monomers.tsv}; TB=${2:?2bit}; CL=${3:?chrom list}; W=${4:?workdir}; T=${5:-16}
for t in mafft hmmbuild nhmmer twoBitToFa; do command -v $t > /dev/null || { echo "cen113_profile_search: $t not on PATH" >&2; exit 1; }; done
mkdir -p "$W/chrom"
for a in E1_113_p7279 E2_113_p11906; do
  if [ ! -s "$W/$a.hmm" ]; then
    awk -F'\t' -v a="$a" 'NR > 1 && $1 == a && $2 == "base_monomer" && $4 >= 95 && $4 <= 130 {print $3 "\t" $17}' "$M" > "$W/$a.mono.tsv"
    n=$(wc -l < "$W/$a.mono.tsv"); [ "$n" -ge 200 ] || { echo "cen113_profile_search: $a has $n monomers of 95-130 bp" >&2; exit 1; }
    awk -v n="$n" -v k=200 'BEGIN {step = n / k} {i++; if (i >= next_pick) {printf ">m%s\n%s\n", $1, $2; next_pick += step}}' "$W/$a.mono.tsv" \
      | awk '/^>/ {c++} c <= 200' > "$W/$a.sample.fa"
    [ "$(grep -c '^>' "$W/$a.sample.fa")" -eq 200 ] || { echo "cen113_profile_search: $a sample is not 200" >&2; exit 1; }
    mafft --auto --thread "$T" --quiet "$W/$a.sample.fa" > "$W/$a.aln.fa"
    hmmbuild --dna -n "cat113_$a" "$W/$a.hmm" "$W/$a.aln.fa" > "$W/$a.hmmbuild.log"
  fi
done
cat "$W/E1_113_p7279.hmm" "$W/E2_113_p11906.hmm" > "$W/cat113.hmm"
while read -r c; do
  [ -n "$c" ] || continue
  [ -s "$W/chrom/$c.tbl" ] && continue
  # nhmmer (HMMER 3.3.2) guesses the target alphabet from the first residues and fails on a chromosome that starts
  # with the telomere (CCCTAA has no G), --dna or not; a 1-kb random ACGT record goes first, its hits are dropped
  { echo ">alphabet_guard"; awk 'BEGIN {srand(1); for (i = 0; i < 1000; i++) printf "%s", substr("ACGT", int(rand() * 4) + 1, 1); print ""}'
    twoBitToFa -seq="$c" "$TB" stdout; } > "$W/chrom/$c.fa"
  grep -c '^>' "$W/chrom/$c.fa" | grep -qx 2 || { echo "cen113_profile_search: could not cut $c" >&2; exit 1; }
  nhmmer --dna --cpu "$T" --noali -E 1e-5 --tblout "$W/chrom/$c.tbl.part" "$W/cat113.hmm" "$W/chrom/$c.fa" > "$W/chrom/$c.nhmmer.log" \
    || { echo "cen113_profile_search: nhmmer failed on $c" >&2; exit 1; }
  /bin/mv -f "$W/chrom/$c.tbl.part" "$W/chrom/$c.tbl"
  /bin/rm -f "$W/chrom/$c.fa"
done < "$CL"
# hits: chrom, start, end (both strands as low-high), model
printf 'chrom\tstart\tend\tlength\thits\ttop_model\n' > "$W/clusters.tsv"
printf 'chrom\thits\thit_bp\tclusters_ge_10kb\tlargest_cluster\n' > "$W/per_chrom.tsv"
while read -r c; do
  [ -n "$c" ] || continue
  awk -v c="$c" '!/^#/ && $1 != "alphabet_guard" {s = ($7 < $8) ? $7 : $8; e = ($7 < $8) ? $8 : $7; print c "\t" s "\t" e "\t" $3}' "$W/chrom/$c.tbl" \
    | sort -k2,2n > "$W/chrom/$c.hits.tsv"
  awk -F'\t' -v c="$c" -v out="$W/clusters.tsv" '
    function flush() {if (n == 0) return; L = ce - cs + 1; tm = ""; best = 0
                      for (m in mm) if (mm[m] > best) {best = mm[m]; tm = m}
                      if (L >= 10000) {printf "%s\t%d\t%d\t%d\t%d\t%s\n", c, cs, ce, L, n, tm >> out; nbig++
                                       if (L > big) {big = L; bigs = cs "-" ce}}
                      delete mm; n = 0}
    {if (n > 0 && $2 > ce + 5000) flush()
     if (n == 0) {cs = $2; ce = $3} else if ($3 > ce) ce = $3
     n++; mm[$4]++; tot++
     if ($2 > ue) {bp += $3 - $2 + 1; ue = $3} else if ($3 > ue) {bp += $3 - ue; ue = $3}}
    END {flush(); printf "%s\t%d\t%d\t%d\t%s\n", c, tot, bp, nbig, (big ? bigs " (" big " bp)" : "-")}' \
    "$W/chrom/$c.hits.tsv" >> "$W/per_chrom.tsv"
done < "$CL"
cat "$W/per_chrom.tsv"
