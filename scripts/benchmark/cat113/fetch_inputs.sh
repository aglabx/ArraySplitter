#!/usr/bin/env bash
# Build the three cat benchmark FASTAs from the NCBI GenBank assemblies.
#
#   fetch_inputs.sh OUTDIR
#
# Downloads F.catus_OcraZ1 pat1.0 (GCA_051911875.1) and mat1.0 (GCA_051911995.1) unless already in OUTDIR,
# checks them against NCBI's md5checksums.txt, cuts the arrays with samtools faidx and writes:
#   OUTDIR/cat7.fa   7-array panel (data/cat7_arrays.bed, all on pat1.0), header = label
#   OUTDIR/pat.fa    41 pat1.0 113-bp fields (data/fields_pat.tsv)
#   OUTDIR/mat.fa    51 mat1.0 113-bp fields (data/fields_mat.tsv)
# Field headers: NAME|ACCESSION:START-END|LENGTHbp|calls=N|calls_minus=M (START 1-based, as in the tsv).
# Sequence is upper case, 50 bases per line, so the files are byte-identical to the ones the paper used;
# the run stops if any MD5 differs from data/inputs.md5.
set -euo pipefail
O=${1:?outdir}
H=$(cd "$(dirname "$0")" && pwd)
for t in curl samtools md5sum gzip; do command -v $t > /dev/null || { echo "fetch_inputs: $t not on PATH" >&2; exit 1; }; done
mkdir -p "$O"
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/911
declare -A DIR=([pat]=875/GCA_051911875.1_F.catus_OcraZ1_pat1.0 [mat]=995/GCA_051911995.1_F.catus_OcraZ1_mat1.0)

for h in pat mat; do
  d=${DIR[$h]}; gz=$(basename "$d")_genomic.fna.gz
  [ -s "$O/$gz" ] || curl -fsSL -o "$O/$gz" "$NCBI/$d/$gz"
  curl -fsSL -o "$O/$h.md5checksums.txt" "$NCBI/$d/md5checksums.txt"
  want=$(awk -v f="./$gz" '$2 == f {print $1}' "$O/$h.md5checksums.txt")
  [ -n "$want" ] || { echo "fetch_inputs: $gz not in NCBI md5checksums.txt" >&2; exit 1; }
  got=$(md5sum "$O/$gz" | cut -c1-32)
  [ "$got" = "$want" ] || { echo "fetch_inputs: $gz md5 $got, NCBI says $want" >&2; exit 1; }
  [ -s "$O/$h.fna" ] || gzip -dc "$O/$gz" > "$O/$h.fna"
  [ -s "$O/$h.fna.fai" ] || samtools faidx "$O/$h.fna"
done

# one region -> upper case, 50 per line
cut1() { samtools faidx "$1" "$2" | grep -v '^>' | tr -d '\n' | tr 'acgtn' 'ACGTN' | fold -w 50; echo; }

: > "$O/cat7.fa"
while IFS=$'\t' read -r acc s e label; do
  [ -n "$acc" ] || continue
  echo ">$label" >> "$O/cat7.fa"
  cut1 "$O/pat.fna" "$acc:$((s + 1))-$e" >> "$O/cat7.fa"
done < "$H/data/cat7_arrays.bed"

for h in pat mat; do
  : > "$O/$h.fa"
  tail -n +2 "$H/data/fields_$h.tsv" | while IFS=$'\t' read -r name acc s e len calls minus; do
    echo ">$name|$acc:$s-$e|${len}bp|calls=$calls|calls_minus=$minus" >> "$O/$h.fa"
    cut1 "$O/$h.fna" "$acc:$s-$e" >> "$O/$h.fa"
  done
done

# fold leaves no empty line, but a sequence whose length is a multiple of 50 would get one from the echo
for f in cat7 pat mat; do sed -i '/^$/d' "$O/$f.fa"; done
(cd "$O" && md5sum -c "$H/data/inputs.md5") || { echo "fetch_inputs: built FASTA differs from the paper's inputs" >&2; exit 1; }
echo "fetch_inputs: cat7.fa pat.fa mat.fa written to $O and match data/inputs.md5"
