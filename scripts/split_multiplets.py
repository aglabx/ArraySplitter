#!/usr/bin/env python3
"""
Split multiplet segments using alignment to known monomers.

Strategy:
1. Use anchor positions to cut array into segments
2. Segments of ~1×P are monomers (templates)
3. Segments of ~2×P, 3×P, etc. are multiplets
4. Align monomer template(s) against multiplets to find internal boundaries
"""

import sys
import editdistance


def read_fasta_single(filepath):
    seq = []
    header = ""
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:]
            else:
                seq.append(line.strip().upper())
    return header, "".join(seq)


def find_kmer_positions(seq, kmer):
    positions = []
    k = len(kmer)
    for i in range(len(seq) - k + 1):
        if seq[i:i+k] == kmer:
            positions.append(i)
    return positions


def align_monomer_to_segment(monomer, segment, period):
    """
    Slide monomer along segment, compute ED at each position.
    Return list of (position, ed) for best non-overlapping placements.
    """
    mono_len = len(monomer)
    results = []

    for start in range(0, len(segment) - mono_len + 1):
        subseq = segment[start:start + mono_len]
        ed = editdistance.eval(monomer, subseq)
        results.append((start, ed))

    return results


def find_best_splits(segment, monomer, period, n_parts):
    """
    Find n_parts monomer boundaries within a segment.
    Uses sliding window ED to find positions where a monomer-sized
    piece best matches the template.
    """
    mono_len = len(monomer)
    seg_len = len(segment)

    # Compute ED at every position
    ed_profile = align_monomer_to_segment(monomer, segment, period)

    # Find n_parts non-overlapping minima, spaced ~period apart
    # Strategy: for each expected start position (0, P, 2P, ...),
    # search within ±window for the local ED minimum
    window = period // 4  # search ±25% of period

    boundaries = [0]  # first monomer always starts at 0
    for part_idx in range(1, n_parts):
        expected_pos = int(part_idx * seg_len / n_parts)
        search_start = max(0, expected_pos - window)
        search_end = min(len(ed_profile), expected_pos + window)

        best_pos = expected_pos
        best_ed = float('inf')
        for pos, ed in ed_profile[search_start:search_end]:
            if ed < best_ed:
                best_ed = ed
                best_pos = pos

        # Actually search in the ed_profile list
        best_pos = expected_pos
        best_ed = float('inf')
        for i in range(search_start, min(search_end, len(ed_profile))):
            pos, ed = ed_profile[i]
            if ed < best_ed:
                best_ed = ed
                best_pos = pos

        boundaries.append(best_pos)

    boundaries.append(seg_len)
    return boundaries


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.fasta> [period] [anchor]")
        sys.exit(1)

    header, seq = read_fasta_single(sys.argv[1])
    period = int(sys.argv[2]) if len(sys.argv) > 2 else 172
    anchor = sys.argv[3] if len(sys.argv) > 3 else "GAAAGAA"

    seq_len = len(seq)
    n_expected = seq_len // period

    print(f"Array: {header}")
    print(f"Length: {seq_len} bp")
    print(f"Period: {period} bp")
    print(f"Expected monomers: {n_expected}")
    print(f"Anchor: {anchor}")
    print()

    # Step 1: Find anchor positions
    positions = find_kmer_positions(seq, anchor)
    print(f"Anchor occurrences: {len(positions)}")

    # Step 2: Cut into segments
    # Add 0 and seq_len as boundaries
    cut_points = [0] + positions + [seq_len]
    # Remove duplicates and sort
    cut_points = sorted(set(cut_points))

    segments = []
    for i in range(len(cut_points) - 1):
        start = cut_points[i]
        end = cut_points[i + 1]
        seg_len_i = end - start
        n_mono = round(seg_len_i / period)
        if n_mono == 0:
            n_mono = 1
        segments.append({
            "start": start,
            "end": end,
            "length": seg_len_i,
            "n_mono": n_mono,
            "type": "monomer" if n_mono == 1 else f"multiplet_{n_mono}x",
        })

    # Classify segments
    monomers_1x = [s for s in segments if s["n_mono"] == 1]
    multiplets = [s for s in segments if s["n_mono"] > 1]

    print(f"\nSegments: {len(segments)} total")
    print(f"  1×P (monomers): {len(monomers_1x)}")
    for n in range(2, 6):
        mult_n = [s for s in segments if s["n_mono"] == n]
        if mult_n:
            print(f"  {n}×P (multiplets): {len(mult_n)}")

    # Step 3: Build consensus from 1×P monomers
    print(f"\n{'='*80}")
    print(f"Building template from {len(monomers_1x)} single-copy monomers")
    print(f"{'='*80}")

    mono_seqs = [seq[s["start"]:s["end"]] for s in monomers_1x]
    mono_lengths = [len(m) for m in mono_seqs]
    print(f"Monomer lengths: min={min(mono_lengths)}, max={max(mono_lengths)}, "
          f"mean={sum(mono_lengths)/len(mono_lengths):.1f}")

    # Use median-length monomer as template (or the one closest to period)
    mono_seqs_sorted = sorted(mono_seqs, key=lambda m: abs(len(m) - period))
    template = mono_seqs_sorted[0]
    print(f"Template length: {len(template)} bp (closest to period)")

    # Compute ED between all monomers
    eds = []
    for m in mono_seqs:
        ed = editdistance.eval(template, m)
        eds.append(ed)
    print(f"ED template vs monomers: min={min(eds)}, max={max(eds)}, "
          f"mean={sum(eds)/len(eds):.1f}")

    # Step 4: Split multiplets using template alignment
    print(f"\n{'='*80}")
    print(f"Splitting {len(multiplets)} multiplets using alignment")
    print(f"{'='*80}")

    all_monomers = []  # (start, end, source_type)

    for s in segments:
        if s["n_mono"] == 1:
            all_monomers.append((s["start"], s["end"], "anchor"))
        else:
            segment_seq = seq[s["start"]:s["end"]]
            n_parts = s["n_mono"]

            print(f"\n  Multiplet at {s['start']}-{s['end']} ({s['length']}bp, {n_parts}×P):")

            # Compute ED profile: slide template along segment
            mono_len = len(template)
            ed_profile = []
            for pos in range(len(segment_seq) - mono_len + 1):
                subseq = segment_seq[pos:pos + mono_len]
                ed = editdistance.eval(template, subseq)
                ed_profile.append((pos, ed))

            # Find n_parts local minima, spaced ~period apart
            window = period // 4
            boundaries = [0]

            for part_idx in range(1, n_parts):
                expected_pos = int(part_idx * len(segment_seq) / n_parts)
                search_start = max(0, expected_pos - window)
                search_end = min(len(ed_profile), expected_pos + window)

                best_pos = expected_pos
                best_ed = float('inf')
                for i in range(search_start, search_end):
                    pos, ed = ed_profile[i]
                    if ed < best_ed:
                        best_ed = ed
                        best_pos = pos

                boundaries.append(best_pos)
                print(f"    Split {part_idx}: expected={expected_pos}, "
                      f"found={best_pos} (ED={best_ed}, delta={best_pos - expected_pos:+d})")

            boundaries.append(len(segment_seq))

            # Add split monomers
            for i in range(len(boundaries) - 1):
                mono_start = s["start"] + boundaries[i]
                mono_end = s["start"] + boundaries[i + 1]
                all_monomers.append((mono_start, mono_end, f"split_{n_parts}x"))

    # Step 5: Summary
    print(f"\n{'='*80}")
    print(f"Final decomposition: {len(all_monomers)} monomers")
    print(f"{'='*80}")

    all_monomers.sort(key=lambda x: x[0])
    lengths = [end - start for start, end, _ in all_monomers]
    print(f"Lengths: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.1f}")
    print(f"Expected: {n_expected} monomers of ~{period}bp")
    print()

    print(f"{'idx':>4} {'start':>6} {'end':>6} {'length':>6} {'source':>10} {'ED_to_tmpl':>10}")
    print(f"{'-'*4} {'-'*6} {'-'*6} {'-'*6} {'-'*10} {'-'*10}")
    for i, (start, end, source) in enumerate(all_monomers):
        mono_seq = seq[start:end]
        ed = editdistance.eval(template, mono_seq)
        length = end - start
        flag = ""
        if length < period * 0.7 or length > period * 1.3:
            flag = " !!!"
        print(f"{i:>4} {start:>6} {end:>6} {length:>6} {source:>10} {ed:>10}{flag}")


if __name__ == "__main__":
    main()
