#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ArraySplitter Classify - Group arrays into families based on cut sequences and decomposition patterns.

The hypothesis: Arrays with the same cut sequence and similar decomposition patterns
likely belong to the same repeat family.
"""

import argparse
import os
from collections import Counter, defaultdict
from statistics import mean, median, stdev
import json

from tqdm import tqdm
import editdistance as ed


def extract_pattern_features_from_lengths_file(header, lengths_line):
    """
    Extract features from lengths file format.
    
    Args:
        header: Header line with metadata (e.g., ">seq_id cut=ATG orientation=fwd n_monomers=50 range=165-175 avg=171.2")
        lengths_line: Space-separated fragment lengths
        
    Returns:
        dict: Features including length distribution, variability metrics
    """
    # Parse header
    parts = header.split()
    seq_id = parts[0][1:]  # Remove '>'
    
    # Extract metadata from header
    metadata = {}
    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key] = value
    
    cut_seq = metadata.get('cut', '')
    n_monomers = int(metadata.get('n_monomers', 0))
    
    if not cut_seq or n_monomers == 0:
        return None
    
    # Parse lengths
    all_lengths = [int(x) for x in lengths_line.strip().split()]
    
    # Extract monomer lengths based on n_monomers count
    # The header tells us how many internal monomers there are
    # Typical structure: [left_flank] [monomer1] [monomer2] ... [monomerN] [right_flank]
    
    if n_monomers == len(all_lengths):
        # All are monomers, no flanks
        monomer_lengths = all_lengths
    elif n_monomers < len(all_lengths):
        # Determine which are monomers vs flanks
        # Use the average length from header to identify likely flanks
        avg_length = float(metadata.get('avg', 0))
        
        # Simple heuristic: monomers are closer to average length
        # First fragment is likely left flank if much different from average
        # Last fragment is likely right flank if much shorter than average
        
        if len(all_lengths) == n_monomers + 1:
            # One flank (either left or right)
            if abs(all_lengths[0] - avg_length) > abs(all_lengths[-1] - avg_length):
                # First is flank
                monomer_lengths = all_lengths[1:]
            else:
                # Last is flank
                monomer_lengths = all_lengths[:-1]
        elif len(all_lengths) == n_monomers + 2:
            # Both flanks present
            monomer_lengths = all_lengths[1:-1]
        else:
            # Fallback: take middle n_monomers elements
            start = (len(all_lengths) - n_monomers) // 2
            monomer_lengths = all_lengths[start:start + n_monomers]
    else:
        # Shouldn't happen, but handle gracefully
        monomer_lengths = all_lengths
    
    # Calculate statistics
    features = {
        'cut_sequence': cut_seq if cut_seq else 'UNKNOWN',
        'num_monomers': len(monomer_lengths),
        'mean_length': mean(monomer_lengths),
        'median_length': median(monomer_lengths),
        'min_length': min(monomer_lengths),
        'max_length': max(monomer_lengths),
        'length_variability': stdev(monomer_lengths) if len(monomer_lengths) > 1 else 0,
        'length_range': max(monomer_lengths) - min(monomer_lengths),
    }
    
    # Add length distribution
    length_counts = Counter(monomer_lengths)
    features['length_distribution'] = dict(length_counts)
    features['unique_lengths'] = len(length_counts)
    
    return features


def calculate_pattern_similarity(features1, features2):
    """
    Calculate similarity score between two decomposition patterns.
    
    Returns:
        float: Similarity score (0-1, where 1 is identical)
    """
    if not features1 or not features2:
        return 0.0
    
    # Different cut sequences = different families
    if features1['cut_sequence'] != features2['cut_sequence']:
        return 0.0
    
    scores = []
    
    # Compare mean lengths (most important)
    length_diff = abs(features1['mean_length'] - features2['mean_length'])
    mean_avg = (features1['mean_length'] + features2['mean_length']) / 2
    length_score = 1 - min(length_diff / mean_avg, 1.0)
    scores.append(length_score * 2)  # Weight this more
    
    # Compare variability
    var_diff = abs(features1['length_variability'] - features2['length_variability'])
    var_avg = (features1['length_variability'] + features2['length_variability']) / 2
    if var_avg > 0:
        var_score = 1 - min(var_diff / var_avg, 1.0)
        scores.append(var_score)
    
    # Compare length range
    range_diff = abs(features1['length_range'] - features2['length_range'])
    range_avg = (features1['length_range'] + features2['length_range']) / 2
    if range_avg > 0:
        range_score = 1 - min(range_diff / range_avg, 1.0)
        scores.append(range_score)
    
    # Compare number of unique lengths
    unique_diff = abs(features1['unique_lengths'] - features2['unique_lengths'])
    unique_avg = (features1['unique_lengths'] + features2['unique_lengths']) / 2
    if unique_avg > 0:
        unique_score = 1 - min(unique_diff / unique_avg, 1.0)
        scores.append(unique_score * 0.5)  # Weight this less
    
    return sum(scores) / len(scores) if scores else 0.0


def cluster_arrays(array_features, similarity_threshold=0.8):
    """
    Cluster arrays into families based on pattern similarity.
    
    Returns:
        dict: Family assignments {array_id: family_id}
    """
    # Group by cut sequence first
    cut_groups = defaultdict(list)
    for array_id, features in array_features.items():
        if features:
            cut_groups[features['cut_sequence']].append(array_id)
    
    # Within each cut group, cluster by pattern similarity
    family_assignments = {}
    family_id = 0
    
    for cut_seq, array_ids in cut_groups.items():
        if len(array_ids) == 1:
            # Single array with this cut
            family_assignments[array_ids[0]] = f"family_{family_id:04d}"
            family_id += 1
            continue
        
        # Cluster arrays with same cut by pattern similarity
        clusters = []
        
        for array_id in array_ids:
            assigned = False
            features = array_features[array_id]
            
            # Try to assign to existing cluster
            for cluster in clusters:
                # Compare with representative (first member)
                rep_id = cluster[0]
                rep_features = array_features[rep_id]
                
                similarity = calculate_pattern_similarity(features, rep_features)
                
                if similarity >= similarity_threshold:
                    cluster.append(array_id)
                    assigned = True
                    break
            
            if not assigned:
                # Create new cluster
                clusters.append([array_id])
        
        # Assign family IDs to clusters
        for cluster in clusters:
            current_family = f"family_{family_id:04d}"
            for array_id in cluster:
                family_assignments[array_id] = current_family
            family_id += 1
    
    return family_assignments


def classify_arrays(input_file, output_prefix, similarity_threshold=0.8, verbose=False):
    """
    Main classification function using lengths file.
    
    Args:
        input_file: Path to .lengths file from decomposition
        output_prefix: Prefix for output files
        similarity_threshold: Threshold for clustering (0-1)
        verbose: Verbose output
    """
    print("ArraySplitter Classify")
    print("="*60)
    print(f"Input file: {input_file}")
    print(f"Similarity threshold: {similarity_threshold}")
    
    # Check if input is a lengths file
    if not input_file.endswith('.lengths'):
        print("Warning: Input file should be a .lengths file from decomposition step")
    
    # Read lengths file
    array_features = {}
    array_info = {}
    
    print("\nReading lengths file...")
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Process pairs of lines (header + lengths)
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i].strip()
            if i + 1 < len(lines):
                lengths_line = lines[i + 1].strip()
                
                # Extract sequence ID
                seq_id = header.split()[0][1:]  # Remove '>'
                
                # Extract features
                features = extract_pattern_features_from_lengths_file(header, lengths_line)
                
                if features:
                    array_features[seq_id] = features
                    
                    # Parse additional info from header
                    parts = header.split()
                    metadata = {}
                    for part in parts[1:]:
                        if '=' in part:
                            key, value = part.split('=', 1)
                            metadata[key] = value
                    
                    # Calculate total length from all fragments
                    all_lengths = [int(x) for x in lengths_line.strip().split()]
                    total_length = sum(all_lengths)
                    
                    array_info[seq_id] = {
                        'length': total_length,
                        'cut_sequence': metadata.get('cut', ''),
                        'orientation': metadata.get('orientation', ''),
                        'n_monomers': int(metadata.get('n_monomers', 0)),
                        'avg_length': float(metadata.get('avg', 0))
                    }
                else:
                    if verbose:
                        print(f"Warning: Could not extract features for {seq_id}")
                
                i += 2
            else:
                i += 1
        else:
            i += 1
    
    print(f"Found {len(array_features)} arrays")
    print(f"\nSuccessfully analyzed {len(array_features)} arrays")
    
    # Cluster arrays
    print("\nClustering arrays into families...")
    family_assignments = cluster_arrays(array_features, similarity_threshold)
    
    # Count families
    family_counts = Counter(family_assignments.values())
    print(f"\nFound {len(family_counts)} families")
    
    # Write results
    output_file = f"{output_prefix}.families.tsv"
    stats_file = f"{output_prefix}.family_stats.tsv"
    json_file = f"{output_prefix}.features.json"
    
    print(f"\nWriting results to:")
    print(f"  {output_file}")
    print(f"  {stats_file}")
    print(f"  {json_file}")
    
    # Write family assignments
    with open(output_file, 'w') as f:
        f.write("array_id\tfamily\tcut_sequence\tlength\tmean_monomer_length\tnum_monomers\n")
        
        for array_id in sorted(family_assignments.keys()):
            family = family_assignments[array_id]
            info = array_info.get(array_id, {})
            features = array_features.get(array_id, {})
            
            f.write(f"{array_id}\t{family}\t{info.get('cut_sequence', 'NA')}\t")
            f.write(f"{info.get('length', 0)}\t")
            f.write(f"{features.get('mean_length', 0):.1f}\t")
            f.write(f"{features.get('num_monomers', 0)}\n")
    
    # Write family statistics
    with open(stats_file, 'w') as f:
        f.write("family\tnum_arrays\tcut_sequence\tmean_length\tstd_length\tmean_monomers\n")
        
        for family in sorted(family_counts.keys()):
            # Get all arrays in this family
            family_arrays = [aid for aid, fam in family_assignments.items() if fam == family]
            
            # Calculate family statistics
            cut_seq = array_features[family_arrays[0]]['cut_sequence']
            lengths = [array_info[aid]['length'] for aid in family_arrays]
            monomer_counts = [array_features[aid]['num_monomers'] for aid in family_arrays]
            
            f.write(f"{family}\t{len(family_arrays)}\t{cut_seq}\t")
            f.write(f"{mean(lengths):.1f}\t")
            f.write(f"{stdev(lengths) if len(lengths) > 1 else 0:.1f}\t")
            f.write(f"{mean(monomer_counts):.1f}\n")
    
    # Write detailed features as JSON
    with open(json_file, 'w') as f:
        output_data = {
            'array_features': array_features,
            'family_assignments': family_assignments,
            'array_info': array_info
        }
        json.dump(output_data, f, indent=2)
    
    # Print summary
    print("\nClassification summary:")
    print(f"  Total arrays: {len(monomers_by_array)}")
    print(f"  Successfully classified: {len(family_assignments)}")
    print(f"  Number of families: {len(family_counts)}")
    
    # Show top families
    print("\nTop 10 families by size:")
    for family, count in family_counts.most_common(10):
        family_arrays = [aid for aid, fam in family_assignments.items() if fam == family]
        cut_seq = array_features[family_arrays[0]]['cut_sequence']
        print(f"  {family}: {count} arrays (cut: {cut_seq})")


def run_it():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Classify satellite DNA arrays into families based on decomposition patterns"
    )
    parser.add_argument("-i", "--input", help="Input .monomers.tsv file from decomposition", required=True)
    parser.add_argument("-o", "--output", help="Output prefix", required=True)
    parser.add_argument(
        "-s", "--similarity", 
        help="Similarity threshold for clustering (0-1, default: 0.8)", 
        type=float, 
        default=0.8
    )
    parser.add_argument(
        "-v", "--verbose", 
        help="Verbose output", 
        action="store_true"
    )
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.input):
        print(f"Error: Input file {args.input} not found")
        exit(1)
    
    if not args.input.endswith('.monomers.tsv'):
        print(f"Warning: Input file should be a .monomers.tsv file from ArraySplitter decomposition")
    
    classify_arrays(
        args.input, 
        args.output, 
        similarity_threshold=args.similarity,
        verbose=args.verbose
    )


if __name__ == "__main__":
    run_it()