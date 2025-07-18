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
import csv

from tqdm import tqdm
import editdistance as ed


def extract_pattern_features_from_monomers(monomers_data):
    """
    Extract features from already decomposed monomers data.
    
    Args:
        monomers_data: List of monomer records for a single array
        
    Returns:
        dict: Features including length distribution, variability metrics
    """
    # Get monomer lengths (excluding flanks)
    monomer_lengths = []
    cut_seq = None
    
    for monomer in monomers_data:
        if monomer['type'] == 'monomer' and monomer['is_flank'] == 'False':
            monomer_lengths.append(monomer['length'])
            # Extract cut sequence from the first internal monomer
            if cut_seq is None and monomer['sequence']:
                # Find the common prefix among all monomers to determine cut sequence
                cut_seq = monomer['sequence'][:10]  # Start with first 10 chars
    
    if not monomer_lengths or not cut_seq:
        return None
    
    # Refine cut sequence by finding common prefix
    sequences = [m['sequence'] for m in monomers_data 
                 if m['type'] == 'monomer' and m['is_flank'] == 'False' and m['sequence']]
    if sequences:
        # Find longest common prefix
        cut_seq = sequences[0]
        for seq in sequences[1:]:
            # Find common prefix
            i = 0
            while i < len(cut_seq) and i < len(seq) and cut_seq[i] == seq[i]:
                i += 1
            cut_seq = cut_seq[:i]
            if not cut_seq:
                break
    
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
    Main classification function.
    
    Args:
        input_file: Path to .monomers.tsv file from decomposition
        output_prefix: Prefix for output files
        similarity_threshold: Threshold for clustering (0-1)
        verbose: Verbose output
    """
    print("ArraySplitter Classify")
    print("="*60)
    print(f"Input file: {input_file}")
    print(f"Similarity threshold: {similarity_threshold}")
    
    # Read monomers data
    monomers_by_array = defaultdict(list)
    array_lengths = {}
    
    print("\nReading monomers data...")
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            array_id = row['sequence_id']
            monomers_by_array[array_id].append({
                'orientation': row['orientation'],
                'index': int(row['index']),
                'type': row['type'],
                'length': int(row['length']),
                'is_flank': row['is_flank'],
                'sequence': row['sequence']
            })
            # Calculate total array length
            if array_id not in array_lengths:
                array_lengths[array_id] = 0
            array_lengths[array_id] += int(row['length'])
    
    print(f"Found {len(monomers_by_array)} arrays")
    
    # Extract features from decomposed data
    array_features = {}
    array_info = {}
    
    print("\nAnalyzing arrays...")
    for array_id, monomers in tqdm(monomers_by_array.items()):
        try:
            # Sort monomers by index to ensure correct order
            monomers.sort(key=lambda x: x['index'])
            
            # Extract features
            features = extract_pattern_features_from_monomers(monomers)
            
            if features:
                array_features[array_id] = features
                array_info[array_id] = {
                    'length': array_lengths[array_id],
                    'cut_sequence': features['cut_sequence'],
                    'orientation': monomers[0]['orientation'] if monomers else 'unknown'
                }
            else:
                if verbose:
                    print(f"Warning: Could not extract features for {array_id}")
        except Exception as e:
            if verbose:
                print(f"Error processing {array_id}: {e}")
            continue
    
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