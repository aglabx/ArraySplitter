#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse

from .core_functions.tools.rotation import rotate_arrays

from .core_functions.io.fasta_reader import sc_iter_fasta_file

def main(fasta_file, output_file, starting_kmer=None):

    arrays = []
    for ii, (header, seq) in enumerate(sc_iter_fasta_file(fasta_file)):
        arrays.append(seq)

    new_arrays = rotate_arrays(arrays, starting_kmer=starting_kmer)
    
    with open(output_file, "w") as fw:
        for ii, (header, seq) in enumerate(sc_iter_fasta_file(fasta_file)):
            fw.write(f'>{header}\n{" ".join(new_arrays[ii])}\n')
            
def run_it():
    parser = argparse.ArgumentParser(description='Rotate arrays to start from the same position.')
    parser.add_argument('-i', '--fasta', type=str, help='Fasta file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file', required=True)
    parser.add_argument('-s', '--start', type=str, help='Starting kmer for rotation [None]', default=None)
    args = parser.parse_args()
    main(args.fasta, args.output, starting_kmer=args.start)

if __name__ == '__main__':
    run_it()