#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import Counter
import argparse
import logging
from .core_functions.io.fasta_reader import sc_iter_fasta_file

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main(fasta_file, output_prefix, distance=0, function='hd'):

    if distance < 0:
        raise ValueError(f'Distance should be >= 0, but {distance} given')
    
    if function not in ['hd', 'ed']:
        raise ValueError(f'Function should be hd or ed, but {function} given')
    
    logging.info(f'Extract monomers from {fasta_file}')
    monomers = Counter()
    for header, seq in sc_iter_fasta_file(fasta_file):
        for monomer in seq.split():
            monomers[monomer] += 1

    logging.info(f'Found {len(monomers)} uniq monomers with total {sum(monomers.values())} monomers')

    logging.info(f'Filter monomers with distance <= {distance} using {function} function')
    # if function == 'hd':
    #     monomers = {k: v for k, v in monomers.items() if v <= distance}
    # else:
    #     raise NotImplementedError('Edit distance is not implemented yet')
    logging.info(f'Found {len(monomers)} monomers with distance <= {distance}')

    output_file = f'{output_prefix}.monomers.tsv'
    logging.info(f'Save monomers to {output_file}')
    with open(output_file, "w") as fw:
        for monomer, count in monomers.most_common(1000000000):
            fw.write(f'{count}\t{len(monomer)}\t{monomer}\n')
            
def run_it():
    parser = argparse.ArgumentParser(description='Extract monomers from array over edit distance.')
    parser.add_argument('-i', '--fasta', type=str, help='Fasta file with arrays', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output prefix', required=True)
    parser.add_argument('-f', '--function', type=str, help='Distance function hamming (hd) or edit (ed) [hd]', default='hd')
    parser.add_argument('-d', '--distance', type=int, help='Distance value, you will get monomers over <= distance [0]', default=0)
    args = parser.parse_args()
    main(args.fasta, args.output, distance=args.distance, function=args.function)

if __name__ == '__main__':
    run_it()