#!/usr/bin/env python
from __future__ import print_function
from maxent_sampling import maxent_motif, maxent_motifs
from maxent_sampling import spoof_maxent_motif, spoof_maxent_motifs
from uniform_sampling import uniform_motif, uniform_motifs            
from uniform_sampling import spoof_uniform_motif, spoof_uniform_motifs
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sample Motifs')
    parser.add_argument('-L', type=int, 
                   help='Number of sites in motif')
    parser.add_argument('-N', type=int,
                   help='Site length')
    parser.add_argument('-I', type=float,metavar='IC',
                   help='Desired IC')
    parser.add_argument('--in_fasta', type=str,metavar='FASTA INPUT',
                   help='Read motif from fasta file')
    parser.add_argument('--out_fasta', type=str,metavar='FASTA OUTPUT',
                   help='Write sampled motifs to this file')
    parser.add_argument('--TU', action='store_true',default='False',
                        help='Use Truncated Uniform distribution (Default MaxEnt)')
    parser.add_argument('--epsilon', type=float, default=0.1,
                        help='IC tolerance if using Truncated Uniform (Default 0.1 bits)')
    parser.add_argument('--num_motifs', type=int, default=1,
                   help='Number of motifs to sample')
    args = parser.parse_args()
    if args.in_fasta:
        motif = read_fasta(parser.in_fasta)
        if args.TU:
            motifs = spoof_uniform_motifs(motif,num_motifs=args.num_motifs)
        else:
            motifs = spoof_maxent_motifs(motif,num_motifs=args.num_motifs)
    else:
        N, L, IC, epsilon, num_motifs = (args.N, args.L, args.I,
                                         args.epsilon, args.num_motifs)
        if args.TU:
            motifs = uniform_motifs(N, L, IC, num_motifs, epsilon)
        else:
            motifs = uniform_motifs(N, L, IC, num_motifs)
    f = open(args.out_fasta,'w') if args.out_fasta else sys.stdout
    for i, motif in enumerate(motifs):
        for j, site in enumerate(motif):
            print(">Motif %s Site %s" % (i,j),file=f)
            print(site,file=f)
    f.close()
