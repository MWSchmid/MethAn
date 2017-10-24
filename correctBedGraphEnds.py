#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""correctBedGraphEnds.py

This script can be used to correct the ends in a bedGraph if they are falling off the chromosome.

"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

def correctBedGraph(bedGraphFile, chromSizesFile):
    chromSizes = {}
    with open(chromSizesFile, 'rb') as infile:
        for line in infile:
            (chrom, size) = line.decode("ascii").rstrip('\n').split('\t')
            chromSizes[chrom] = int(size)
    with open(bedGraphFile, "rb") as infile:
        for line in infile:
            (chrom, start, end, value) = line.decode("ascii").rstrip('\n').split('\t')
            if int(end) > chromSizes[chrom]:
                end = str(chromSizes[chrom])
            sys.stdout.write('\t'.join([chrom, start, end, value])+'\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="correctBedGraphEnds.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to correct the ends in
            a bedGraph if they are falling off the chromosome.
            
            LC_COLLATE=C
            correctBedGraphEnds.py $BEDGRAPH $CHROMSIZES | sort -k1,1 -k2,2n > $NEWBEDGRAPH
            
            The results are printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("bedGraph", type=str,
                        help="""
                        A file with the methylation data.
                        """)

    parser.add_argument("chromSizes", type=str,
                        help="""
                        A file with chrom ids and size (tab-sep).
                        """)

    args = parser.parse_args()
    
    correctBedGraph(args.bedGraph, args.chromSizes)

