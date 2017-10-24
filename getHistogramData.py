#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""getHistogramData.py

This script can be used to summarize methylation percentages for a histogram.

The input file must have the first five columns (tab-separated):
chromosome
position
context
total coverage
percent methylated (0..100)

Output columns are:
percentageBin (0..100, by 1)
occurence
context
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

def getHistogramData(fileName):
    stats = dict([(x, [0]*101) for x in ["CpG", "CG", "CHG", "CHH"]])
    lineCounter = 0
    with open(fileName, "rb") as infile:
        for line in infile:
            lineCounter += 1
            if (lineCounter % 1000000) == 0:
                logging.info("Processed {0} million lines.".format(lineCounter/1000000))
            fields = line.decode("ascii").rstrip('\n').split('\t')
            context = fields[2]
            totCov = int(fields[3])
            pCentMeth = int(round(float(fields[4])))
            if totCov >= 5 and totCov <= 100:
                stats[context][pCentMeth] += 1
    for context, conStats in stats.items():
        occSum = sum(conStats)
        if occSum == 0:
            continue
        for i, v in enumerate(conStats):
            sys.stdout.write(str(i)+'\t'+str(v)+'\t'+context+'\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="getHistogramData.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to summarize methylation percentages for a histogram.

            The input file must have the first five columns (tab-separated):
            chromosome
            position
            context
            total coverage
            percent methylated (0..100)

            Output columns are:
            percentageBin (0..100, by 1)
            occurence
            context
            
            The results are printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("methFile", type=str,
                        help="""
                        A file with the methylation data.
                        """)
    
    args = parser.parse_args()
    
    getHistogramData(args.methFile)

