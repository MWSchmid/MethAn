#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""getHistogramDataFromMappedFile.py

This script can be used to summarize methylation percentages for a histogram.

The input file must have the first seven columns (tab-separated):
chromosome
strand
position
methylated coverage
unmethylated coverage
context
featureID,featureClass|featureID,featureClass

(output from MethAnMap and MethAnDirectMap)

Output columns are:
percentageBin (0..100, by 1)
occurence
context
feature
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

def getHistogramDataFromMappedFile(fileName):
    stats = dict([(x, dict()) for x in ["CpG", "CG", "CHG", "CHH"]])
    lineCounter = 0
    with open(fileName, "rb") as infile:
        for line in infile:
            lineCounter += 1
            if (lineCounter % 1000000) == 0:
                logging.info("Processed {0} million lines.".format(lineCounter/1000000))
            fields = line.decode("ascii").rstrip('\n').split('\t')
            context = fields[5]
            meCov = int(fields[3])
            unCov = int(fields[4])
            totCov = float(meCov+unCov)
            pCentMeth = int(round(meCov/totCov*100))
            featurePart = fields[6]
            if totCov >= 5 and totCov <= 100:
                for curFeatureEntry in featurePart.split('|'):
                    curFeatureType = curFeatureEntry.split(',')[1]
                    if curFeatureType not in stats[context]:
                        stats[context][curFeatureType] = [0]*101
                    stats[context][curFeatureType][pCentMeth] += 1
    for context, conStats in stats.items():
        for featureType, featureStats in conStats.items():
            occSum = sum(featureStats)
            if occSum == 0:
                continue
            for i, v in enumerate(featureStats):
                sys.stdout.write(str(i)+'\t'+str(v)+'\t'+context+'\t'+featureType+'\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="getHistogramDataFromMappedFile.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to summarize methylation percentages for a histogram.

            The input file must have the first seven columns (tab-separated):
            chromosome
            strand
            position
            methylated coverage
            unmethylated coverage
            context
            featureID,featureClass|featureID,featureClass

            (output from MethAnMap and MethAnDirectMap)

            Output columns are:
            percentageBin (0..100, by 1)
            occurence
            context
            feature
            
            The results are printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("methFile", type=str,
                        help="""
                        A file with the methylation data.
                        """)
    
    args = parser.parse_args()
    
    getHistogramDataFromMappedFile(args.methFile)

