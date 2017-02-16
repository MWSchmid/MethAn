#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""getSimpleCoverageStats.py

This script can be used to get some coverage statistics per chromosome and context.

The input file must have the first five columns (tab-separated):
chromosome
position
context
total coverage
percent methylated (0..100)

Output columns are:
chrom
contect
counts
lowCounts
totCovSum
meCovSum
unCovSum
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')
import math

class coverageInfo():
    def __init__(self, chrom, ctxt, totCov, pCentMeth):
        self.chrom = chrom
        self.anyContext = ["CpG", "CG", "CHG", "CHH"]
        self.count = dict([(x, 0) for x in self.anyContext])
        self.lowCount = dict([(x, 0) for x in self.anyContext])
        self.highCount = dict([(x, 0) for x in self.anyContext])
        self.totCovSum = dict([(x, 0) for x in self.anyContext])
        self.meCovSum = dict([(x, 0) for x in self.anyContext])
        self.unCovSum = dict([(x, 0) for x in self.anyContext])
        self.addValue(ctxt, totCov, pCentMeth)
    
    def __str__(self):
        outLines = []
        for ctxt in self.anyContext:
            outLines.append('\t'.join([self.chrom, ctxt, str(self.lowCount[ctxt]+self.highCount[ctxt]+self.count[ctxt]), str(self.lowCount[ctxt]), str(self.highCount[ctxt]), str(self.count[ctxt]), str(self.totCovSum[ctxt]), str(self.meCovSum[ctxt]), str(self.unCovSum[ctxt])]))
        out = '\n'.join(outLines)
        return out
    
    def addValue(self, ctxt, totCov, pCentMeth):
        if totCov < 5:
            self.lowCount[ctxt] += 1
        else if totCov > 100:
            self.highCount[ctxt] += 1
        else:
            self.count[ctxt] += 1
            meCov = math.floor(totCov*(pCentMeth/100.0))
            self.totCovSum[ctxt] += totCov
            self.meCovSum[ctxt] += meCov
            self.unCovSum[ctxt] += totCov-meCov

def getCoverageStats(fileName):
    stats = {}
    lineCounter = 0
    with open(fileName, "rb") as infile:
        for line in infile:
            lineCounter += 1
            if (lineCounter % 1000000) == 0:
                logging.info("Processed {0} million lines.".format(lineCounter/1000000))
            fields = line.decode("ascii").rstrip('\n').split('\t')
            chrom = fields[0]
            ctxt = fields[2]
            totCov = int(fields[3])
            pCentMeth = float(fields[4])
            try:
                stats[chrom].addValue(ctxt, totCov, pCentMeth)
            except KeyError:
                stats[chrom] = coverageInfo(chrom, ctxt, totCov, pCentMeth)
    for chromName, covItem in stats.items():
        sys.stdout.write(str(covItem)+'\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="getSimpleCoverageStats.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to get some coverage statistics
            per chromosome and context.

            The input file must have the first five columns (tab-separated):
                chromosome
                position
                context
                total coverage
                percent methylated (0..100)
            
            Output columns are:
                chrom
                contect
                totalCounts
                lowCounts (cov <5)
                highCounts (cov >100)
                okCounts
                totCovSum
                meCovSum
                unCovSum
            
            The results are printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("methFile", type=str,
                        help="""
                        A file with the methylation data.
                        """)
    
    args = parser.parse_args()
    
    getCoverageStats(args.methFile)

