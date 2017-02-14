#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""binMethylationData.py

This script can be used to summarize methylation data into larger windows.

The input file must have the first five columns (tab-separated):
chromosome
position
context
total coverage
percent methylated (0..100)

Output columns are:
chrom
start
end
context
percent methylated
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

class methBin():
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.anyContext = ["CpG", "CG", "CHG", "CHH"]
        self.count = dict([(x, 0) for x in self.anyContext])
        self.meth = dict([(x, 0) for x in self.anyContext])
    
    def __str__(self):
        outLines = []
        for ctxt in self.anyContext:
            if self.count[ctxt] > 0:
                aveMeth = self.meth[ctxt]/self.count[ctxt]
            else:
                aveMeth = "NA"
            outLines.append('\t'.join([self.chrom, str(self.start), str(self.end), ctxt, str(aveMeth)]))
        out = '\n'.join(outLines)
        return out
    
    def hasValues(self):
        for ctxt in self.anyContext:
            if self.count[ctxt] > 0:
                return True
        return False
    
    def tryToAddPosition(self, chrom, position, ctxt, percMeth):
        if (chrom != self.chrom) or (position < self.start) or (position >= self.end): # half open
            return False
        self.count[ctxt] += 1
        self.meth[ctxt] += percMeth
        return True

def binMethylationData(fileName, minCov, binSize):
    curBin = methBin("none", 0, 0)
    lineCounter = 0
    with open(fileName, "rb") as infile:
        for line in infile:
            lineCounter += 1
            if (lineCounter % 1000000) == 0:
                logging.info("Processed {0} million lines.".format(lineCounter/1000000))
            fields = line.decode("ascii").rstrip('\n').split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ctxt = fields[2]
            totCov = int(fields[3])
            if totCov < minCov:
                continue
            pCentMeth = float(fields[4])
            if curBin.tryToAddPosition(chrom, pos, ctxt, pCentMeth):
                continue
            if curBin.hasValues():
                sys.stdout.write(str(curBin)+'\n')
            newStart = math.floor(pos/float(binSize))*binSize
            curBin = methBin(chrom, newStart, newStart+binSize)
            if not curBin.tryToAddPosition(chrom, pos, ctxt, pCentMeth):
                logging.critical("This should not happen...")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="binMethylationData.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to summarize methylation data into larger windows.

            The input file must have the first five columns (tab-separated):
                chromosome
                position
                context
                total coverage
                percent methylated (0..100)
            
            IMPORTANT: THE FILE SHOULD BE SORTED BY CHROMOSOME AND POSITION!!!

            Output columns are:
                chrom
                start
                end
                context
                percent methylated
            
            IMPORTANT: THE LAST BINS MAY FALL OFF THE REFERENCE (the script does
            not know about the size of the chromosome).
            
            The results are printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("methFile", type=str,
                        help="""
                        A file with the methylation data.
                        """)

    parser.add_argument("-c", "--minCov", type=int, required=False, default = 5,
                        help="""
                        Minimal coverage at a position to be considered (default: 5).
                        """)

    parser.add_argument("-s", "--binSize", type=int, required=False, default = 1000,
                        help="""
                        Size of the bin (default: 1000).
                        """)
    
    args = parser.parse_args()
    
    binMethylationData(args.methFile, args.minCov, args.binSize)

