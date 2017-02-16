#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mergeGFF.py

This script can be used to merge GFF files. The script assumes that the original GFFs are sorted correctly.

Merging means that the rows of the original GFFs are combined into one file and correctly sorted. Features are NOT merged!
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

class genomeFeature:
    
    def __init__(self, line):
        (self.chrom, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, self.rest) = line.rstrip('\n').split('\t')
        self.start = int(self.start)
        self.rest = self.rest.rstrip(';')
        self.restDict = dict([x.split('=') for x in self.rest.split(';')])
        self.children = []
    
    def __str__(self):
        outList = ['\t'.join([self.chrom, self.source, self.feature, str(self.start), self.end, self.score, self.strand, self.phase, self.rest])]
        outList.extend([str(child) for child in self.children])
        return '\n'.join(outList)
    
    def addPotentialChild(self, toAdd):
        if not "ID" in self.restDict.keys(): # exons sometimes have no ID - and they won't have children as well
            for okWithoutID in ["exon", "UTR", "CDS"]:
                if okWithoutID in self.feature:
                    return False
            logging.warning("Expected an ID for this entry: " + str(self))
            return False
        if not "Parent" in toAdd.restDict.keys():
            return False
        for child in self.children:
            if child.addPotentialChild(toAdd):
                return True
        if toAdd.restDict["Parent"] == self.restDict["ID"]:
            self.children.append(toAdd)
            return True
        return False

def loadGFF(infileName):
    out = []
    with open(infileName, 'rb') as infile:
        line = infile.readline().decode("ascii")
        while line[0] == "#":
            line = infile.readline().decode("ascii")
        curFeature = genomeFeature(line)
        out.append(curFeature)
        for line in infile:
            curFeature = genomeFeature(line.decode("ascii"))
            if not out[-1].addPotentialChild(curFeature):
                out.append(curFeature)
    return out

def mergeGFFs(fileList):
    allEntries = []
    for infileName in fileList:
        allEntries.extend(loadGFF(infileName))
    byChrom = {}
    for entry in allEntries:
        try:
            byChrom[entry.chrom].append(entry)
        except KeyError:
            byChrom[entry.chrom] = [entry]
    for chrom in sorted(byChrom.keys()):
        byChrom[chrom].sort(key=operator.attrgetter("start"))
        for entry in byChrom[chrom]:
            sys.stdout.write(str(entry) + '\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="mergeGFF.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to merge GFF files. The script
            assumes that the original GFFs are sorted correctly.
            
            Merging means that the rows of the original GFFs are
            combined into one file and correctly sorted.
            
            Features will NOT be merged!
            
            The resulting GFF is printed to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("gffFile", nargs='+', type=str,
                        help="""
                        Two or more GFF files. The file may contain headers
                        specified with #.
                        """)
    
    args = parser.parse_args()
    
    if len(args.gffFile) > 1:
        mergeGFFs(args.gffFile)
    else:
        logging.critical("You did not provide two or more GFFs.")
        sys.exit(66)

