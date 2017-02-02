#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""summarizeMethylationByGene.py

This script can be used to summarize methylation info by gene.

Takes the mapped nucleotides from MethAnMap or MethAnDirectMap.

Prints to std::out.

"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import operator
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

class nucCounter:
    
    def __init__(self, name, hierarchy, totCov, mePerc):
        self.name = name
        self.feature = hierarchy.pop(0)
        (self.num, self.cov, self.met) = (0, 0, 0)
        self.children = {}
        self.addNucleotide(hierarchy, totCov, mePerc)
    
    def __str__(self):
        out = '\n'.join([str(child) for child in self.children])
        out += '\n' + '\t'.join([self.name, self.feature, str(self.num), str(self.cov), str(self.met/self.num)])
        return out
   
    def addNucleotide(self, hierarchy, totCov, mePerc):
        self.num += 1
        self.cov += totCov
        self.met += mePerc
        if len(hierarchy) > 0:
            curLevel = hierarchy.pop()
            try:
                self.children[curLevel].addNucleotide(hierarchy, totCov, mePerc)
            except KeyError:
                self.children[curLevel] = nucCounter(self.name, hierarchy, totCov, mePerc)

def summarizeByGene(mappedNucleotides):
    geneDict = {}
    with open(mappedNucleotides, 'rb') as infile:
        for line in infile:
            fields = line.decode("ascii").rstrip('\n').split('\t')
            (chrom, strand, pos, umeCov, meCov, ctxt, mapping) = fields[:7]
            umeCov = int(umeCov)
            meCov = int(meCov)
            totCov = umeCov+meCov
            if totCov == 0:
                continue
            mePerc = (meCov/totCov)*100
            allGenes = list(set([x.split(',') for x in mapping.split('|')]))
            for curGeneName,curGeneFeature in allGenes:
                if curGeneName == "none":
                    continue
                try:
                    geneDict[curGeneName].addNucleotide(["topLevel", curGeneFeature, ctxt], totCov, mePerc)
                except KeyError:
                    geneDict[curGeneName] = nucCounter(curGeneName, ["topLevel", curGeneFeature, ctxt], totCov, mePerc)
    for curGeneName, curGeneName in geneDict.items():
        sys.stdout.write(str(curGeneName) + '\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="summarizeMethylationByGene.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used to summarize methylation info by gene.

            Takes the mapped nucleotides from MethAnMap or MethAnDirectMap.

            Prints to stdout.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("mappedNucleotides", nargs='+', type=str,
                        help="""
                        A file with mapped nucleotides. Output from
                        MethAnMap or MethAnDirectMap.
                        """)
    
    args = parser.parse_args()
    
    summarizeByGene(mappedNucleotides)
