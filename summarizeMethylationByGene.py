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
        self.feature = ""
        (self.num, self.cov, self.met) = (0, 0, 0)
        self.children = {}
        self.addNucleotide(hierarchy, totCov, mePerc)
    
    def __str__(self):
        out = [str(self.children[child])+'\t'+self.feature for child in self.children]
        if self.feature == "topLevel":
            out.append('\t'.join([self.name, self.feature, str(self.num), str(self.cov), str(self.met/self.num), "self"]))
        else:
            out.append('\t'.join([self.name, self.feature, str(self.num), str(self.cov), str(self.met/self.num)]))
        return '\n'.join(out)
   
    def addNucleotide(self, hierarchy, totCov, mePerc):
        self.feature = hierarchy.pop(0)
        self.num += 1
        self.cov += totCov
        self.met += mePerc
        if len(hierarchy) > 0:
            curLevel = hierarchy[0]
            try:
                self.children[curLevel].addNucleotide(hierarchy, totCov, mePerc)
            except KeyError:
                self.children[curLevel] = nucCounter(self.name, hierarchy, totCov, mePerc)

class nucCounterAlternative:
    
    def __init__(self, name, feature, context, totCov, mePerc):
        self.name = name
        self.featureAndContext = {}
        self.addNucleotide(feature, context, totCov, mePerc)
    
    def __str__(self):
        out = ['\t'.join([self.name, featureContextPair, str(vals[0]), str(vals[1]/vals[0]), str(vals[2]/vals[0])]) for featureContextPair, vals in self.featureAndContext.items()]
        return '\n'.join(out)
   
    def addNucleotide(self, feature, context, totCov, mePerc):
        featureContextPair = '\t'.join([feature, context])
        try:
            self.featureAndContext[featureContextPair][0] += 1
            self.featureAndContext[featureContextPair][1] += totCov
            self.featureAndContext[featureContextPair][2] += mePerc
        except KeyError:
            self.featureAndContext[featureContextPair] = [1, totCov, mePerc]

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
            allGenes = [x.split(',') for x in mapping.split('|')]
            for curGeneName,curGeneFeature in allGenes:
                if curGeneName == "none":
                    continue
                try:
                    #geneDict[curGeneName].addNucleotide(["topLevel", curGeneFeature, ctxt], totCov, mePerc)
                    geneDict[curGeneName].addNucleotide(curGeneFeature, ctxt, totCov, mePerc)
                except KeyError:
                    #geneDict[curGeneName] = nucCounter(curGeneName, ["topLevel", curGeneFeature, ctxt], totCov, mePerc)
                    geneDict[curGeneName] = nucCounterAlternative(curGeneName, curGeneFeature, ctxt, totCov, mePerc)
                    
    for curGeneName, curGeneItem in geneDict.items():
        sys.stdout.write(str(curGeneItem) + '\n')

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
    
    parser.add_argument("mappedNucleotides", type=str,
                        help="""
                        A file with mapped nucleotides. Output from
                        MethAnMap or MethAnDirectMap.
                        """)
    
    args = parser.parse_args()
    
    summarizeByGene(args.mappedNucleotides)
