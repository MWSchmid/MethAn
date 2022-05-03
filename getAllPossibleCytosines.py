#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""getAllPossibleCytosines.py
This script can be used get chrom/pos/context for all possible cytosines.
            
Results are printed to std::out.
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import sys
import logging
import textwrap
import os
import gzip
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

class anySequence:
    """A class to handle sequences"""

    def __init__(self, name=None, seq=None):
        self.name = name
        self.seq = seq
        #self.clean()

    def clean(self):
        seq = ''
        for c in self.seq:
            if c in self.alphabet:
                seq += c
                self.seq = seq
            else:
                logging.warning(''.join(["unknown character in", self.name, "-", c]))

    def __str__(self):
        return '>'+self.name+'\n'+self.seq

    def __getitem__(self, i):
        return self.seq[i]

    def getname(self):
        return self.name

    def getseq(self):
        return self.seq

    def getlen(self):
        return len(self.seq)

class DNAsequence(anySequence):
    """A subclass to handle DNA sequences"""

    alphabet = "atcgrykmswbdhvnATCGRYKMSWBDHVN"

    def __init__(self, name=None, seq=None):
        anySequence.__init__(self, name, seq.upper())

    def gc(self):
        """GC percent"""
        count_c = self.seq.count('G')
        count_g = self.seq.count('G')
        return float(count_c + count_g) / len(self.seq)
    
    def revcompl(self):
        """reverse complement"""
        trans = {"A":"T",
                 "T":"A",
                 "C":"G",
                 "G":"C",
                 "R":"Y",
                 "Y":"R",
                 "K":"M",
                 "M":"K",
                 "S":"W",
                 "W":"S",
                 "B":"V",
                 "D":"H",
                 "H":"D",
                 "V":"B",
                 "N":"N"}
        revcompseq = ''.join([trans[x] for x in self.seq[::-1]])
        return revcompseq

    def translate(self, frame=0, mergeFrameToName=False):
        """translation into a protein. frame: 0, 1, 2, -1, -2, -3"""
        Standard_Genetic_Code = {'TTT':'F','TCT':'S','TAT':'Y','TGT':'C','TTC':'F','TCC':'S','TAC':'Y','TGC':'C','TTA':'L','TCA':'S','TAA':'-','TGA':'-','TTG':'L','TCG':'S','TAG':'-','TGG':'W','CTT':'L','CCT':'P','CAT':'H','CGT':'R','CTC':'L','CCC':'P','CAC':'H','CGC':'R','CTA':'L','CCA':'P','CAA':'Q','CGA':'R','CTG':'L','CCG':'P','CAG':'Q','CGG':'R','ATT':'I','ACT':'T','AAT':'N','AGT':'S','ATC':'I','ACC':'T','AAC':'N','AGC':'S','ATA':'I','ACA':'T','AAA':'K','AGA':'R','ATG':'M','ACG':'T','AAG':'K','AGG':'R','GTT':'V','GCT':'A','GAT':'D','GGT':'G','GTC':'V','GCC':'A','GAC':'D','GGC':'G','GTA':'V','GCA':'A','GAA':'E','GGA':'G','GTG':'V','GCG':'A','GAG':'E','GGG':'G',"GTN":"V", "TCN":"S", "CCN":"P", "ACN":"T", "GCN":"A", "GGN":"G","CGN":"R"}
        originalframe = frame
        if frame < 0 :
            seq = self.revcompl()
            frame = abs(frame) - 1
        else:
            seq = self.seq
        if frame > 2:
            return ''
        protseq = []
        for i in range(frame,len(seq) - 2,3):
            codon = seq[i:i+3]
            try:
                protseq.append(Standard_Genetic_Code[codon])
            except KeyError:
                protseq.append("X")
        protseq = ''.join(protseq)
        if mergeFrameToName:
            newName = self.name + "_" + str(originalframe)
        else:
            newName = self.name + " | translated_frame=" + str(originalframe)
        return PROTEINsequence(name = newName, seq = protseq)

    def setname(self, name):
        self.name = name

class PROTEINsequence(anySequence):
    '''a subclass for proteins'''

    alphabet = "-acedgfihkmlnqpsrtwvyxACEDGFIHKMLNQPSRTWVYX"

    def __init__(self, name=None, seq=None):
        anySequence.__init__(self, name, seq.upper())

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

def loadSequenceFromFasta(infileName):
    out = {}
    with myopen(infileName, 'rb') as infile:
        for line in infile:
            line = line.decode("ascii").rstrip('\n')
            if line[0] == ">":
                seqName = line[1:].strip().split(' ')[0]
                out[seqName] = []
            else:
                out[seqName].append(line)
    for seqName, seqList in out.items():
        out[seqName] = DNAsequence(seqName, ''.join(seqList))
    return out

def printChromPosContext(chrom, seq):
    seq = seq.upper()
    for i, b in enumerate(seq):
        if b != 'C':
            continue
        try:
            subSeq = seq[i:(i+3)]
        except IndexError:
            continue
        numOKchars = sum([subSeq.count(x) for x in "ACGTH"])
        if numOKchars != 3:
            continue
        if seq[i+1] == 'G': # CG
            sys.stdout.write('\t'.join([chrom, str(i), "CG"])+'\n')
            continue
        if seq[i+2] == 'G': # CHG
            sys.stdout.write('\t'.join([chrom, str(i), "CHG"])+'\n')
            continue
        sys.stdout.write('\t'.join([chrom, str(i), "CHH"])+'\n') #CHH

def getAllPossibleCytosines(fastaFile):
    sequences = loadSequenceFromFasta(fastaFile)
    for seqName, seq in sequences.items():
        printChromPosContext(seqName, seq.getseq())
        printChromPosContext(seqName, seq.revcompl())

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="getAllPossibleCytosines.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            This script can be used get chrom/pos/context for
            all possible cytosines.

            Results are printed to std::out.
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("fastaFile", nargs='+', type=str,
                        help="""
                        One or more fasta files.
                        """)
    
    args = parser.parse_args()
    
    for fastaFile in args.fastaFile:
        fastaFileName = os.path.basename(fastaFile)
        getAllPossibleCytosines(fastaFile)
    
    
        

