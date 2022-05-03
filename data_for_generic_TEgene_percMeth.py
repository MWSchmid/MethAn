#!/usr/bin/env python
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re

usage = """

usage:
python %s stats outFile

WARNING:

only works if saturated! OR USE t-onesubplot TO FIX IT, then it will first get the percentages taking the numbers into account and smooth at the very end

REQUIRED:
\tstats: a file holding chrom, pos, gene:
\t\tChr1\t5000\tAT1G01010,exon|AT1G01015,exon\tpercMeth
\tfigure: the file where the plot shall be saved

OPTIONAL:
\t-fraction <int> number of bins for gene body
\t\tdefault: 100
\t-bordersize <int> size of the flanking regions
\t\tdefault: 1000
\t-window_len <int> size of the smoothing window
\t\tdefault: 50 
\t-window <str> type of smoothing window [flat, hanning, hamming, bartlett, blackman]
\t\tdefault: blackman
\t-png <given> will create png instead of svg
\t\tdefault: false
\t-endcorrection <int> to avoid strange behaviour at the ends, this number of bases will be appended for the smoothing
\t\tdefault: 100
\t-skip <given> will omit ambiguous positions
\t\tdefault: false

""" % sys.argv[0]


if len(sys.argv) < 3:
    sys.exit(usage)

statsFile = sys.argv[1]
outFile = sys.argv[2]

if "-fraction" in sys.argv:
    fraction = int(sys.argv[sys.argv.index("-fraction") + 1])
else:
    fraction = 100

if "-bordersize" in sys.argv:
    bordersize= int(sys.argv[sys.argv.index("-bordersize") + 1])
else:
    bordersize = 1000

if "-window_len" in sys.argv:
    window_len = int(sys.argv[sys.argv.index("-window_len") + 1])
else:
    window_len = 50

if "-window" in sys.argv:
    window = sys.argv[sys.argv.index("-window") + 1]
else:
    window = "blackman"

if "-png" in sys.argv:
    usePNG = True
else:
    usePNG = False

if "-endcorrection" in sys.argv:
    endCorrection = int(sys.argv[sys.argv.index("-endcorrection") + 1])
else:
    endCorrection = 100

if "-skip" in sys.argv:
    skipOV = True
else:
    skipOV = False

if window not in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
    sys.exit(usage)

# a metagene class
class metaGene:
    
    def __init__(self, upstream, body, downstream, full, rawFull):
        self.upstream = upstream
        self.body = body
        self.downstream = downstream
        self.full = full
        self.rawFull = rawFull
    
    def __str__(self):
        a = list(self.upstream)
        b = list(self.body)
        c = list(self.downstream)
        d = list(self.full)
        e = list(self.rawFull)
        out = '|'.join([str(x) for x in a]) + '\n' + '|'.join([str(x) for x in b]) + '\n' + '|'.join([str(x) for x in c]) + '\n' + '|'.join([str(x) for x in d]) + '\n' + '|'.join([str(x) for x in e]) + '\n'
        out = ','.join([str(x) for x in e])
        return (out)

# a mini gene class
class miniGene:
    
    def __init__(self, chrom, strand, start, end):
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
    
    def getStartEndStep(self, chromsizes, bordersize, fraction, endCorrection):
        '''returns start, end, stepping and reverse for: upstream (5prime), body and downstream (3prime). This of course depends on the strand'''
        use = True
        # the left sie
        ls = self.start - bordersize
        le = self.start
        if (ls-endCorrection) < 0:
            ls = 0
            use = False
        # the right side
        rs = self.end
        re = self.end + bordersize
        if (re+endCorrection) > chromsizes[self.chrom]:
            re = chromsizes[self.chrom]
            use = False
        # the body
        steplen = int(self.end-self.start)/int(fraction) # int() just for clarity - must be integer division
        bs = int(steplen)/int(2) # if one would always start at 0, one of the ends would be preferentially treated
        be = fraction*steplen
        if steplen == 0:
            use = False
        # return depending on strand
        if self.strand == '+':
            reverse = False
        elif self.strand == '-':
            reverse = True
        else:
            print >> sys.stderr, "unexpected strand"
        return ( ls, le, bs, be, steplen, rs, re, reverse , use )
    


# a genome handler class
class genomeHandler:
    
    def __init__(self, skipOV):
        self.chromsizes = { '1': 30427671, '2': 19698289, '3': 23459830, '4': 18585056, '5': 26975502, 'Pt': 154478, 'Mt': 366924}
        self.cov = dict([(chrom,np.zeros(self.chromsizes[chrom])) for chrom in self.chromsizes])
        self.num = dict([(chrom,np.zeros(self.chromsizes[chrom])) for chrom in self.chromsizes])
        self.toCheck = set([])
        self.skipOV = skipOV
        self.skippedPositions = 0
        self.totalPositions = 0
        self.anno = {}
        self.loadAnnotation()
    
    def loadAnnotation(self):
        reg_meta_agi = re.compile('AT\w{2,2}\d{5,5}', re.IGNORECASE)
        chromTranslator = {"Chr1":"1", "Chr2":"2", "Chr3":"3", "Chr4":"4", "Chr5":"5", "ChrC":"Mt", "ChrM":"Pt"}
        with open("/media/mwschmid/archiveNoBackup/genomes/arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff", 'r') as infile:
            for line in infile:
                if "transposable_element_gene" not in line:
                    continue
                fields = line[:-1].split('\t')
                chrom = chromTranslator[fields[0]]
                start = int(fields[3]) - 1 # to make it zero-based
                end = int(fields[4]) # - 1 # to make it zero-based -> not applicable here because we access it via the list
                strand = fields[6]
                agi = reg_meta_agi.findall(fields[8])[0]
                gene = miniGene(chrom, strand, start, end)
                self.anno.update({agi:gene})
    
    def loadPositions(self, statsFile):
        reg_meta_agi = re.compile('AT\w{2,2}\d{5,5}', re.IGNORECASE)
        with open(statsFile, 'r') as infile:
            for line in infile:
                (chrom, spos, agifield, percMeth) = line[:-1].split('\t')
                agis = reg_meta_agi.findall(agifield)
                self.totalPositions += 1
                if self.skipOV and (len(agis) > 1):
                    self.skippedPositions += 1
                    continue
                self.cov[chrom][int(spos)] += float(percMeth)
                self.num[chrom][int(spos)] += 1
                for agi in agis: # version where overlapping genes are taken NOTE given the priorities, overlaps are only visible if both genes were hit in the same spot -> one would have to run it with equal priorities
                    self.toCheck.add(agi)
    
    def smooth_cov(self, cov_array , window_len = 50, window = 'blackman'): 
        # taken from http://www.scipy.org/Cookbook/SignalSmooth -> more details there. note: window_len should be an odd integer. window = flat is moving average
        x = cov_array
        if x.size < window_len:
            window_len = x.size
            #raise ValueError, "Input vector needs to be bigger than window size."
        if x.ndim != 1:
            raise ValueError, "smooth only accepts 1 dimension arrays."
        if window_len < 3:
            return list(x)
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
        if window == 'flat': 
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')
        y = np.convolve(w / w.sum(), s, mode = 'same')
        return y[window_len - 1:-window_len + 1]
    
    def getMetaGene(self, fraction, bordersize, window_len, window, endCorrection):
        upstream = np.zeros(bordersize)
        body = np.zeros(fraction)
        downstream = np.zeros(bordersize)
        rawUpstream = np.zeros(fraction)
        rawUpstreamNum = np.zeros(fraction)
        rawBody = np.zeros(fraction)
        rawBodyNum = np.zeros(fraction)
        rawDownstream = np.zeros(fraction)
        rawDownstreamNum = np.zeros(fraction)
        divideBy = 0
        for agi in list(self.toCheck):
            if agi not in self.anno.keys():
                continue
            (ls, le, bs, be, steplen, rs, re, reverse, use) = self.anno[agi].getStartEndStep(self.chromsizes, bordersize, fraction, endCorrection)
            if not use:
                continue #skip the ones that are too close to the chromosome border - also the ones that are shorter than the fraction
            divideBy += 1
            # get a copy of up and downstream
            # NOTE - DONT do this - gives a skewed picture due to the overlap that only goes to 1 kb: 
            # smoothUpstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(ls-endCorrection):(le+endCorrection)], window_len, window)
            # upstreamCov = np.array([smoothUpstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
            #rawUpstreamCov = np.array(self.cov[self.anno[agi].chrom][ls:le])
            #rawUpstreamCovNum = np.array(self.num[self.anno[agi].chrom][ls:le])
            steplenForRaw = (le-ls)/float(fraction)
            rawUpstreamCov = np.array([sum(self.cov[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(ls, le, steplenForRaw)])
            rawUpstreamCovNum = np.array([sum(self.num[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(ls, le, steplenForRaw)])
            rawUpstreamCov = rawUpstreamCov[0:fraction] # sometimes it's not 100
            rawUpstreamCovNum = rawUpstreamCovNum[0:fraction] # sometimes it's not 100
            smoothUpstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][ls:(le+endCorrection)], window_len, window)
            upstreamCov = np.array([smoothUpstreamCov[i] for i in xrange(0,bordersize,1)])
            # NOTE - DONT do this - gives a skewed picture due to the overlap that only goes to 1 kb: 
            # smoothDownstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(rs-endCorrection):(re+endCorrection)], window_len, window)
            # downstreamCov = np.array([smoothDownstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
            #rawDownstreamCov = np.array(self.cov[self.anno[agi].chrom][rs:re])
            #rawDownstreamCovNum = np.array(self.num[self.anno[agi].chrom][rs:re])
            steplenForRaw = (re-rs)/float(fraction)
            rawDownstreamCov = np.array([sum(self.cov[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(rs, re, steplenForRaw)])
            rawDownstreamCovNum = np.array([sum(self.num[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(rs, re, steplenForRaw)])
            rawDownstreamCov = rawDownstreamCov[0:fraction] # sometimes it's not 100
            rawDownstreamCovNum = rawDownstreamCovNum[0:fraction] # sometimes it's not 100
            smoothDownstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(rs-endCorrection):re], window_len, window)
            downstreamCov = np.array([smoothDownstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
            # process body
            smoothBodyCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(le-endCorrection):(rs+endCorrection)], window_len, window)
            bodyCov = np.array([smoothBodyCov[i] for i in xrange((bs+endCorrection),(be+endCorrection),steplen)])
            steplenForRaw = (rs-le)/float(fraction)
            rawBodyCov = np.array([sum(self.cov[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(le, rs, steplenForRaw)])
            rawBodyCovNum = np.array([sum(self.num[self.anno[agi].chrom][int(i):int(i+steplenForRaw)]) for i in np.arange(le, rs, steplenForRaw)])
            rawBodyCov = rawBodyCov[0:fraction] # sometimes it's not 100
            rawBodyCovNum = rawBodyCovNum[0:fraction] # sometimes it's not 100
            #print >> sys.stderr, rawBodyCov
            if reverse:
                upstream += downstreamCov[::-1]
                downstream += upstreamCov[::-1]
                body += bodyCov[::-1]
                temp = rawUpstreamCov[::-1]
                tempNum = rawUpstreamCovNum[::-1]
                for i in xrange(0, len(temp)):
                    if tempNum[i] > 0:
                        rawUpstream[i] += temp[i]
                        rawUpstreamNum[i] += tempNum[i]
                temp = rawDownstreamCov[::-1]
                tempNum = rawDownstreamCovNum[::-1]
                for i in xrange(0, len(temp)):
                    if tempNum[i] > 0:
                        rawDownstream[i] += temp[i]
                        rawDownstreamNum[i] += tempNum[i]
                temp = rawBodyCov[::-1]
                tempNum = rawBodyCovNum[::-1]
                for i in xrange(0, len(temp)):
                    if tempNum[i] > 0:
                        rawBody[i] += temp[i]
                        rawBodyNum[i] += tempNum[i]
            else:
                upstream += upstreamCov
                downstream += downstreamCov
                body += bodyCov
                for i in xrange(0, len(rawUpstreamCov)):
                    if rawUpstreamCovNum[i] > 0:
                        rawUpstream[i] += rawUpstreamCov[i]
                        rawUpstreamNum[i] += rawUpstreamCovNum[i]
                for i in xrange(0, len(rawDownstreamCov)):
                    if rawDownstreamCovNum[i] > 0:
                        rawDownstream[i] += rawDownstreamCov[i]
                        rawDownstreamNum[i] += rawDownstreamCovNum[i]
                for i in xrange(0, len(rawBodyCov)):
                    if rawBodyCovNum[i] > 0:
                        rawBody[i] += rawBodyCov[i]
                        rawBodyNum[i] += rawBodyCovNum[i]
            #print >> sys.stderr, rawBody
        rawUpstream /= rawUpstreamNum
        rawBody /= rawBodyNum
        rawDownstream /= rawDownstreamNum
        rawFull = np.concatenate((rawUpstream, rawBody, rawDownstream))
        full = self.smooth_cov(rawFull, window_len, window)
        upstream /= divideBy
        body /= divideBy
        downstream /= divideBy
        return metaGene(upstream, body, downstream, full, rawFull)

### runner
genome = genomeHandler(skipOV)
genome.loadPositions(statsFile)
commonGene = genome.getMetaGene(fraction, bordersize, window_len, window, endCorrection)
with open(outFile, "w") as outfile:
    print >> outfile, str(commonGene)

