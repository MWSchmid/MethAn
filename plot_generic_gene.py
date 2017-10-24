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
python %s stats figure

REQUIRED:
\tstats: a file holding chrom, pos, gene:
\t\tChr1\t5000\tAT1G01010,exon|AT1G01015,exon
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
\t-onesubplot <given> will produce a figure with a single subplot (instead of three)
\t\tdefault: false

""" % sys.argv[0]


if len(sys.argv) < 3:
	sys.exit(usage)

statsFile = sys.argv[1]
figFile = sys.argv[2]

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

if "-onesubplot" in sys.argv:
	oneSubPlot = True
else:
	oneSubPlot = False

if window not in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
	sys.exit(usage)

# a metagene class
class metaGene:
	
	def __init__(self, upstream, body, downstream, full):
		self.upstream = upstream
		self.body = body
		self.downstream = downstream
		self.full = full
	
	def __str__(self):
		a = list(self.upstream)
		b = list(self.body)
		c = list(self.downstream)
		d = list(self.full)
		out = '|'.join([str(x) for x in a]) + '\n' + '|'.join([str(x) for x in b]) + '\n' + '|'.join([str(x) for x in c]) + '\n' + '|'.join([str(x) for x in d]) + '\n'
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
		self.toCheck = set([])
		self.skipOV = skipOV
		self.skippedPositions = 0
		self.totalPositions = 0
		self.anno = {}
		self.loadAnnotation()
	
	def loadAnnotation(self):
		reg_meta_agi = re.compile('AT\w{2,2}\d{5,5}', re.IGNORECASE)
		with open("/home/marc/MethAn/data/TAIR10.v2.sorted.gff", 'r') as infile:
			for line in infile:
				if "protein_coding_gene" not in line:
					continue
				fields = line[:-1].split('\t')
				chrom = fields[0]
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
				(chrom, spos, agifield) = line[:-1].split('\t')
				agis = reg_meta_agi.findall(agifield)
				self.totalPositions += 1
				if self.skipOV and (len(agis) > 1):
					self.skippedPositions += 1
					continue
				self.cov[chrom][int(spos)] += 1
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
		rawUpstream = np.zeros(bordersize)
		rawDownstream = np.zeros(bordersize)
		for agi in list(self.toCheck):
			if agi not in self.anno.keys():
				continue
			(ls, le, bs, be, steplen, rs, re, reverse, use) = self.anno[agi].getStartEndStep(self.chromsizes, bordersize, fraction, endCorrection)
			if not use:
				continue #skip the ones that are too close to the chromosome border - also the ones that are shorter than the fraction
			# get a copy of up and downstream
			# NOTE - DONT do this - gives a skewed picture due to the overlap that only goes to 1 kb: 
			# smoothUpstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(ls-endCorrection):(le+endCorrection)], window_len, window)
			# upstreamCov = np.array([smoothUpstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
			rawUpstreamCov = np.array(self.cov[self.anno[agi].chrom][ls:le])
			smoothUpstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][ls:(le+endCorrection)], window_len, window)
			upstreamCov = np.array([smoothUpstreamCov[i] for i in xrange(0,bordersize,1)])
			# NOTE - DONT do this - gives a skewed picture due to the overlap that only goes to 1 kb: 
			# smoothDownstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(rs-endCorrection):(re+endCorrection)], window_len, window)
			# downstreamCov = np.array([smoothDownstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
			rawDownstreamCov = np.array(self.cov[self.anno[agi].chrom][rs:re])
			smoothDownstreamCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(rs-endCorrection):re], window_len, window)
			downstreamCov = np.array([smoothDownstreamCov[i] for i in xrange(endCorrection,(bordersize+endCorrection),1)])
			# process body
			smoothBodyCov = self.smooth_cov(self.cov[self.anno[agi].chrom][(le-endCorrection):(rs+endCorrection)], window_len, window)
			bodyCov = np.array([smoothBodyCov[i] for i in xrange((bs+endCorrection),(be+endCorrection),steplen)])
			# add it
			if reverse:
				upstream += downstreamCov[::-1]
				downstream += upstreamCov[::-1]
				body += bodyCov[::-1]
				rawUpstream += rawUpstreamCov[::-1]
				rawDownstream += rawDownstreamCov[::-1]
			else:
				upstream += upstreamCov
				downstream += downstreamCov
				body += bodyCov
				rawUpstream += rawUpstreamCov
				rawDownstream += rawDownstreamCov
		rawFull = np.concatenate((rawUpstream, body, rawDownstream))
		full = self.smooth_cov(rawFull, window_len, window)
		return metaGene(upstream, body, downstream, full)
	
# a plotter class
class metaGenePlotter:
	
	def __init__(self, commonGene, figFile, fraction, bordersize, usePNG):
		self.commonGene = commonGene
		self.fraction = fraction
		self.bordersize = bordersize
		self.figFile = figFile
		self.usePNG = usePNG
	
	def adjust_spines(self, ax, spines):#from matplotlib
		for loc, spine in ax.spines.iteritems():
			if loc in spines:
				spine.set_position(('outward', 10)) # outward by 10 points
			else:
				spine.set_color('none') # don't draw spine
		# turn off ticks where there is no spine
		if 'left' in spines:
			ax.yaxis.set_ticks_position('left')
		else:
			ax.yaxis.set_ticks([]) # no yaxis ticks
		if 'bottom' in spines:
			ax.xaxis.set_ticks_position('bottom')
		else:
			ax.xaxis.set_ticks([])  # no xaxis ticks

	def getPanel(self, ax, cov, ymax, wP):
		if wP == "body":
			xpos = np.arange(0, self.fraction, 1)
		elif wP == "down":
			xpos = np.arange(0, bordersize, 1)
		else:
			xpos = np.arange(-bordersize+1, 1, 1)
		ax.plot(xpos, cov, color = 'k', linestyle = '-', linewidth = 2)
		if wP == "body":
			self.adjust_spines(ax, ["bottom"])
			ax.set_xlabel("gene body (fraction)", size = "large")
			ax.set_xlim(0,self.fraction)
		elif wP == "down":
			self.adjust_spines(ax, ["bottom"])
			ax.set_xlabel("bp downstream", size = "large")
			ax.set_xlim(0,bordersize)
		else:
			self.adjust_spines(ax, ["left", "bottom"])
			ax.set_ylabel("count", size = "large")
			ax.set_xlabel("bp upstream", size = "large")
			ax.set_xlim(-bordersize,0)
		ax.set_ylim(0, ymax)
	
	def plotCommonGene(self):
		fig = plt.figure(figsize = (60, 10))
		# get the max value for y
		maxes = np.array([max(self.commonGene.upstream), max(self.commonGene.body), max(self.commonGene.downstream)])
		ymax = max(maxes)
		ax_up = fig.add_subplot(1,3,1)
		self.getPanel(ax_up, self.commonGene.upstream, ymax, "up")
		ax_body = fig.add_subplot(1,3,2)
		self.getPanel(ax_body, self.commonGene.body, ymax, "body")
		ax_down = fig.add_subplot(1,3,3)
		self.getPanel(ax_down, self.commonGene.downstream, ymax, "down")
		if self.usePNG:
			fig.savefig(self.figFile, format = "png", dpi = 300)
		else:
			fig.savefig(self.figFile, format = "svg", dpi = 300)
		plt.close("all")
	
	def plotCommonGeneSingle(self):
		fig = plt.figure(figsize = (30, 10))
		# get the max value for y
		ymax = max(self.commonGene.full)
		ax = fig.add_subplot(1,1,1)
		xpos = np.arange(0, (2*self.bordersize+self.fraction), 1)
		ax.plot(xpos, self.commonGene.full, color = 'k', linestyle = '-', linewidth = 2)
		# mark the breakpoints
		ax.axvline(self.bordersize)
		ax.axvline(self.bordersize+self.fraction)
		# arange some stuff
		self.adjust_spines(ax, ["left", "bottom"])
		ax.set_ylabel("count", size = "large")
		ax.set_xlabel("position", size = "large")
		ax.set_xlim(0, xpos.size)
		ax.set_ylim(0, ymax)
		if self.usePNG:
			fig.savefig(self.figFile, format = "png", dpi = 300)
		else:
			fig.savefig(self.figFile, format = "svg", dpi = 300)
		plt.close("all")


### runner
genome = genomeHandler(skipOV)
genome.loadPositions(statsFile)
commonGene = genome.getMetaGene(fraction, bordersize, window_len, window, endCorrection)
plotter = metaGenePlotter(commonGene, figFile, fraction, bordersize, usePNG)
if oneSubPlot:
	plotter.plotCommonGeneSingle()
else:
	plotter.plotCommonGene()
print "skipped %d ambiguous positions of %d positions (total)" % (genome.skippedPositions, genome.totalPositions)

