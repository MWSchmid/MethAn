import sys
import time

usage = """

usage:
python %s gff feature stats outfile [-useOverlap]

REQUIRED:
\tgff: a gff file with the annotation
\tfeature: the name of the feature in the gff file that shall be used
\tstats: a file holding chrom and pos:
\t\tChr1\t5000

""" % (sys.argv[0])

if len(sys.argv) < 5:
	sys.exit(usage)

gffFile = sys.argv[1]
feature = sys.argv[2]
statsFile = sys.argv[3]
resultsFile = sys.argv[4]
if "-useOverlap" in sys.argv:
	useOverlap = True
else:
	useOverlap = False

# a mini gene class
class miniGene:
	
	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.start = start
		self.end = end
	
	def queryPos(self, chrom, pos):
		dist = 100000
		if (chrom != self.chrom):
			#print >> sys.stderr, "chrombreak"
			#print >> sys.stderr, '\t'.join([chrom, str(pos), self.chrom, str(self.start), str(self.end)])
			status = "break"
		else:
			if (pos < self.start):
				dist = self.start - pos
				status = "break"
			elif (pos > self.end):
				dist = pos - self.end
				status = "continue"
			else:
				dist = 0
				status = "overlap"
				sDist = pos - self.start
				eDist = self.end - pos
				if sDist < eDist:
					dist = sDist
				else:
					dist = eDist
		return (status,dist)

# a class taking care of the annotation
class genomeHandler:
	
	def __init__(self, gffFile, feature):
		self.chromsizes = { '1': 30427671, '2': 19698289, '3': 23459830, '4': 18585056, '5': 26975502, 'Pt': 154478, 'Mt': 366924 }
		self.feature = feature
		self.gffFile = gffFile
		self.offsetstep = 1000
		self.anno = []
		self.loadAnnotation()
		self.index = {}
		self.loadIndex()
	
	def decode_gff(self, line):
		fields = line.split('\t')
		chrom = fields[0]
		source = fields[1]
		feature = fields[2]
		start = int(fields[3]) - 1 # to make it zero based
		end = int(fields[4]) -1 # to make it zero based
		try:
			score = float(fields[5])
		except ValueError:
			score = 0
		strand = fields[6]
		frame = fields[7]
		rest = fields[8]
		return (chrom, source, feature, start, end, score, strand, frame)
	
	def loadAnnotation(self):
		with open(self.gffFile, 'r') as infile:
			for line in infile:
				(chrom, source, feature, start, end, score, strand, frame) = self.decode_gff(line)
				if self.feature != feature:
					continue
				gene = miniGene(chrom, start, end)
				self.anno.append(gene)
	
	def getOffset(self, int_in):
		x = int_in/self.offsetstep
		return x
	
	def loadIndex(self):
		# similar to c++ common_functions -> create_index() - IMPORTANT THAT SORTED BY START (should be true) and: not via a pair (nested dict instead)
		for i in range(0, len(self.anno), 1):
			chrom = self.anno[i].chrom
			if not self.index.has_key(chrom): # this is just solved this way to make it a bit clearer down in the loops
				self.index.update({chrom:{0:i}})
			start = self.anno[i].start
			end = self.anno[i].end
			start_step = self.getOffset(start)
			end_step = self.getOffset(end)
			step_diff = end_step - start_step
			if step_diff == 0:
				if not self.index[chrom].has_key(start_step):
					self.index[chrom].update({start_step:i})
			else:
				for step_i in range(start_step, (end_step+1), 1): # the last step needs to be inside!
					if not self.index[chrom].has_key(step_i):
						self.index[chrom].update({step_i:i})
	
	def queryOverlap(self, chrom, pos):
		has_chrom_error = False
		has_overlap = False
		min_dist = 100000
		if not self.index.has_key(chrom):
			has_chrom_error = True
		else:
			start_step = self.getOffset(pos)
			search_start = self.index[chrom][0]
			while start_step > 0: # find the next lower offset that has a feature inside
				start_step -= 1
				try:
					search_start = self.index[chrom][start_step]
					break
				except KeyError:
					continue
			for i in xrange(search_start,len(self.anno),1):
				(status,dist) = self.anno[i].queryPos(chrom, pos)
				if dist < min_dist:
					min_dist = dist
				if status == "break":
					break
				elif status == "continue":
					continue
				elif status == "overlap":
					has_overlap = True
					min_dist = -1
					break
				else:
					sys.exit("unexpected ending of conditions")
		return (min_dist, has_overlap)

# load the anno
genome = genomeHandler(gffFile, feature)

# process the infile
#startTime = time.time()
#lineCounter = 0
withOverlap = 0
distances = [0]*100001
with open(statsFile, 'r') as infile:
	for line in infile:
		#lineCounter += 1
		#if (lineCounter % 500000) == 0:
		#	print >> sys.stderr, "time elapsed: %4.2f" % (time.time() - startTime)
		(chrom, strPos) = line[:-1].split('\t')
		pos = int(strPos)
		(minDist, hasOverlap) = genome.queryOverlap(chrom, pos)
		if hasOverlap:
			withOverlap += 1
			continue
		distances[minDist] += 1

with open(resultsFile, 'w') as outfile:
	if useOverlap:
		outfile.write('\t'.join(["#within", str(withOverlap)])+'\n')
	for i in xrange(0, len(distances), 1):
		outfile.write('\t'.join([str(i), str(distances[i])])+'\n')

if not useOverlap:
	print >> sys.stderr, "%d overlaps were found and skipped" % withOverlap