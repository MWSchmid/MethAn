#!/usr/bin/env python
import sys

usage = """

python %s mappingFile gf2dmc g2gf g2num

input:
mappingFile: The mapped nucleotides - output from MethAnMap.

output:
gf2dmc: gene-feature to number of cytosines
g2gf: gene to gene-feature
g2num: gene to number of cytosines

""" % sys.argv[0]

if len(sys.argv) < 5:
	sys.exit(usage)

infileName = sys.argv[1]
gf2dmcName = sys.argv[2]
g2gfName = sys.argv[3]
g2numName = sys.argv[4]

gf2dmc = {}
g2gf = {}
g2dmc = {}
with open(infileName, 'rb') as infile:
	for line in infile:
		fields = line[:-1].split()
		chrom, strand, pos, A, B, context, mapping = fields[0:7]
		dmc = '_'.join([chrom, pos])
		gfs = mapping.split('|')
		for gf in gfs:
			gene = gf.split(',')[0]
			if gene == "none":
				continue
			try:
				if gf not in g2gf[gene]:
					g2gf[gene].append(gf)
			except KeyError:
				g2gf[gene] = [gf]
			try:
				if dmc not in gf2dmc[gf]:
					gf2dmc[gf].append(dmc)
			except KeyError:
				gf2dmc[gf] = [dmc]
			try:
				if dmc not in g2dmc[gene]:
					g2dmc[gene].append(dmc)
			except KeyError:
				g2dmc[gene] = [dmc]

with open(gf2dmcName, 'wb') as outfile:
	for gf, dmc in gf2dmc.items():
		print >> outfile, gf + ';' + '|'.join(dmc)

with open(g2gfName, 'wb') as outfile:
	for gene, gf in g2gf.items():
		print >> outfile, gene + ';' + '|'.join(gf)

with open(g2numName, 'wb') as outfile:
	for gene, dmc in g2dmc.items():
		print >> outfile, gene + ';' + str(len(dmc))
