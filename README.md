# MethAn
A collection of methylome analysis tools

Note:
* the tools are provided "as-they-are"
* I used them for some methylome analysis in Arabidopsis - thus, the Python and R scripts:
    * contain some hard-coded parameters (e.g. chromosome names are 1, 2, 3, 4, 5, Mt and Pt)
    * are not optimized for large genomes (plot_generic_gene.py creates for each chromosome a list as long as the chromosome)
* MethAnMap, the binary for mapping positions to the annotation is an exception in this respect. It deals with any genome (as long as the annotation and the genome sequence match).
* MethAnPre is as well generalized.

## MethAnPre

The scripts in [MethAnPre](MethAnPre) can be used for trimming, alignment, duplicate removal and methylation extraction. The scripts require:

* [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [picard](https://broadinstitute.github.io/picard)
* [Pile-O-Meth](https://bioconda.github.io/recipes/pileometh/README.html)
* [samtools](https://github.com/samtools/samtools)
* [tabix and bgzip](http://www.htslib.org/doc/tabix.html)

### Install requirements for scripts in MethAnPre

Download the [requirements](MethAnPre/requirements.txt?raw=true) to install some of the dependencies with conda.

```SH
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install unzip build-essential zlibc zlib1g zlib1g-dev tabix git

# install miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# add conda to path
export PATH="/home/ubuntu/miniconda3/bin:$PATH"

# logout/login again
conda config --add channels r
conda config --add channels bioconda

# install whatever you can with conda
conda install --file requirements.txt

# bismark 0.16.3
wget http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/bismark_v0.16.3.tar.gz
tar xzf bismark_v0.16.3.tar.gz
find ~/bismark_v0.16.3 -maxdepth 1 -type f -perm /a+x -exec sudo cp {} /usr/local/bin \;

# trim_galore 0.4.1
wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
unzip trim_galore_v0.4.1.zip
sudo cp trim_galore_zip/trim_galore /usr/local/bin/

# trimmomatic 0.36
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip

# PileOMeth
git clone https://github.com/dpryan79/PileOMeth.git
cd PileOMeth
make
sudo cp PileOMeth /usr/local/bin/
cd
```

**TODO** PileOMeth is now MethylDackel - TODO - check and update this README

**NOTE/TODO** We recently observed quite some differences with trimmomatic and trim_galore - check why?

### Download and index a genome, trim, align, de-duplicate, extract methylation

```SH
# some paths and variables
GENOME_FOLDER="/path/to/the/folder/with/the/genome"
GENOME_FASTA="At.fasta"
INDIR="/path/to/the/folder/containing/the/reads"
OUTDIR="/path/to/the/folder/where/alignments/and/counts/will/be/stored"

# add conda bin to path
export PATH="/home/ubuntu/miniconda3/bin:$PATH"

# download a genome and index it (TAIR10 as an example)
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-30/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.30.dna.genome.fa.gz
mkdir -p $GENOME_FOLDER
gunzip -c Arabidopsis_thaliana.TAIR10.30.dna.genome.fa.gz > $GENOME_FOLDER/$GENOME_FASTA
cd $GENOME_FOLDER
bismark_genome_preparation --bowtie2 --genomic_composition ./
cd

# array of samples (i.e. file prefixes, _R1/2.fq.gz will be added later on)
MYSAMPLES=("leaf_1" "leaf_2" "leaf_3" "root_1" "root_2" "root_3")

### Trimming with trim_galore

# SINGLE-END reads
for PREFIX in "${MYSAMPLES[@]}"; do
FASTQFILE="${PREFIX}_R1.fq.gz"
bismark_SE.sh $INDIR $OUTDIR $PREFIX $FASTQFILE $GENOME_FOLDER
bamToBedGraph.sh $OUTDIR $OUTDIR $PREFIX "${PREFIX}.dupMarked.sorted.bam" "$GENOME_FOLDER/$GENOME_FASTA"
done

# PAIRED-END reads
for PREFIX in "${MYSAMPLES[@]}"; do
FASTQFILE="${PREFIX}_R1.fq.gz"
FASTQFILEREVERSE="${PREFIX}_R2.fq.gz"
bismark_PE.sh $INDIR $OUTDIR $PREFIX $FASTQFILE $FASTQFILEREVERSE $GENOME_FOLDER
bamToBedGraph.sh $OUTDIR $OUTDIR $PREFIX "${PREFIX}.dupMarked.sorted.bam" "$GENOME_FOLDER/$GENOME_FASTA"
done

### Trimming with trimmomatic

# SINGLE-END reads
for PREFIX in "${MYSAMPLES[@]}"; do
FASTQFILE="${PREFIX}_R1.fq.gz"
bismark_SE_trimmomatic.sh $INDIR $OUTDIR $PREFIX $FASTQFILE $GENOME_FOLDER
bamToBedGraph.sh $OUTDIR $OUTDIR $PREFIX "${PREFIX}.dupMarked.sorted.bam" "$GENOME_FOLDER/$GENOME_FASTA"
done

# PAIRED-END reads
for PREFIX in "${MYSAMPLES[@]}"; do
FASTQFILE="${PREFIX}_R1.fq.gz"
FASTQFILEREVERSE="${PREFIX}_R2.fq.gz"
bismark_PE_trimmomatic.sh $INDIR $OUTDIR $PREFIX $FASTQFILE $FASTQFILEREVERSE $GENOME_FOLDER
bamToBedGraph.sh $OUTDIR $OUTDIR $PREFIX "${PREFIX}.dupMarked.sorted.bam" "$GENOME_FOLDER/$GENOME_FASTA"
done
```

### Merge BedGraphs and annotation into a long-format table

**TODO**

### Run a linear model for each cytosine

**TODO**

## MethAnMap
MethAnMap annotates cytosine positions to the annotation. A binary built on Kubuntu 16.04 is in the repository or can be downloaded [here](MethAnMap/MethAnMap?raw=true). Below some examples.

Required input:
* A table with nucleotide positions with the follwing first five columns (no header):
    * chromosome (same id as in annotation)
    * position (zero-based)
    * context (CG/CHG/CHH)
    * coverage
    * percent methylation
* The annotation in Rcount-XML format:
    * can be generated with [Rcount-format](https://github.com/MWSchmid/Rcount).
    * there is an annotation with 2 kb flanking regions [available here](data/Araport11_annotation_with_2kb_flank.xml.zip?raw=true). It was generated with a [reformatted Araport 11 GFF](data/Araport11_ref_noChr.gff.zip?raw=true) and is compatible with the TAIR10 reference assembly (just make sure that the chromosomes are named 1, 2, 3, 4, 5, Mt and Pt). The priorities for the annotation were set like shown [here](data/Araport11_annotation_with_2kb_flank.png). Finally, there is also a [TEGtoTEfile available here](data/TEGtoTE.txt.zip?raw=true) (see below for what it can be used).

    
Note that the mapping statistics will make use of the priorities set with Rcount-format.

```SH
### Required files ###

# Input
myNucleotides="/path/to/a/file/with/nucleotides.txt"
myAnnotation="/path/to/an/annotation.xml"

# Output (the myRegionFile is optional)
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.mapped"
myRegions="/path/to/a/file/with/the/regions.txt"
myMapStats="/path/to/a/file/with/the/mapping/nucleotides.stats"

### Examples ###

# Map all cytosines in the file
MethAnMap -S ${myNucleotides} -A ${myAnnotation} -O ${myMappedNucleotides} > ${myMapStats}

# Map all cytosines with a coverage of at least 10
MethAnMap -m 10 -S ${myNucleotides} -A ${myAnnotation} -O ${myMappedNucleotides} > ${myMapStats}

# Map only a random subset of all cytosines (e.g. 50/20/20 % of the CG/CHG/CHH nucleotides)
MethAnMap -Z 0.5,0.2,0.2 -S ${myNucleotides} -A ${myAnnotation} -O ${myMappedNucleotides} > ${myMapStats}

# Note regarding the mapping:
# If you are only interested in the mapping statistics, you can set -O to "SKIP"
MethAnMap -S ${myNucleotides} -A ${myAnnotation} -O SKIP > ${myMapStats}
MethAnMap -m 10 -S ${myNucleotides} -A ${myAnnotation} -O SKIP > ${myMapStats}
MethAnMap -Z 0.5,0.2,0.2 -S ${myNucleotides} -A ${myAnnotation} -O SKIP > ${myMapStats}

# Map nucleotides and get regions (with at least 3 cytosines and a maximal
# distance of 250 bp between two neighboring cytosines)
MethAnMap -S ${myNucleotides} -A ${myAnnotation} -O ${myMappedNucleotides} -R ${myRegions} > ${myMapStats}

# Get only the regions [works only for random subsets - i.e. if -Z is set]
MethAnMap -Z 0.5,0.2,0.2 -S ${myNucleotides} -A SKIP -O SKIP -R ${myRegions} > ${myMapStats}
```

In case you would like to build MethAnMap:

```SH
# Get the MethAn and Rcount repositories
git clone https://github.com/MWSchmid/Rcount
git clone https://github.com/MWSchmid/MethAn

# Open MethAnMap.pro with QtCreator and build the program
```

## Extract some useful files for the analysis in R and plot meta-genes

```SH
# ===============================================================================
# Create mapping files for R (optional but faster than R-standalone)
# Input
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.txt"

# Output
gf2dmc="/path/to/file/with/gene-feature/to/number/of/cytosines.txt"
g2gf="/path/to/file/with/gene/to/gene-feature.txt"
g2num="/path/to/file/with/gene/to/number/of/cytosines.txt"

# Run makeMappingFiles.py
python makeMappingFiles.py ${myMappedNucleotides} ${gf2dmc} ${g2gf} ${g2num}

# ===============================================================================
# Plot metagenes
# Important - in plot_generic_gene.py, you need to adjust (class genomeHandler): 
#   chromosome name and sizes (chromsizes)
#   the path to the GFF file (in the function loadAnnotation)
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.txt"
workingDirectory="/path/to/a/working/directory"
borderSize=1000
geneFraction=100
windowType=blackman
windowSize=100
endCorrection=100

cd ${workingDirectory}
awk -v OFS='\t' '{print $1,$3,$7 > "forMetagene_ALL.txt"; print $1,$3,$7 > "forMetagene_"$6".txt"}' \
    ${myMappedNucleotides}
cd /path/to/MethAn

for CONTEXT in ALL CG CHG CHH; do
inputFile="${workingDirectory}/forMetagene_${CONTEXT}.txt"
outputFile="${workingDirectory}/metagene_${CONTEXT}.svg"
python plot_generic_gene.py inputFile outputFile -fraction ${geneFraction} -window ${windowType} \
    -window_len ${windowSize} -bordersize ${borderSize} -endcorrection ${endCorrection}
done

# ===============================================================================
# Get distances to a specific feature of the annotation
# Important - in get_distance.py, you need to adjust (class genomeHandler): 
#   chromosome name and sizes (chromsizes)
# Note:
# you can use -useOverlap to include Cs which are within a feature. The number of
# overlapping Cs will be in the first row (starts with hash-tag to skip reading it later on)
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.txt"
workingDirectory="/path/to/a/working/directory"

cd ${workingDirectory}
awk -v OFS='\t' '{print $1,$3 > "forMetagene_ALL.txt"; print $1,$3 > "forDistance_"$6".txt"}' \
    ${myMappedNucleotides}
cd /path/to/MethAn

gffFile="/path/to/a/gff/file.gff"
gffFeature="one-GFF-feature"
for CONTEXT in ALL CG CHG CHH; do
inputFile="${workingDirectory}/forDistance_${CONTEXT}.txt"
outputFile="${workingDirectory}/distance_to_${gffFeature}_${CONTEXT}.txt"
python get_distance.py ${gffFile} ${gffFeature} ${inputFile} ${outputFile} -useOverlap
done

# ===============================================================================
# Get distances to a specific feature of the annotation for random sets of Cs
# We need to first generate some random sets - given the following input:
mySigNucleotides="/path/to/a/file/with/SIGNIFICANT_NUCLEOTIDES.txt"
myNucleotides="/path/to/a/file/with/ALL_POSSIBLE_NUCLEOTIDES.txt"
myAnnotation="/path/to/an/annotation.xml" # like above
numSets="the number of sets"
randomDir="/path/to/a/folder/where/random/sets/will/be/stored"
outSuffix="a suffix without file extension"

createRandomNucleotideSets.sh -d ${myNucleotides} ${mySigNucleotides} ${myAnnotation} \
    ${numSets} ${randomDir} ${outSuffix}

# and now the distances to the feature
gffFile="/path/to/a/gff/file.gff"
gffFeature="one-GFF-feature"
for RSN in $(seq 1 "${numberOfSets}"); do
RS="random_${RSN}"
echo "processing ${RS}..."
for CONTEXT in ALL CG CHG CHH; do
inputFile="${randomDir}/${RS}${outSuffix}_forDistance_${CONTEXT}.txt"
outputFile="${randomDir}/${RS}${outSuffix}_${gffFeature}_${CONTEXT}.txt"
python get_distance.py ${gffFile} ${gffFeature} ${inputFile} ${outputFile} -useOverlap
done
done

# ===============================================================================
# Note - the distance files should be named like this:
# observed data: <aPrefix>_<feature>_<context>[_<gain/loss>].txt
# random data: random_<number>_<aPrefix>_<feature>_<context>[_<gain/loss>].txt
# [_<gain/loss>] is optional
```

## R functions

```{R}
library("gplots")
library("graphics")
library("RColorBrewer")
source("/path/to/MethAn/MethAnR.R")
dataDir <- "/path/to/directory/with/nucleotide/tables"
randomDir <- "/path/to/directory/with/random/distance/files"
rDir <- "/path/to/directory/with/the/results"

######################################################################################
## data reading
######################################################################################

# a table with chrom and pos for all possible nucleotides
# (see ALL_POSSIBLE_NUCLEOTIDES.txt above)
myDataFullTable <- read.table(file.path(dataDir, "allPossiblePositions.txt"), sep='\t',
                              stringsAsFactors=F, col.names=c("chrom", "pos"))

# the table with the candidate positions (mapped)
# define columns here - this ones are per default in the file - if the original 
# SIGNIFICANT_NUCLEOTIDES.txt had more columns, you can add them here at the end
myDataCols <- c("chrom", "strand", "pos", "meth", "unmeth", "context", "mapping")
myData <- read.table(file.path(dataDir, "significantPositions.mapped"), sep='\t', header=F,
                     stringsAsFactors=F, row.names=NULL, col.names=myDataCols)
rownames(myData) <- paste0(myData$chrom, '_', myData$pos)

######################################################################################
## test different region parameters
######################################################################################

res <- f.analyze.regions(myData, rDir, "aPrefix", seq(5, 100, by = 5),
                         c(10,20,25,30,40,50,75,100,200,300), useLog = TRUE)

######################################################################################
## creat regions with a given set of parameters
######################################################################################

myRegs <- f.create.regions(myData, minNum = 50, maxDist = 100)

# if you want to extract all the cytosines which are within the regions
myCsInRegs <- f.extract.Cs.within.regions(myData, myRegs)
myDataSub <- myData[myCsInRegs,]

######################################################################################
## plot cytosine and region density along chromosomes
######################################################################################

f.plot.pos.density(myData, rDir)
f.plot.region.density(myRegs, rDir)


######################################################################################
## draw a histogram with the distances between significant and random Cs
######################################################################################

f.draw.distance.between.DMCs.combined(myData, myDataFullTable, rDir, 10, "combined_")

#####################################################################################
## draw distances to siRNA and transposons (for this you need to get the mappings first)
# NOTE THAT THIS DOES NOT TAKE THE CYTOSINES WITHIN THE FEATURES INTO ACCOUNT
######################################################################################

infilePrefix <- "see outSuffix above and info on how to name distance files"

f.draw.distance.to.something.combined(dataDir, rDir, "si_gene", TRUE, "smooth_",
                                      infilePrefix, sum, 10, randomDir)
f.draw.distance.to.something.combined(dataDir, rDir, "transposable_element", TRUE, "smooth_",
                                      infilePrefix, sum, 10, randomDir)

f.draw.distance.to.something.boxplot(dataDir, rDir, "si_gene", "boxplot_",
                                     infilePrefix, randomDir)
f.draw.distance.to.something.boxplot(dataDir, rDir, "transposable_element", "boxplot_",
                                     infilePrefix, randomDir)

######################################################################################
## test distances to siRNA and transposons (t-test or empirical)
# NOTE THAT THIS DOES NOT TAKE THE CYTOSINES WITHIN THE FEATURES INTO ACCOUNT
######################################################################################

testRes <- f.test.distance.to.something(rDir, "si_gene", infilePrefix, 10, randomDir)
testRes <- f.test.distance.to.something(rDir, "transposable_element", infilePrefix, 10, randomDir)

# the empirical test requires a lot of random sets...
testResEmpSiRNA <- f.test.distance.to.something.empirical(rDir, "si_gene", 500, randomDir)
testResEmpTE <- f.test.distance.to.something.empirical(rDir, "transposable_element", 500, randomDir)
write.table(testResEmpSiRNA, file.path(rDir, "distance_test_siRNA.txt"), sep='\t', quote=F)
write.table(testResEmpTE, file.path(rDir, "distance_test_TEs.txt"), sep='\t', quote=F)

######################################################################################
## plot number of cytosines per gene
######################################################################################

gf2dmcFile <- "/path/to/file/with/gene-feature/to/number/of/cytosines.txt"
g2gfFile <- "/path/to/file/with/gene/to/gene-feature.txt"
g2numFile <- "/path/to/file/with/gene/to/number/of/cytosines.txt"

genes2xyzWithTEG <- f.get.genes2xyz.via.files(gf2dmcFile, g2gfFile, g2numFile)

# remove transposable_element_gene from the data (each transposable_element_gene is
# within a transposable_element anyway)
TEGtoTEfile <- "/path/to/a/file/with/two/columns/TEG-AGI/and/TE-AGIs.txt" # tab sep
genes2xyzWithoutTEG <- f.remove.tegs.from.genes2xyz(genes2xyzWithTEG, TEGtoTEfile)

f.plot.number.DMCs.per.gene(genes2xyzWithTEG, rDir, "a prefix")
f.plot.number.DMCs.per.gene(genes2xyzWithoutTEG, rDir, "a prefix")

######################################################################################
## check mapping stats
######################################################################################

mappingStats <- read.table(file.path(dataDir, "significantPositions.stats"), sep='\t',
                           header=T, row.names=NULL, fill=T, stringsAsFactors = FALSE)
mappingStats <- mappingStats[,c("context", "feature", "methylatedCs", "unmethylatedCs",
                                "methylatedCoverage", "unmethylatedCoverage")]
colnames(mappingStats) <- c("context", "feature", "meCnum", "deCnum", "methCov", "demethCov")
mappingStats$totalC <- mappingStats$meCnum + mappingStats$deCnum
mappingStats$perc <- round(mappingStats$totalC/sum(mappingStats$totalC)*100,3)

svg(file.path(rDir, "dmcPieCharts.svg"), height = 10, width = 10)
layout(matrix(1:4, ncol = 2))
sry <- aggregate(cbind(totalC, meCnum, deCnum) ~ context, data = mappingStats, sum)
pie(c(13, 15, 72), labels = c("CG", "CHG", "CHH"), col = c("gray10", "gray40", "gray80"), main = "Cs")
pie(sry$totalC, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "all Cs")
pie(sry$meCnum, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "methylated Cs")
pie(sry$deCnum, labels = sry$context, col = c("gray10", "gray40", "gray80"), main = "unmethylated Cs")
dev.off()

contextColors <- data.frame(
  CG = brewer.pal(3,"YlOrBr"),
  CHG = brewer.pal(3,"PuBuGn"),
  CHH = brewer.pal(3,"BuGn"),
  row.names = c("totalC", "meCnum", "deCnum"), stringsAsFactors = FALSE
)
segPlots <- paste(rep(c("CG", "CHG", "CHH"), each = 2), rep(c("meCnum", "deCnum"), 3), sep = '_')
segCols <- unique(mappingStats$feature)
numPlots <- matrix(0, nrow = length(segPlots), ncol = length(segCols), dimnames = list(segPlots, segCols))
for (ctxt in c("CG", "CHG", "CHH")) {
  temp <- subset(mappingStats, context == ctxt)
  for (pc in c("meCnum", "deCnum")) {
    numPlots[paste(ctxt, pc, sep = '_'), temp$feature] <- temp[[pc]]
  }
}
numPlots <- as.data.frame(numPlots)
numPlots$transposon <- numPlots$transposable_element + numPlots$transposable_element_gene

# rearrange manually for nice order and plot spider graphs (see colnames(numPlots))
numPlots <- numPlots[,c("intergenic", "transposon", "five_flank", "five_prime_UTR", "exon",
                        "splice", "three_prime_UTR", "three_flank", "pseudogene",
                        "miRNA_primary_transcript")]
labs <- c("Intergenic", "Transposon", "5'Upstream", "5'UTR", "Exon", "Intron",
          "3'UTR", "3'Downstream", "Pseudogene", "miRNA")

percPlots <- numPlots/apply(numPlots, 1, sum)*100
acrossPercPlots <- numPlots/sum(numPlots)*100
svg(file.path(rDir, "spiderPlots.svg"), height = 15, width = 10)
layout(matrix(1:9, ncol = 3, byrow = TRUE))
for (ctxt in c("CG", "CHG", "CHH")) {
  cR <- paste(ctxt, c("meCnum", "deCnum"), sep = '_')
  meCol <- contextColors["meCnum", ctxt]
  deCol <- contextColors["deCnum", ctxt]
  radial.plot(percPlots[cR,], labels = labs, main = ctxt, radial.lim = seq(0,60,10), 
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
  radial.plot(acrossPercPlots[cR,], labels = labs, main = ctxt, radial.lim = seq(0,60,10),
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
  radial.plot(numPlots[cR,], labels = labs, main = ctxt, 
              start = pi/2, line.col = c(meCol,deCol), poly.col = c(meCol, NA),
              rp.type = "p", show.grid.labels = 3, lwd = 3, point.symbols = 16,
              point.col = c(meCol,deCol), show.centroid = FALSE)
}
dev.off()
```
