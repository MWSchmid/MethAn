# MethAn
A collection of methylome analysis tools

Note:
* this tools are provided "as-they-are"
* I used them for some methylome analysis in Arabidopsis - thus, the Python and R scripts:
    * contain some hard-coded parameters (e.g. chromosome names are 1, 2, 3, 4, 5, Mt and Pt)
    * are not optimized for large genomes (plot_generic_gene.py creates for each chromosome a list as long as the chromosome)
* MethAnMap, the binary for mapping positions to the annotation is an exception in this respect. It deals with any genome (as long as the annotation and the genome sequence match).

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

Note that the mapping statistics will make use of the priorities set with Rcount-format.

```SH
### Required files ###

# Input
myNucleotides="/path/to/a/file/with/nucleotides.txt"
myAnnotation="/path/to/an/annotation.xml"

# Output (the myRegionFile is optional)
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.txt"
myRegions="/path/to/a/file/with/the/regions.txt"
myMapStats="/path/to/a/file/with/the/mapping/statistics.txt"

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
awk -v OFS='\t' '{print $1,$3,$7 > "forMetagene_ALL.txt"; print $1,$3,$7 > "forMetagene_"$6".txt"}' ${myMappedNucleotides}
cd /path/to/MethAn

for CONTEXT in ALL CG CHG CHH; do
inputFile="${workingDirectory}/forMetagene_${CONTEXT}.txt"
outputFile="${workingDirectory}/metagene_${CONTEXT}.svg"
python plot_generic_gene.py inputFile outputFile -fraction ${geneFraction} -window ${windowType} -window_len ${windowSize} -bordersize ${borderSize} -endcorrection ${endCorrection}
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
awk -v OFS='\t' '{print $1,$3 > "forMetagene_ALL.txt"; print $1,$3 > "forDistance_"$6".txt"}' ${myMappedNucleotides}
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

createRandomNucleotideSets.sh -d ${myNucleotides} ${mySigNucleotides} ${myAnnotation} ${numSets} ${randomDir} ${outSuffix}

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
```

## R functions

... add ...
