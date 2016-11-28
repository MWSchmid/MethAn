# MethAn
A collection of methylome analysis tools

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

## Extract some useful files for the analysis in R
```SH
# Input
myMappedNucleotides="/path/to/file/with/mapped/nucleotides.txt"

# Output
gf2dmc="/path/to/file/with/gene-feature/to/number/of/cytosines.txt"
g2gf="/path/to/file/with/gene/to/gene-feature.txt"
g2num="/path/to/file/with/gene/to/number/of/cytosines.txt"

# Run makeMappingFiles.py
python makeMappingFiles.py ${myMappedNucleotides} ${gf2dmc} ${g2gf} ${g2num}
```

## R functions

... add ...
