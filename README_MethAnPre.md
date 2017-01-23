The scripts in [MethAnPre](MethAnPre) can be used for trimming, alignment, duplicate removal and methylation extraction. The scripts require:

* [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [picard](https://broadinstitute.github.io/picard)
* [Pile-O-Meth](https://bioconda.github.io/recipes/pileometh/README.html)
* [samtools](https://github.com/samtools/samtools)
* [tabix and bgzip](http://www.htslib.org/doc/tabix.html)


# Install requirements for MethAnPre

```SH
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install unzip build-essential zlibc zlib1g zlib1g-dev tabix git

# install miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# add conda bin to path
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

**TODO** PileOMeth is now MethylDackel - check and update this README


# Download and index a genome, trim, align, de-duplicate, extract methylation

```SH
# some paths and variables
GENOME_FOLDER="/path/to/the/folder/with/the/genome"
GENOME_FASTA="At.fasta"
INDIR="/path/to/the/folder/containing/the/reads"
OUTDIR="/path/to/the/folder/where/alignments/and/counts/will/be/stored"

# add conda and MethAnPre to the path env-var
export PATH="$HOME/miniconda3/bin:$HOME/MethAn/MethAnPre:$PATH"

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


# Merge BedGraphs and annotation into a long-format table

Merge BED graphs produced by bismark_SE/PE.sh and bamToBedGraph.sh and a sample annotation file into a single long-format table (for lm() in R). The format of the sample annotation table (taken from the mergeBedGraphs.py help):

* <required columns>
    * sampleID
    * bedFile
    * group (should be the smallest grouping possible)

* <optional columns>
    * other factors which will be used in the model later on

* <example>
    * sampleID,bedFile,group,organ,sex
    * lm1,lm1.bed,liver_male,liver,male
    * lm2,lm2.bed,liver_male,liver,male
    * lm3,lm3.bed,liver_male,liver,male
    * lf1,lf1.bed,liver_female,liver,female
    * lf2,lf2.bed,liver_female,liver,female
    * lf3,lf3.bed,liver_female,liver,female
    * hm1,hm1.bed,heart_male,heart,male
    * hm2,hm2.bed,heart_male,heart,male
    * hm3,hm3.bed,heart_male,heart,male
    * hf1,hf1.bed,heart_female,heart,female
    * hf2,hf2.bed,heart_female,heart,female
    * hf3,hf3.bed,heart_female,heart,female

* Notes:
    * group = <organ>_<sex> (generally the most detailed grouping possible)
    * the bed files do not have to contain the sampleID in their names, but you should provide the full path

The example below assumes that you have three tables which are identical except for the file paths (different because of the three separate and different contexts):

* sampleTableForMerge_CpG.csv
* sampleTableForMerge_CHG.csv
* sampleTableForMerge_CHH.csv

For each context (or chromosome and context if kept separate), the script will produce two files:
* one only with positions having a coverage of at least 5 within at least one replicate per group (name as specified on the command line - e.g. merged_CpG.txt).
* one only with positions having a coverage of at least 5 in any sample (name as specified on the command line plus ending ".noGroupFilter" - e.g. merged_CpG.txt.noGroupFilter)


```SH
# extract all chromosome identifiers
GENOME_FOLDER="/path/to/the/folder/with/the/genome"
GENOME_FASTA="At.fasta"
grep ">" $GENOME_FOLDER/$GENOME_FASTA | awk '{print $1}' | sed 's/>//g' > allChroms.txt

# add MethAnPre to the path env-var
export PATH="$HOME/MethAn/MethAnPre:$PATH"

# either: keep chromosomes separate
ALLCHROMS=($(<allChroms.txt))
for CHROM in "${ALLCHROMS[@]}"; do
for CTXT in CpG CHG CHH; do
mergeBedGraphs.py $CHROM $CTXT sampleTableForMerge_${CTXT}.csv merged_${CTXT}_{CHROM}.txt
done
done

# or: merge all chromosomes into one file
for CTXT in CpG CHG CHH; do
mergeBedGraphs.py allChroms.txt $CTXT sampleTableForMerge_${CTXT}.csv merged_${CTXT}.txt
done
```


# Run a linear model for each cytosine

**TODO**
