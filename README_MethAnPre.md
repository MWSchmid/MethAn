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

**TODO**


# Run a linear model for each cytosine

**TODO**
