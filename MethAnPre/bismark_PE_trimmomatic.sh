#!/bin/bash
#
# bismark_PE_trimmomatic.sh -- Align PE reads with bismark.
#
# Authors: Marc W. Schmid <marcschmid@gmx.ch>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

me=$(basename "$0")

## defaults
inputDir=""
outputDir=""
prefix=""
inputFile=""
inputFileReverse=""
genomeFolder=""
directional=""
pbat=""
threads=8
memory=3 # samtools frequently underestimates RAM usage
maxMem=4
verbose=""
trimmomatic="$HOME/Trimmomatic-0.36/trimmomatic-0.36.jar"
adapterSeqs="$HOME/Trimmomatic-0.36/adapters/TruSeq3-PE.fa"

## Exit status codes (following <sysexits.h>)
EX_OK=0			# successful termination
EX_USAGE=64		# command line usage error
EX_DATAERR=65		# data format error
EX_NOINPUT=66		# cannot open input
EX_NOUSER=67		# addressee unknown
EX_NOHOST=68		# host name unknown
EX_UNAVAILABLE=69	# service unavailable
EX_SOFTWARE=70		# internal software error
EX_OSERR=71		# system error (e.g., can't fork)
EX_OSFILE=72		# critical OS file missing
EX_CANTCREAT=73		# can't create (user) output file
EX_IOERR=74		# input/output error
EX_TEMPFAIL=75		# temp failure; user is invited to retry
EX_PROTOCOL=76		# remote error in protocol
EX_NOPERM=77		# permission denied
EX_CONFIG=78		# configuration error

## helper functions

function die () {
    rc="$1"
    shift
    (echo -n "$me: ERROR: ";
        if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
    exit $rc
}

## usage info

usage () {
    cat <<__EOF__
Usage:
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [FASTQ-forward] [FASTQ-reverse] [GENOMEFOLDER]
This script aligns Illumina reads (PE) with bismark (bowtie2). You should use a machine with at least 8 cores.
Arguments:
INDIR: Directory with the input file (<extension/type>).
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. Files will be:
<OUTPREFIX>.dupMark.sorted.bam
<OUTPREFIX>.dupMark.sorted.bam.bai
<OUTPREFIX>_markdup_metrics.txt
FASTQ-forward: Name of the .fq(.gz) file with the forward reads.
FASTQ-reverse: Name of the .fq(.gz) file with the reverse reads.
GENOMEFOLDER: Folder with the genome (see bismark_genome_preparation).
Options:
  -v        enable verbose logging (no effect)
  -h        print this help text
  -t        number of available threads (no effect)
  -m        amount of memory to be allocated (no effect)
  -n        library is notDirectional (off, see bismark user guide)
  -p        it's a PBAT library (off, see bismark user guide)
  -j		Path to trimmomatic-xx.jar
  -a		Path to adapter sequences
__EOF__
}

warn () {
  (echo -n "$me: WARNING: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die $EX_UNAVAILABLE "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

input_exists () {
  echo -n "Checking input file ${1}..."
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_NOINPUT "Could not find input file: ${1}."
  fi
}

output_exists () {
  echo -n "Checking output file ${1}... "
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_OSFILE "Could not find output file: ${1}."
  fi
}

remove_if_present () {
  if [ -e $1 ]; then
    rm $1
  fi
}

## parse command-line

short_opts='hvt:m:npj:a:'
long_opts='help,verbose,threads,memory,notDirectional,pbat,jarTrimmomatic,adapterSeqs'

getopt -T > /dev/null
rc=$?
if [ "$rc" -eq 4 ]; then
    # GNU getopt
    args=$(getopt --name "$me" --shell sh -l "$long_opts" -o "$short_opts" -- "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    # use 'eval' to remove getopt quoting
    eval set -- $args
else
    # old-style getopt, use compatibility syntax
    args=$(getopt "$short_opts" "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    set -- $args
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --threads|-t)  shift; threads=$1 ;;
        --memory|-m)   shift; memory=$1 ;;	
        --verbose|-v)  verbose=" --verbose" ;;
        --help|-h)     usage; exit 0 ;;
        --notDirectional|-d) directional=" --non_directional" ;;
        --pbat|-p)     pbat=" --pbat";;
    	--adapterSeqs|-a)	shift; adapterSeqs=$1;;
    	--jarTrimmomatic|-j)	shift; trimmomatic=$1 ;;
        --)            shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 6 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

inputDir=$1
shift
outputDir=$1
shift
prefix=$1
shift
inputFile=$1
shift
inputFileReverse=$1
shift
genomeFolder=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

export PATH="/home/ubuntu/miniconda3/bin:$PATH"

require_command bowtie2
require_command bismark
require_command samtools
require_command picard

# checking input and trimmomatic
input_exists ${inputDir}/${inputFile}
input_exists ${inputDir}/${inputFileReverse}
input_exists ${trimmomatic}

### run script

command="java -Xmx${maxMem}g -jar ${trimmomatic} PE -threads ${threads} -phred33\
	${inputDir}/${inputFile} ${inputDir}/${inputFileReverse}\
	${outputDir}/${inputFile%%.*}_paired.tr.fq.gz ${outputDir}/${inputFile%%.*}_unpaired.tr.fq.gz\
	${outputDir}/${inputFileReverse%%.*}_paired.tr.fq.gz ${outputDir}/${inputFileReverse%%.*}_unpaired.tr.fq.gz\
	ILLUMINACLIP:${adapterSeqs}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
remove_if_present ${outputDir}/"${inputFile%%.*}_unpaired.tr.fq.gz"
remove_if_present ${outputDir}/"${inputFileReverse%%.*}_unpaired.tr.fq.gz"

command="gunzip ${outputDir}/${inputFile%%.*}_paired.tr.fq.gz"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

command="gunzip ${outputDir}/${inputFileReverse%%.*}_paired.tr.fq.gz"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

inputFile="${inputFile%%.*}_paired.tr.fq"
inputFileReverse="${inputFileReverse%%.*}_paired.tr.fq"
input_exists ${outputDir}/${inputFile}
input_exists ${outputDir}/${inputFileReverse}

command="bismark --bowtie2${directional}${pbat} --seedmms 1 -o ${outputDir} --basename ${prefix} ${genomeFolder} -1 ${outputDir}/${inputFile} -2 ${outputDir}/${inputFileReverse}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
remove_if_present ${outputDir}/${inputFile}
remove_if_present ${outputDir}/${inputFileReverse}
output_exists "${outputDir}/${prefix}_pe.bam"

command="samtools sort -O bam -@ ${threads} -T tempoSort -m ${memory}G -o ${outputDir}/${prefix}.sorted.bam ${outputDir}/${prefix}_pe.bam"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
remove_if_present ${outputDir}/${prefix}_pe.bam
output_exists "${outputDir}/${prefix}.sorted.bam"

command="samtools index ${outputDir}/${prefix}.sorted.bam"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
output_exists "${outputDir}/${prefix}.sorted.bam"
output_exists "${outputDir}/${prefix}.sorted.bam.bai"

command="picard MarkDuplicates -Xmx${maxMem}g I=${outputDir}/${prefix}.sorted.bam O=${outputDir}/${prefix}.dupMarked.sorted.bam M=${outputDir}/${prefix}_markdup_metrics.txt ASSUME_SORTED=true"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
output_exists "${outputDir}/${prefix}.dupMarked.sorted.bam"
output_exists "${outputDir}/${prefix}_markdup_metrics.txt"

#command="samtools sort -O bam -@ ${threads} -T tempoSort -m ${memory}G -o ${outputDir}/${prefix}.dupMarked.sorted.bam ${outputDir}/${prefix}.dupMarked.bam"
#echo "=== ${me}: Running: ${command}"
#eval $command
#rc=$?
#echo "=== ${me}: Command ended with exit code $rc"
#output_exists "${outputDir}/${prefix}.dupMarked.sorted.bam"

command="samtools index ${outputDir}/${prefix}.dupMarked.sorted.bam"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

## Checking final output
output_exists "${outputDir}/${prefix}.dupMarked.sorted.bam"
output_exists "${outputDir}/${prefix}.dupMarked.sorted.bam.bai"
output_exists "${outputDir}/${prefix}_markdup_metrics.txt"

## Remove all remaining intermediary files
remove_if_present ${outputDir}/${prefix}.sorted.bam
remove_if_present ${outputDir}/${prefix}.sorted.bam.bai
#remove_if_present ${outputDir}/${prefix}.dupMarked.bam

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
