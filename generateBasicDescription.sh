#!/bin/bash
#
# generateBasicDescription
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
outputDir=""
annotation=""
allNucleotides=""
selectedNucleotides=""
chromFile=""
gffFile=""
TEfile=""
species=""
verbose=""
threads=1
memory=4
maxMem=4

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
  $me [options] [OUTPUTDIR] [ANNOTATION] [ALL NUCLEOTIDES] [SELECTED NUCLEOTIDES] [CHROMFILE] [GFFFILE] [TEFILE] [SPECIES]
Generate some descriptive tables and figures for a set of genomic positions:
Arguments:
OUTPUTDIR: directory where all the output will be stored
ANNOTATION: an Rcount XML annotation file
ALL NUCLEOTIDES: a table containing a reference set of nucleotides (e.g. all which passed the first filter)
SELECTED NUCLEOTIDES: a table with selected nucleotides - must be a subset of all nucleotides
CHROMFILE (for generic genes): a file with chromName<tab>chromSize
GFFFILE (for generic genes): a gff file (e.g. the one used for the Rcount XML annotation)
TEFILE (for overview in R): a file with two columns, tab-sep. First column should contain transposon/repeat identifiers
SPECIES: either At for Arabidopsis thaliana or Mp for Marchantia polymorpha
Options:
  -v        enable verbose logging (no effect)
  -h        print this help text
  -t        number of available threads (no effect)
  -m        amount of memory to be allocated (no effect)
DESCRIBE
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

short_opts='hvt:m:'
long_opts='help,verbose,threads,memory'

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
        --)            shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 8 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

outputDir=$1
shift
annotation=$1
shift
allNucleotides=$1
shift
selectedNucleotides=$1
shift
chromFile=$1
shift
gffFile=$1
shift
TEfile=$1
shift
species=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

# check if the required programs are around
require_command MethAnDirectMap
require_command makeMappingFiles.py
require_command plot_generic_gene_anyOrganism.py
require_command generate_overview.R
require_command generate_mapping_stats_figures.R

# checking input
input_exists ${annotation}
input_exists ${allNucleotides} 
input_exists ${selectedNucleotides} 
input_exists ${chromFile} 
input_exists ${gffFile} 
input_exists ${TEfile} 
input_exists "MethAnR.R"

# create output directory
mkdir -p ${outputDir}

# define output files
allMappedNucleotides="${outputDir}/all.mappedNucleotides"
allMapStats="${outputDir}/all.mapStats"
selectedMappedNucleotides="${outputDir}/selected.mappedNucleotides"
selectedMapStats="${outputDir}/selected.mapStats"
gf2dmc="${outputDir}/gf2dmc.txt"
g2gf="${outputDir}/g2gf.txt"
g2num="${outputDir}/g2num.txt"

# Map all cytosines in the files and check output
command="MethAnDirectMap -S ${allNucleotides} -A ${annotation} -O ${allMappedNucleotides} > ${allMapStats}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
command="MethAnDirectMap -S ${selectedNucleotides} -A ${annotation} -O ${selectedMappedNucleotides} > ${selectedMapStats}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
output_exists ${allMappedNucleotides}
output_exists ${allMapStats}
output_exists ${selectedMappedNucleotides}
output_exists ${selectedMapStats}

# Generate mapping files
command="makeMappingFiles.py ${selectedMappedNucleotides} ${gf2dmc} ${g2gf} ${g2num}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
output_exists ${gf2dmc}
output_exists ${g2gf}
output_exists ${g2num}

# Plot metagenes
borderSize=1000
geneFraction=100
windowType=blackman
windowSize=100
endCorrection=100
for CONTEXT in ALL CpG CHG CHH; do
inputFile="${outputDir}/forMetagene_${CONTEXT}.txt"
awk -v OFS='\t' -v CON=$CONTEXT '{if ($6==CON) {print $1,$3,$7} else {if(CON=="ALL"){print $1,$3,$7}}}' ${selectedMappedNucleotides} > ${inputFile}
outputFile="${outputDir}/metagene_${CONTEXT}.svg"
input_exists ${inputFile}
command="plot_generic_gene_anyOrganism.py ${chromFile} ${gffFile} ${inputFile} ${outputFile} -fraction ${geneFraction} -window ${windowType} -window_len ${windowSize} -bordersize ${borderSize} -endcorrection ${endCorrection}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
output_exists ${outputFile}
remove_if_present ${inputFile}
done

# Mapping stats in R
command="generate_mapping_stats_figures.R ${outputDir} all ${allMapStats} ${species}" 
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

# Mapping stats in R
command="generate_mapping_stats_figures.R ${outputDir} selected ${selectedMapStats} ${species}" 
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

# Overview in R
allPosFile="${outputDir}/tempAllPos.txt"
awk -v OFS='\t' '{print $1,$2}' ${allNucleotides} > ${allPosFile}
command="generate_overview.R ${outputDir} default ${allPosFile} ${selectedMappedNucleotides} ${gf2dmc} ${g2gf} ${g2num} ${TEfile}" 
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
# add file checks

# remove the temporary positions
remove_if_present ${allPosFile}

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
