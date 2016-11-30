#!/bin/bash
#
# createRandomNucleotideSets.sh -- create sets of random nucleotides
# sampling is done on a per-context probability
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
myFolder=$(dirname "$0")

## defaults
allNucleotides=""
sigNucleotides=""
annotationFile=""
numberOfSets=10
outputDir=""
outSuffix=""
doDistance="no"

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
  $me [options] [ALLNUCLEOTIDES] [SIGNUCLEOTIDES] [ANNOTATIONFILE] [NUMBEROFSETS] [OUTPUTDIR] [OUTSUFFIX]

This script generates random sets of nucleotides. Sampling probabilities are context specific.
Use -d to create the files for the get_distane.py script.

Arguments:
ALLNUCLEOTIDES: file from which the nucleotides will be sampled (normally all possible nucleotides)
SIGNUCLEOTIDES: file from which the probabilities will be extracted (e.g. significant nucleotides)
ANNOTATIONFILE: an annotation file (Rcount XML)
NUMBEROFSETS: number of random sets
OUTPUTDIR: directory where the random sets will be stored
OUTSUFFIX: a suffix for the output file names (without extension)

see MethAnMap for details on the nucleotide files.

It's assumed that this script is located in the same folder as the folder "MethAnMap".

Options:
  -h            Print this help text
  -d            Create the files for the distance-to-feature script.
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

short_opts='hd'
long_opts='help,dist'

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
        --dist|-d)     doDistance="yes" ;;
        --help|-h)     usage; exit 0 ;;
        --)            shift; break ;;
    esac
    shift
done

## sanity checks

# the required arguments must be present 
if [ $# -lt 6 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

allNucleotides=$1
shift
sigNucleotides=$1
shift
annotationFile=$1
shift
numberOfSets=$1
shift
outputDir=$1
shift
outSuffix=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command $myFolder/MethAnMap/MethAnMap

# checking input
input_exists $allNucleotides
input_exists $sigNucleotides
input_exists $annotationFile

# get the probabilities for the sampling
TOTAL_CG=$(awk '{if ($3 == "CG") {SIGSUM+=1}} END {print SIGSUM}' "$allNucleotides")
TOTAL_CHG=$(awk '{if ($3 == "CHG") {SIGSUM+=1}} END {print SIGSUM}' "$allNucleotides")
TOTAL_CHH=$(awk '{if ($3 == "CHH") {SIGSUM+=1}} END {print SIGSUM}' "$allNucleotides")
SIGSUM_CG=$(awk '{if ($3 == "CG") {SIGSUM+=1}} END {print SIGSUM}' "$sigNucleotides")
FRAC_CG=$(bc<<<"scale=10;$SIGSUM_CG/$TOTAL_CG")
SIGSUM_CHG=$(awk '{if ($3 == "CHG") {SIGSUM+=1}} END {print SIGSUM}' "$sigNucleotides")
FRAC_CHG=$(bc<<<"scale=10;$SIGSUM_CHG/$TOTAL_CHG")
SIGSUM_CHH=$(awk '{if ($3 == "CHH") {SIGSUM+=1}} END {print SIGSUM}' "$sigNucleotides")
FRAC_CHH=$(bc<<<"scale=10;$SIGSUM_CHH/$TOTAL_CHH")
echo -e "context\tnumber significant\tsampling probability"
echo -e "CG:\t${SIGSUM_CG}\t${FRAC_CG}"
echo -e "CHG:\t${SIGSUM_CHG}\t${FRAC_CHG}"
echo -e "CHH:\t${SIGSUM_CHH}\t${FRAC_CHH}"

# generate random sets - check for the output at each step
for RSN in $(seq 1 "${numberOfSets}"); do
RS="random_${RSN}"
echo "processing ${RS}..."
myMappedNucleotides="$outputDir/${RS}${outSuffix}.mappedNucleotides"
myMapStats="$outputDir/${RS}${outSuffix}.stats"
#$myFolder/MethAnMap/MethAnMap -Z "${FRAC_CG},${FRAC_CHG},${FRAC_CHH}" -S ${allNucleotides} -A ${annotationFile} -O ${myMappedNucleotides} > ${myMapStats}
output_exists $myMappedNucleotides
output_exists $myMapStats
if [[ $doDistance == yes ]]; then
awk -v r="${outputDir}/${RS}${outSuffix}" '{print $1"\t"$3 > r"_forDistance_ALL.txt"; print $1"\t"$3 > r"_forDistance_"$6".txt"}' ${myMappedNucleotides}
fi
done

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."
exit $rc
