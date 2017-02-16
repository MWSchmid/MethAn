#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mergeBedGraphs.py

This script can be used to merge BED graphs produced by bismark_SE/PE.sh and bamToBedGraph.sh and a sample annotation file (see below) into a single long-format table (for lm() in R).

NOTE: minCov >= 5, maxCov <= 100, and minGroupCount >= 2 are hard-coded.

sampleTable.csv has a header and contains the following columns:

<required>
\tsampleID
\tbedFile
\tgroup (should be the smallest grouping possible)

<optional>
\tother factors which will be used in the model later on

<example>
sampleID,bedFile,group,organ,sex
lm1,lm1.bed,liver_male,liver,male
lm2,lm2.bed,liver_male,liver,male
lm3,lm3.bed,liver_male,liver,male
lf1,lf1.bed,liver_female,liver,female
lf2,lf2.bed,liver_female,liver,female
lf3,lf3.bed,liver_female,liver,female
hm1,hm1.bed,heart_male,heart,male
hm2,hm2.bed,heart_male,heart,male
hm3,hm3.bed,heart_male,heart,male
hf1,hf1.bed,heart_female,heart,female
hf2,hf2.bed,heart_female,heart,female
hf3,hf3.bed,heart_female,heart,female

Notes:
\tgroup = <organ>_<sex>
\tthe bed files don't have to contain the sampleID, but you should provide the full path

Acknowledgements:
I would like to thank Önder Kartal. Most of this script is taken from his shannon script. See: gitlab.com/okartal/shannon

ToDo:
add requirements and a readme for the preprocessing.
"""

__author__ = "Marc W Schmid"
__version__ = "0"

import argparse
import collections
import os
import subprocess
import sys
import textwrap
import math

import numpy as np
import pandas as pd

import logging
logging.basicConfig(format="=== %(levelname)s === %(asctime)s === %(message)s", level=logging.DEBUG, datefmt='%Y-%m-%d %H:%M:%S')

def mergeBedGraphs(sampleTableFile, chrom, ctxt, outfileName, query=None, bedFileType="MethAn_bismark"):
    """The main driver function.
    """
    
    # read the sample table
    logging.info("Reading sample table.")
    try:
        sampleTab = pd.read_csv(sampleTableFile, comment='#', sep=',', header=0)
    except FileNotFoundError as e:
        logging.critical("Could not find the sampleTable.csv file.")
        sys.exit(66)

    # filter it if requested
    if query:
        logging.info("Selecting samples according to query.")
        sampleTab.query(query, inplace=True)
    else:
        logging.info("Selecting all samples.")
    
    # merge the bedGraphs
    logging.info("Merging BED files.")
    okCounter = 0
    chromCounter = 1
    if not os.path.isfile(chrom):
        worked = mergeFiles(sampleTab, chrom, ctxt, outfileName, bedFileType)
        if worked:
            okCounter += 1
    else:
        with open(chrom, "rb") as infile:
            chromList = [line[:-1].decode("ascii") for line in infile]
        chromCounter = len(chromList)
        for chrom in chromList:
            worked = mergeFiles(sampleTab, chrom, ctxt, outfileName, bedFileType)
            if worked:
                okCounter += 1
    
    # no return, just print a time information
    logging.info("Finished ("+str(okCounter)+"/"+str(chromCounter)+" with data/total)")
    pass


def mergeFiles(sampleTab, chrom, ctxt, outfileName, fileType="MethAn_bismark", imputation=False):
    """Merge individual BED files and the sample annotation into a long-format table.
    """

    inputFiles = sampleTab["bedFile"]
    labels = sampleTab["sampleID"]
    groups = sampleTab["group"]
    factorNames = list(sampleTab.columns.values)[3:]
    
    # initialize the output table
    logging.info("Initializing table.")
    
    column = collections.OrderedDict()
    if fileType == "MethAn_bismark":
        column[0] = "chrom"
        column[1] = "pos"
        #column[2] = "end"
        #column[3] = "pcentmeth"
        column[4] = 'M' # methylated coverage
        column[5] = 'U' # unmethylated coverage
        types = [str, np.int64, np.int32, np.int32]
        dtype = dict(zip(column.values(), types))
    else:
        logging.critical("User-provided file type is not implemented (see --bedFileType).")
        sys.exit(64)
        
    # Get the query region:
    # load is average number of expected sites in the region, this determines
    # the memory requirements
    logging.info("Preparing for sequential reading.")
    LOAD = 1e6
    
    supChrom = chromSupremum(inputFiles, chrom)
    if supChrom == None:
        logging.info("Skipping because chromosome is missing.")
        return False
    supNsites = nsitesSupremum(inputFiles, chrom)
    if supNsites == None or supNsites == 0:
        logging.info("Skipping because there are no entries.")
        return False
    stepSize = math.ceil(supChrom/supNsites * LOAD)

    if stepSize < supChrom:
        step = stepSize
        logging.info("step size: "+format(step))
    else:
        step = supChrom
        logging.info("step size: "+format(step)+" (max. for contig {0}).".format(chrom))

    posStart = list(range(0, supChrom, step + 1))
    posEnd = list(range(step, supChrom, step + 1)) + [supChrom]
    regions = zip(posStart, posEnd)

    logging.info("Merging data.")
    for interval in regions:
        region = chrom + ':{0}-{1}'.format(*interval)
        
        # load into data frame
        logging.info("Loading "+region)
        tabixQuery = (subprocess.Popen(['tabix', f, region],
                                        stdout=subprocess.PIPE,
                                        universal_newlines=True)
                       for f in inputFiles)
        dataframes = (pd.read_table(query.stdout, comment='#', header=None,
                                    usecols=list(column.keys()),
                                    names=list(column.values()), dtype=dtype)
                      for query in tabixQuery)
        
        # add annotation
        logging.info("Adding annotation to "+region)
        rowNums = range(len(labels))
        reformDFs = []
        for rowNum, sample, group, df in zip(rowNums, labels, groups, dataframes):
            if ctxt != "FILE":
                df["ctxt"] = ctxt
            df["coverage"] = df['M']+df['U']
            df["pcentmeth"] = df['M']/df["coverage"]
            #df["pcentmeth"] = 0.0
            df["sample"] = sample
            df["group"] = group
            del df['M']
            del df['U']
            curSample = sampleTab["sampleID"].get_value(rowNum)
            for fn in factorNames:
                df[fn] = sampleTab[fn].get_value(rowNum)
            reformDFs.append(df)
            
        # merge
        logging.info("Merging "+region)
        mergedTab = pd.concat(reformDFs, ignore_index = True)

        # sort
        logging.info("Sorting "+region)
        mergedTab.sort_values(["pos"], ascending=True, inplace=True)
        
        # filter
        logging.info("Filtering "+region)
        #mergedTab = mergedTab[mergedTab["cov"]>=5]
        #mergedTabGrouped = mergedTab.groupby(by=["pos", "group"])
        #groupCounts = mergedTabGrouped.apply(lambda g: len(g))
        minGroupCount = 2
        groupCounts = mergedTab.query("coverage>=5 & coverage<=100").groupby(by=["pos", "group"])["coverage"].size()
        posCounts = groupCounts[groupCounts>=minGroupCount].count(level="pos")
        okPositions = list(posCounts[posCounts == len(set(groups))].index)
        selectedPositions = mergedTab[mergedTab["pos"].isin(okPositions)]
        
        # add the output - write the header in in case the file does not exist
        logging.info("Writing "+region)
        writeHeader = not os.path.isfile(outfileName)
        selectedPositions.query("coverage>=5 & coverage<=100").to_csv(outfileName, header=writeHeader, sep='\t', index=False, mode='a')
        
        # write also all positions that have a coverage of at least 5 and max 100
        outfileNameNoGroupFilter = outfileName+".noGroupFilter"
        writeHeader = not os.path.isfile(outfileNameNoGroupFilter)
        mergedTab.query("coverage>=5 & coverage<=100").to_csv(outfileNameNoGroupFilter, header=writeHeader, sep='\t', index=False, mode='a')
        
    return True


def chromSupremum(tabixfiles, chrom):
    """Return the least upper bound for the chrom end coordinate.
    """

    end_coordinate = list()

    for f in tabixfiles:
        tabix = subprocess.Popen(["tabix", f, chrom], stdout=subprocess.PIPE)
        tail = subprocess.Popen(["tail", "-1"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        cut = subprocess.Popen(["cut", "-f3"], stdin=tail.stdout, stdout=subprocess.PIPE)
        tabix.stdout.close()  # Allow first process to receive a SIGPIPE if process 2 exits.
        try:
            base_position = int(cut.communicate()[0])
            end_coordinate.append(base_position)
        except ValueError:
            continue
    
    try:
        out = np.max(end_coordinate)
    except ValueError:
        out = None
    return out


def nsitesSupremum(tabixfiles, chrom):
    """Return the least upper bound for the number of covered sites.
    """
    
    sites = list()

    for f in tabixfiles:
        tabix = subprocess.Popen(["tabix", f, chrom], stdout=subprocess.PIPE)
        wcl = subprocess.Popen(["wc", "-l"], stdin=tabix.stdout, stdout=subprocess.PIPE)
        tabix.stdout.close()  # Allow tabix to receive a SIGPIPE if wcl exits.
        try:
            site_count = int(wcl.communicate()[0])
            sites.append(site_count)
        except ValueError:
            continue

    try:
        out = np.max(sites)
    except ValueError:
        out = None
    return out


def impute(data, method='pseudocount'):
    """This should not be used - it's just a placeholder.
    """
    logging.critical("DO NOT USE IMPUTATION.")
    sys.exit(64)
    if method == 'pseudocount':
        value = 1
    return value


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog="mergeBedGraphs.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            Merge BED graphs produced by bismark_SE/PE.sh and
            bamToBedGraph.sh and a sample annotation file into
            a single long-format table (for lm() in R).
            
            NOTE: minCov >= 5, maxCov <= 100, and minGroupCount >= 2 are hard-coded.
            ======================================================
            """),
        epilog=textwrap.dedent("""\
            ======================================================
            Format of the sampleTable:
            
            <required columns>
            * sampleID
            * bedFile
            * group (should be the smallest grouping possible)

            <optional columns>
            * other factors which will be used in the model later on

            <example>
            sampleID,bedFile,group,organ,sex
            lm1,lm1.bed,liver_male,liver,male
            lm2,lm2.bed,liver_male,liver,male
            lm3,lm3.bed,liver_male,liver,male
            lf1,lf1.bed,liver_female,liver,female
            lf2,lf2.bed,liver_female,liver,female
            lf3,lf3.bed,liver_female,liver,female
            hm1,hm1.bed,heart_male,heart,male
            hm2,hm2.bed,heart_male,heart,male
            hm3,hm3.bed,heart_male,heart,male
            hf1,hf1.bed,heart_female,heart,female
            hf2,hf2.bed,heart_female,heart,female
            hf3,hf3.bed,heart_female,heart,female

            Notes:
            * group = <organ>_<sex>
            * the bed files do not have to contain the sampleID,
              but you should provide the full path
              
            ======================================================
            Acknowledgements:
            I would like to thank Önder Kartal. Substantial parts
            of this script are taken from him.
            See: gitlab.com/okartal/shannon
            """))

    parser.add_argument("-v", "--version", action="version",
                        version='%(prog)s {0}'.format(__version__))
    
    parser.add_argument("chromosome", metavar="chromOrFile", type=str,
                        help="""
                        Chromosome name or path to a file with several
                        chromosome names (one per line).
                        """)

    parser.add_argument("nucleotideContext", metavar="ctxt", type=str, 
                        help="Nucleotide context (CG, CHG, CHH, FILE).")

    parser.add_argument("sampleTable", type=str,
                        help="""
                        A csv table with metadata for each sample/record.
                        Lines starting with # are ignored. The first line
                        is interpreted as the header. Details and example below.
                        """)

    parser.add_argument("outputFile", type=str,
                        help="Name of the output file.")

    parser.add_argument("-q", "--query", metavar='"STR"', type=str, required=False,
                        help="""
                        Query expression to select a subset of samples. The
                        expression has to be in double quotes. Examples: "organ ==
                        'heart'", "age >= 32".
                        """)

    parser.add_argument("-b", "--bedFileType", metavar='[MethAn_bismark]', type=str, required=False,
                        default = "MethAn_bismark",
                        help="""
                        Describe the type of the bed file. Only MethAn_bismark
                        is currently implemented (default).
                        """)
    
    args = parser.parse_args()

    mergeBedGraphs(sampleTableFile=args.sampleTable,
                   chrom=args.chromosome,
                   ctxt=args.nucleotideContext,
                   outfileName=args.outputFile,
                   query=args.query,
                   bedFileType=args.bedFileType)



