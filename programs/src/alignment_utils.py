import sys
import argparse
from collections import defaultdict
import re
from multiprocessing import Pool
from copy import deepcopy
import random
from Bio import SeqIO
import gzip
import math
from block_utils import *


def splitIntoAlignmentsWithShortGaps(alignment, maxGapLength):
    refStart = -1
    refEnd = -1
    contigStart = -1
    contigEnd = -1
    cigarList = []
    splitAlignments = []
    for op, opLen, refInterval, contigInterval in iterateCigar(alignment):
        if ((op == 'D' or op == 'I') and opLen > maxGapLength) or op == 'S' or op == 'H':
            if refStart < refEnd and contigStart < contigEnd:
                splitAlignment = Alignment(None)
                splitAlignment.contigName = alignment.contigName
                splitAlignment.contigLength = alignment.contigLength
                splitAlignment.contigStart = contigStart
                splitAlignment.contigEnd = contigEnd
                splitAlignment.orientation = alignment.orientation
                splitAlignment.chromName = alignment.chromName
                splitAlignment.chromStart = refStart
                splitAlignment.chromEnd = refEnd
                splitAlignment.chromLength = alignment.chromLength
                splitAlignment.isPrimary = alignment.isPrimary
                splitAlignment.cigarList = cigarList
                splitAlignment.editDistance = None
                splitAlignment.leftHardClip = None
                splitAlignment.rightHardClip = None
                # add the new split alignment
                splitAlignments.append(splitAlignment)
                # reset start and end coors
                refStart = -1
                refEnd = -1
                contigStart = -1
                contigEnd = -1
                cigarList = []
        else:
            if refStart == -1 or contigStart == -1:
                refStart = refInterval[0]
                contigStart = contigInterval[0]
            refEnd = refInterval[1]
            contigEnd = contigInterval[1]
            cigarList.append((op, opLen))
    # get the last piece of alignment if it exists
    if refStart < refEnd and contigStart < contigEnd:
        splitAlignment = Alignment(None)
        splitAlignment.contigName = alignment.contigName
        splitAlignment.contigLength = alignment.contigLength
        splitAlignment.contigStart = contigStart
        splitAlignment.contigEnd = contigEnd
        splitAlignment.orientation = alignment.orientation
        splitAlignment.chromName = alignment.chromName
        splitAlignment.chromStart = refStart
        splitAlignment.chromEnd = refEnd
        splitAlignment.chromLength = alignment.chromLength
        splitAlignment.isPrimary = alignment.isPrimary
        splitAlignment.cigarList = cigarList
        splitAlignment.editDistance = None
        splitAlignment.leftHardClip = None
        splitAlignment.rightHardClip = None
        # add the new split alignment
        splitAlignments.append(splitAlignment)

    return splitAlignments


def getRefBlockListPerChromFromAlignments(alignments, maxIndelSize=500, merge=True):
    blockListPerChrom = defaultdict(BlockList)
    for alignment in alignments:
        blockListPerChrom[alignment.chromName].extend(alignment.getRefCoveredBlockList(maxIndelSize))
    if merge:
        for chrom in blockListPerChrom:
            blockListPerChrom[chrom].sort()
            blockListPerChrom[chrom].mergeWithOverlapCount(inplace=True)
    return blockListPerChrom