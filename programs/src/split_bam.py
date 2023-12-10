import sys
import argparse
from collections import defaultdict
from multiprocessing import Pool, Lock
import re
from block_utils import *
import pysam

lockForWritingOutputBam = Lock()

class IteratorForSplitting:
    def __init__(self, inputPysamAlignmentFile, outputPysamAlignmentFile, tagsToKeep, blocks):
        self.inputPysamIterator = inputPysamAlignmentFile.fetch()
        self.inputPysamHeader = inputPysamAlignmentFile.header
        self.outputPysamAlignmentFile = outputPysamAlignmentFile
        self.tagsToKeep = tagsToKeep
        self.blocks = blocks
    def __iter__(self):
        return self
    def __next__(self):
        return self.inputPysamIterator.__next__(), \
            self.outputPysamAlignmentFile, \
            self.tagsToKeep, \
            self.blocks

def getPysamTags(pysamRecord, tagNames):
    newTagList = []
    for tagName, tagValue in pysamRecord.tags:
        if tagName in tagNames:
            newTagList.append((tagName, tagValue))
    return  newTagList

def writeSplitAlignments(inputPysamRecord,
                         splitRefBlocks,
                         splitQueryBlocks,
                         splitCigarLists,
                         outputPysamAlignmentFile,
                         tagsToKeep):
    splitIndex = 0
    for refBlock, queryBlock, cigarList in zip(splitRefBlocks, splitQueryBlocks, splitCigarLists):
        if queryBlock[0] < queryBlock[1]:
            splitPySamRecord = pysam.AlignedSegment()
            splitPySamRecord.query_name = inputPysamRecord.query_name + "_split_" + str(splitIndex)
            splitPySamRecord.query_sequence = inputPysamRecord.query[queryBlock[0]: queryBlock[1]]
            splitPySamRecord.flag = inputPysamRecord.flag
            splitPySamRecord.reference_id = inputPysamRecord.reference_id
            splitPySamRecord.reference_start = refBlock[0]
            splitPySamRecord.mapping_quality = inputPysamRecord.mapping_quality
            splitPySamRecord.cigar = createCigarListCompatibleWithPysam(cigarList)
            splitPySamRecord.next_reference_id = inputPysamRecord.next_reference_id
            splitPySamRecord.next_reference_start = inputPysamRecord.next_reference_start
            splitPySamRecord.template_length = inputPysamRecord.template_length
            splitPySamRecord.query_qualities = inputPysamRecord.qual[queryBlock[0]: queryBlock[1]]
            splitPySamRecord.tags = getPysamTags(inputPysamRecord, tagsToKeep)
            with lockForWritingOutputBam:
                outputPysamAlignmentFile.write(splitPySamRecord)
            splitIndex += 1
    return splitIndex

def splitOneAlignment(inputPysamRecord, outputPysamAlignmentFile, tagsToKeep, blocks):
    alignment = Alignment.createFromPysamRecord(inputPysamRecord)
    # Extract the alignment attributes like the contig name, alignment boundaries, orientation and cigar
    chromName = alignment.chromName
    contigName = alignment.contigName
    orientation = alignment.orientation
    if alignment.isPrimary == False:
        return [chromName, contigName, orientation, [], [], []]
    # rBlocks contains the projections and
    # qBlocks contains the projectable blocks
    if len(blocks[chromName]) == 0: # Continue if there is no block in the chrom
        return [chromName, contigName, orientation, [], [], []]
    projectableBlocks, projectionBlocks, projectionCigarLists = findProjections(
        mode='ref2asm',
        cigarList=alignment.cigarList,
        forwardBlocks=blocks[chromName],
        chromLength=alignment.chromLength,
        chromStart=alignment.chromStart + 1, # make 1-based start
        chromEnd=alignment.chromEnd,
        contigLength=alignment.contigLength,
        contigStart=alignment.contigStart + 1,
        contigEnd=alignment.contigEnd, # make 1-based start
        orientation=alignment.orientation,
        includeEndingIndel=False,
        includePostIndel=False
    )

    splitQueryBlocks = convertContigToQueryCoordinatesList(projectionBlocks, alignment)
    numberOfAlignments = writeSplitAlignments(inputPysamRecord = inputPysamRecord,
                                              splitRefBlocks = projectableBlocks,
                                              splitQueryBlocks = splitQueryBlocks,
                                              splitCigarLists = projectionCigarLists,
                                              outputPysamAlignmentFile = outputPysamAlignmentFile,
                                              tagsToKeep = tagsToKeep)
    return numberOfAlignments

def splitAlignmentsMultiThreaded(inputBam, blocks, outputBam, tagsToKeep, threads):
    inputPysamAlignmentFile = pysam.AlignmentFile(inputBam)
    outputPysamAlignmentFile = pysam.AlignmentFile(outputBam, header=inputPysamAlignmentFile.header)
    iteratorForSplitting = IteratorForSplitting(inputPysamAlignmentFile = inputPysamAlignmentFile,
                                                outputPysamAlignmentFile = outputPysamAlignmentFile,
                                                tagsToKeep = tagsToKeep,
                                                blocks = blocks)
    pool = Pool(threads)
    results = pool.starmap(splitOneAlignment, iteratorForSplitting)
    pool.close()
    return results

def main():
    parser = argparse.ArgumentParser(description='Split the records in BAM file based on the given BED file and create a new BAM file with split records ')
    parser.add_argument('--inputBam', type=str,
                        help='(BAM format) The alignments to be split')
    parser.add_argument('--blocks', type=str,
                        help='(BED format) The desired blocks in the reference coordinates to be used for splitting alignment records. It should be sorted and with no overlaps.')
    parser.add_argument('--outputBam', type=str,
                        help='(BAM format) A path for saving the output BAM file with split records.')
    parser.add_argument('--threads', type=int,
                        help='Number of threads')


    # Fetch the arguments
    args = parser.parse_args()
    inputBam = args.inputBam
    blocksBed = args.blocks
    outputBam = args.outputBam
    threads = args.threads

    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    blocks = defaultdict(list)
    with open(blocksBed,"r") as f:
        for line in f:
            if line.startswith("track name"):
                continue
            attrbs = line.strip().split()
            chromName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            info = attrbs[3:] if len(attrbs) > 3 else [""]
            blocks[chromName].append((start, end, info))

    results = splitAlignmentsMultiThreaded(inputBam = inputBam,
                                           blocks = blocks,
                                           outputBam = outputBam,
                                           tagsToKeep = ['HP'],
                                           threads = threads)
    totalNumberOfAlignments = sum(results)
    print("total number of written alignments  = %d", totalNumberOfAlignments)

if __name__ == "__main__": main()

