import sys
import argparse
from collections import defaultdict
from multiprocessing import Pool, Lock, Manager
import re
from block_utils import *
import pysam
import copy
import os
import datetime


class PysamIteratorForBlockList:
    """ It takes a pysam.AlignmentFile object and a list of regions saved 
        in a dictionary of BlockLists. As an iterator it will iterate over the overlapping reads 
    """
    def __init__(self, inputPysamAlignmentFile, blockListPerContig):
        self.inputPysamAlignmentFile = inputPysamAlignmentFile
        self.blockListPerContig = blockListPerContig
        self.contigList = list(blockListPerContig.keys())
        self.contigIndex = -1
        self.blockIndex = -1
        contigName, start, end = self._get_next_region()
        if contigName is not None:
            self.inputPysamIterator = self.inputPysamAlignmentFile.fetch(contigName, start, end)
        else:
            self.inputPysamIterator = None
    
    def _get_next_region(self):
        self.contigIndex, self.blockIndex = BlockList.getNextIndicesBlockListPerContig(self.blockListPerContig,
                                                                                       self.contigList,
                                                                                       self.contigIndex, 
                                                                                       self.blockIndex)
        if self.contigIndex == None:
            return None, None, None
        contigName = self.contigList[self.contigIndex]
        block = self.blockListPerContig[contigName].blocks[self.blockIndex]
        return contigName, block[0] - 1, block[1]

    def __iter__(self):
        return self

    def __next__(self):
        # this if statement can be true if the given region
        # was empty initially
        if self.contigIndex is None:
            raise StopIteration()
        try:
            # it will not raise StopIteration if we still have reads
            # to iterate over in the current region
            record = self.inputPysamIterator.__next__()
            return record
        except: # we have parsed all the overlapping reads in the current region
            # update coordinates to the next block to iterate over its reads
            contigName, start, end = self._get_next_region()
            if contigName is None: # there is no more block to iterate over
                raise StopIteration()
            self.inputPysamIterator = self.inputPysamAlignmentFile.fetch(contigName, start, end)
            self.__next__()

def getPysamTags(pysamRecord, tagNames):
    """ Return a list of only the requested tags """
    newTagList = []
    for tagName, tagValue in pysamRecord.tags:
        if tagName in tagNames:
            newTagList.append((tagName, tagValue))
    return  newTagList

def getBlockListPerContigFromPysamHeader(pysamHeader):
    """ Returns a dictionary of BlockLists that spans the whole genome """
    blockListPerContig = {}
    for name, length in zip(pysamHeader.references, pysamHeader.lengths):
        blockListPerContig[name] = BlockList([(1, length)])
    return blockListPerContig

def writeSplitAlignments(inputPysamRecord,
                         splitRefBlocks,
                         splitQueryBlocks,
                         splitCigarLists,
                         outputPysamAlignmentFile,
                         tagsToKeep):
    """ Writes split alignments to the output file path. It will use the original pysam record, 
        ref/query coordinates and cigar lists of the split alignments, which are calculated previously 
    """
    splitIndex = 0
    # all input coordindates are 1-based closed 
    for refBlock, queryBlock, cigarList in zip(splitRefBlocks, splitQueryBlocks, splitCigarLists):
        if queryBlock[0] <= queryBlock[1]:
            splitPySamRecord = pysam.AlignedSegment()
            splitPySamRecord.query_name = inputPysamRecord.query_name + "_split_" + str(splitIndex)
            splitPySamRecord.query_sequence = inputPysamRecord.query[queryBlock[0] - 1: queryBlock[1]]
            splitPySamRecord.flag = inputPysamRecord.flag
            splitPySamRecord.reference_id = inputPysamRecord.reference_id
            splitPySamRecord.reference_start = refBlock[0] - 1 # make it 0-based closed
            splitPySamRecord.mapping_quality = inputPysamRecord.mapping_quality
            splitPySamRecord.cigar = createCigarListCompatibleWithPysam(cigarList)
            splitPySamRecord.next_reference_id = inputPysamRecord.next_reference_id
            splitPySamRecord.next_reference_start = inputPysamRecord.next_reference_start
            splitPySamRecord.template_length = inputPysamRecord.template_length
            splitPySamRecord.query_qualities = inputPysamRecord.query_qualities[queryBlock[0] - 1: queryBlock[1]]
            splitPySamRecord.tags = getPysamTags(inputPysamRecord, tagsToKeep)
            # write the record of this split alignment
            outputPysamAlignmentFile.write(splitPySamRecord)
            splitIndex += 1
    return splitIndex

def splitOneAlignment(inputPysamRecord, header, outputPysamAlignmentFile, tagsToKeep, blocks):
    """
        Takes one pysam record and split the alignment into multiple ones using the given blocks
        Returns the number of written split alignments
    """
    alignment = Alignment.createFromPysamRecord(inputPysamRecord, header)
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

def splitAlignmentsInRegionOneShard(args):
    inputSamPath = args[0]
    blockListPerContig = args[1]
    outputSamPath = args[2]
    tagsToKeep = args[3]
    blocks = args[4]
    shared_reads_dict = args[6]
    lock = args[7]
    
    inputPysamAlignmentFile = pysam.AlignmentFile(inputSamPath, "r")
    outputPysamAlignmentFile = pysam.AlignmentFile(outputSamPath, "w", header=inputPysamAlignmentFile.header)
    iterator = PysamIteratorForBlockList(inputPysamAlignmentFile=inputPysamAlignmentFile,
                                         blockListPerContig = blockListPerContig)
    for inputPysamRecord in iterator:
        # if blockListPerContig for different processes are close enough
        # in terms of the coordiantes on the reference
        # we may end up parsing the same read in multiple processes
        # To avoid this we keep track of the parsed reads in a shared
        # dictionary. This dictionary has to be updated once we parsed 
        # a read
        with lock:
            if inputPysamRecord.query_name in shared_reads_dict: 
                continue
            else:
                shared_reads_dict[inputPysamRecord.query_name] = 1
        splitOneAlignment(inputPysamRecord = inputPysamRecord,
                          header=inputPysamAlignmentFile.header,
                          outputPysamAlignmentFile=outputPysamAlignmentFile,
                          tagsToKeep=tagsToKeep,
                          blocks=blocks)
        
    inputPysamAlignmentFile.close()
    outputPysamAlignmentFile.close()

def splitAlignmentsInParallel(inputSamPath, blocks, outputSamPath, tagsToKeep, threads):
    inputPysamAlignmentFile = pysam.AlignmentFile(inputSamPath)
    
    # take the blocks spanning whole genome and split them into regions of roughly equal length
    # the number of regions is equal to the number of threads to make the mount
    # of work that each thread does roughly equal
    wholeGenomeBlockListPerContig = getBlockListPerContigFromPysamHeader(inputPysamAlignmentFile.header)
    regionsForMultiProcessing = BlockList.split(wholeGenomeBlockListPerContig, threads)

    outputSamPathList = []
    with Manager() as manager:
        # create a shared dict for keeping track of the reads
        # that have been parsed so far
        # It is for making sure that we are not outputing the
        # same reads multiple times
        sharedReadsDict = manager.dict()
        lock = manager.Lock()

        pool = Pool(threads)
        for i, regionOneShard in enumerate(regionsForMultiProcessing):
            outputSamPathOneShard = os.path.splitext(outputSamPath)[0] + f'_{i}' + os.path.splitext(outputSamPath)[1]
            outputSamPathList.append(outputSamPathOneShard)
            pool.apply_async(splitAlignmentsInRegionOneShard, args=((inputSamPath, 
                                                                     regionOneShard, 
                                                                     outputSamPathOneShard, 
                                                                     tagsToKeep, 
                                                                     blocks,
                                                                     i, 
                                                                     sharedReadsDict, 
                                                                     lock),))
        pool.close()
        pool.join()
    
    inputPysamAlignmentFile.close()
    return outputSamPathList

def mergeSamFiles(outputSamPath, samPathList, suffix):
    if suffix:
        outputSamPath = os.path.splitext(outputSamPath)[0] + suffix  + os.path.splitext(outputSamPath)[1] 
    mergeParameters = ['-f', outputSamPath] + samPathList
    pysam.merge(*mergeParameters)
    return outputSamPath

def main():
    parser = argparse.ArgumentParser(description='Split the records in BAM file based on the given BED file and create a new BAM file with split records ')
    parser.add_argument('--input', type=str,
                        help='(BAM/SAM format) The alignments to be split')
    parser.add_argument('--blocks', type=str,
                        help='(BED format) The desired blocks in the reference coordinates to be used for splitting alignment records. It should be sorted and with no overlaps.')
    parser.add_argument('--output', type=str,
                        help='(BAM/SAM format) A path for saving the output BAM file with split records.')
    parser.add_argument('--threads', type=int,
                        help='Number of threads')


    # Fetch the arguments
    args = parser.parse_args()
    inputSamPath = args.input
    blocksBed = args.blocks
    outputSamPath = args.output
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

    print(f"[{datetime.datetime.now()}] Input blocks for splitting are parsed.")

    print(f"[{datetime.datetime.now()}] Started splitting alignments using {threads} threads")
    shardSamPathList = splitAlignmentsInParallel(inputSamPath = inputSamPath,
                                           blocks = blocks,
                                           outputSamPath = outputSamPath,
                                           tagsToKeep = ['HP'],
                                           threads = threads)
    print(f"[{datetime.datetime.now()}] Split alignments are written in {len(shardSamPathList)} shard files: {shardSamPathList}")
    
    print(f"[{datetime.datetime.now()}] Started merging files.")
    unsortedSamPath = mergeSamFiles(outputSamPath, shardSamPathList, suffix=".unsorted")
    print(f"[{datetime.datetime.now()}] Merged file (but unsorted): {unsortedSamPath}")

    print(f"[{datetime.datetime.now()}] Started sorting.")
    pysam.sort('-@', '8', '-o', outputSamPath, unsortedSamPath)
    print(f"[{datetime.datetime.now()}] Sorting is done.")

    print(f"[{datetime.datetime.now()}] Started removing temporary sam files.")
    # remove shard sam files
    for samPath in shardSamPathList:
        os.unlink(samPath)
    os.unlink(unsortedSamPath)

    print(f"[{datetime.datetime.now()}] Done! Output SAM file : ", outputSamPath)

if __name__ == "__main__": main()

