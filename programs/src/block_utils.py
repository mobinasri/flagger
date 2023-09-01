import sys
import argparse
from collections import defaultdict
import re
from multiprocessing import Pool
from copy import deepcopy
import random
from Bio import SeqIO

CS_PATTERN = r'(:([0-9]+))|(([+-])([a-z]+)|([\\*]([a-z]+))+)'

def reverseComplement(seq):
    comp={'A':'T',
          'T':'A',
          'C':'G',
          'G':'C',
          'a':'t',
          't':'a',
          'g':'c',
          'c':'g'}
    return "".join([comp[x] for x in seq[::-1]])

def induceSingleBaseErrors(seq, errorRate):
    """
    Takes a sequence and induces single-base errors in it with the given error rate
    :param seq: A string of base letters
    :param errorRate: The ratio of the bases that has be replaced by false bases
    :return: The erroneous sequence
    """
    otherBases = {'A':['G', 'T', 'C'],
                  'T':['A', 'C', 'G'],
                  'C':['A', 'T', 'G'],
                  'G':['C', 'A', 'T'],
                  'a':['t', 'c', 'g'],
                  't':['a', 'c', 'g'],
                  'g':['c', 't', 'a'],
                  'c':['g', 'a', 't']}
    erroneousSeq = []
    for x in seq:
        if random.uniform(0,1) < errorRate:
            erroneousSeq.append(random.choice(otherBases[x]))
        else:
            erroneousSeq.append(x)
    return "".join(erroneousSeq)



def getCigarList(cigarString):
    """
        Returns a list of tuples based on the cigar string. 
        Each tuple has two elements; the first one is showing 
        the operation which could be one of =, X, M, I or D and the 
        second element is an integer which shows the length of 
        the associated operation.
    """
    cigarOps = re.compile("[0-9]+").split(cigarString)[1:]
    # H (hard clipping) or S (Soft clipping) are not included since in the PAF format 
    # the start and end positions are shifted instead of adding H or S to the cigar string
    cigarSizes = [int(size) for size in re.compile("M|I|D|X|=").split(cigarString)[:-1]]
    cigarList = [ (op, size) for op, size in zip(cigarOps, cigarSizes)]
    return cigarList

class BlockList:
    def __init__(self, blocks = None):
        """
        Takes a list of intervals like below
            [(s1, e1), (s2, e2), ... ,(sN, eN)]
        It will copy the list and save the copy in the "blocks" attribute
        A 3rd entry is also added for each tuple, which is initialized to zero
        the purpose of the 3rd entry is mainly for merging and intersecting

        Each interval can also be a triple
            [(s1, e1, c1), (s2, e2, c2), ..., (sN, eN, cN)]
        where c can be an integer data for each interval
        In this case the given 3rd entries will be saved in the "blocks" attribute
        Please note that if the 3rd entries are given
        """
        if blocks == None:
            self.blocks = []
        elif len(blocks) == 0:
            self.blocks = []
        else:
            if len(blocks[0]) == 2:
                self.blocks = sorted([(block[0], block[1], 0) for block in blocks])
            else:
                self.blocks = sorted(blocks)

    def isEqual(self, otherBlockList):
        if len(self.blocks) != len(otherBlockList.blocks): return False
        for int1, int2 in zip(self.blocks, otherBlockList.blocks):
            if int1 != int2: return False
        return True
    def copy(self):
        return BlockList(deepcopy(self.blocks))

    def append(self, block):
        self.blocks.append((block[0], block[1], 0 if len(block) == 2 else block[2]))

    def extend(self, blockList):
        self.blocks.extend(blockList.blocks)

    def setThirdEntry(self, entry):
        for b in self.blocks:
            b[2] = entry

    def sort(self):
        self.blocks.sort()

    def intersect(self, otherBlockList, inplace):
        """
            Note that this function assumes that blocks in self do not have any
            overlap within themselves. Same assumption applies for otherBlockList.
            It will copy the third entries of the self blocks to the intersected blocks
            Arguments:
                otherBlockList: a BlockList to intersect with self blocks
                inplace: If True the intersected blocks will be saved inplace
            Returns:
                If inplace is False it will return a new BlockList with the intersected blocks
        """
        if len(self.blocks) == 0 or len(otherBlockList.blocks) == 0:
            if inplace:
                self.blocks = []
                return 
            else:
                return BlockList([])

        # there is at least one block in both self and otherBlockList
        i1 = 0
        i2 = 0
        s1 = self.blocks[i1][0]
        e1 = self.blocks[i1][1]
        c1 = self.blocks[i1][2]
        s2 = otherBlockList.blocks[i2][0]
        e2 = otherBlockList.blocks[i2][1]
        newBlocks = []

        while True:
            if e2 < s1: # block1 is after block2
                i2 += 1
                if len(otherBlockList.blocks) <= i2:
                    break # break the while loop since there is no more block
                else:
                    # update block2
                    s2 = otherBlockList.blocks[i2][0]
                    e2 = otherBlockList.blocks[i2][1]
                    continue
            elif e1 < s2: # block2 is after block1
                i1 += 1
                if len(self.blocks) <= i1:
                    break # break the while loop
                else:
                    # update block1
                    s1 = self.blocks[i1][0]
                    e1 = self.blocks[i1][1]
                    c1 = self.blocks[i1][2]
                    continue
            else: # there should be an overlap
                newBlocks.append((max(s1, s2), min(e2, e1), c1))
                if e2 < e1:
                    s1 = e2 + 1
                else:
                    i1 += 1
                    if len(self.blocks) <= i1:
                        break # break the while loop
                    else:
                        # update block1
                        s1 = self.blocks[i1][0]
                        e1 = self.blocks[i1][1]
                        c1 = self.blocks[i1][2]
                        continue
        if inplace:
            self.blocks = newBlocks
        else:
            newBlockList = BlockList(newBlocks)
            return newBlockList


    def subtract(self, otherBlockList, inplace):
        """
            Note that this function assumes that blocks in self do not have any
            overlap within themselves. Same assumption applies for otherBlockList.
            Arguments:
                otherBlockList: a BlockList to subtract from self blocks
                inplace: If True the subtracted blocks will be saved inplace
            Returns:
                If inplace is False it will return a new BlockList with the subtracted blocks
        """
        if len(self.blocks) == 0 or len(otherBlockList.blocks) == 0:
            if inplace:
                return
            else:
                return self.copy()
        # there is at least one block in both self and otherBlockList
        i1 = 0
        i2 = 0
        s1 = self.blocks[i1][0]
        e1 = self.blocks[i1][1]
        c1 = self.blocks[i1][2]
        s2 = otherBlockList.blocks[i2][0]
        e2 = otherBlockList.blocks[i2][1]
        newBlocks = []

        while True:
            if e2 < s1: # block1 is after block2
                i2 += 1
                if len(otherBlockList.blocks) <= i2:
                    # add the last block, which could be partial
                    # then add all the remaining blocks completely
                    newBlocks.append((s1, e1, 0))
                    i1 += 1
                    while i1 < len(self.blocks):
                        s1 = self.blocks[i1][0]
                        e1 = self.blocks[i1][1]
                        newBlocks.append((s1, e1, c1))
                        i1 += 1
                    break # break the while loop since there is no more block
                else:
                    # update block2
                    s2 = otherBlockList.blocks[i2][0]
                    e2 = otherBlockList.blocks[i2][1]
                    continue
            elif e1 < s2: # block2 is after block1
                newBlocks.append((s1, e1, c1))
                i1 += 1
                if len(self.blocks) <= i1:
                    break # break the while loop
                else:
                    # update block1
                    s1 = self.blocks[i1][0]
                    e1 = self.blocks[i1][1]
                    continue
            else: # there should be an overlap
                if s1 < s2:
                    newBlocks.append((s1, s2 - 1, c1))
                if e2 < e1:
                    s1 = e2 + 1
                else:
                    i1 += 1
                    if len(self.blocks) <= i1:
                        break # break the while loop
                    else:
                        # update block1
                        s1 = self.blocks[i1][0]
                        e1 = self.blocks[i1][1]
                        c1 = self.blocks[i1][2]
                        continue
        if inplace:
            self.blocks = newBlocks
        else:
            newBlockList = BlockList(newBlocks)
            return newBlockList

    def truncateFromRight(self, lengthToTruncate, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            # skip if the whole block is shorter than given length
            if e - s + 1 <= lengthToTruncate: continue
            newBlocks.append((s, e - lengthToTruncate, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def truncateFromLeft(self, lengthToTruncate, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            # skip if the whole block is shorter than given length
            if e - s + 1 <= lengthToTruncate: continue
            newBlocks.append((s + lengthToTruncate, e, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def extendFromLeft(self, lengthToExtend, leftMostPosition, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            # do not extend before leftMostPosition
            newS = max(s - lengthToExtend, leftMostPosition)
            # skip if block is empty
            if e < newS: continue
            newBlocks.append((newS, e, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def extendFromRight(self, lengthToExtend, rightMostPosition, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            # do not extend after rightMostPosition
            newE = min(e + lengthToExtend, rightMostPosition)
            # skip if block is empty
            if newE < s: continue
            newBlocks.append((s, newE, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def truncateFromBothSides(self, lengthToTruncate, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            # skip if the whole block is shorter than given length
            if e - s + 1 <= (2 * lengthToTruncate): continue
            newBlocks.append((s + lengthToTruncate, e - lengthToTruncate, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def shift(self, lengthToShift, minCoordinate, maxCoordinate, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            newS = s + lengthToShift
            newE = e + lengthToShift
            # skip if the shifted block is out of the range
            if newE < minCoordinate or maxCoordinate < newS: continue
            # truncate if a part of the shifted block is out of the range
            newS = minCoordinate if newS < minCoordinate else newS
            newE = maxCoordinate if maxCoordinate < newE else newE
            newBlocks.append((newS, newE, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def reverse(self, wholeBlockLength, inplace):
        """
            make the coordinates with respect to the end of the whole block/contig
        """
        newBlocks = []
        for s, e, c in self.blocks[::-1]:
            newBlocks.append((wholeBlockLength - e + 1, wholeBlockLength - s + 1, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    def removeBlocksShorterThan(self, minLength, inplace):
        newBlocks = []
        for s, e, c in self.blocks:
            if minLength <= e - s + 1:
                newBlocks.append((s, e, c))
        if inplace:
            self.blocks = newBlocks
        else:
            return BlockList(newBlocks)

    # This function is adapted from ptBlock_merge_blocks_v2() function from Secphase repo v0.4.3
    # https://github.com/mobinasri/secphase/blob/v0.4.3/programs/submodules/ptBlock/ptBlock.c
    def mergeWithOverlapCount(self, inplace):
        blocksMergedFinalized = []
        blocksMergedOngoing = []
        if len(self.blocks) == 0: return blocksMergedFinalized

        for b2 in self.blocks:
            if len(blocksMergedOngoing) == 0: # Initiate bMerged for the first block
                blocksMergedOngoing.append((b2[0], b2[1], 1))
                continue
            e2 = b2[1]
            s2 = b2[0]
            blocksMergedTemp = blocksMergedOngoing
            blocksMergedOngoing = []
            for b1 in blocksMergedTemp:
                e1 = b1[1]
                s1 = b1[0]
                c = b1[2]
                #
                # finalized:
                #
                #  s1       e1
                # [**********]
                #               [----------]
                #                s2       e2
                #
                if e1 < s2 :
                    bMerged = (s1, e1, c)
                    blocksMergedFinalized.append(bMerged)
                elif s1 <= s2:
                    #
                    # finalized:
                    #   s1       e1
                    #  [***-------]
                    #     [----------]
                    #      s2       e2
                    #
                    if s1 < s2: # && s2 <= e1
                        bMerged = (s1, s2-1, c)
                        blocksMergedFinalized.append(bMerged)

                    #
                    #  ongoing:
                    #
                    #    s1       e1                      s1        e1
                    #   [---*******]           OR        [----*****--]
                    #      [*******--]                       [*****]
                    #       s2      e2                        s2   e2
                    #
                    #
                    bMerged = (s2, min(e1, e2), c + 1)
                    blocksMergedOngoing.append(bMerged)
                    #
                    # finalized:
                    #     s1       e1
                    #    [-----*****]
                    #      [---]
                    #       s2 e2
                    #
                    if e2 < e1:
                        bMerged = (e2 + 1, e1, c)
                        blocksMergedOngoing.append(bMerged)
                #
                # ongoing:
                #
                #            s1       e1
                #           [**********]
                #       [----**********---]
                #        s2              e2
                #
                elif e1 <= e2: # && s2 < s1
                    bMerged = (s1, e1, c + 1)
                    blocksMergedOngoing.append(bMerged)
                else: # e2 < e1 && s2 < s1
                    #
                    # ongoing:
                    #
                    #            s1       e1
                    #           [******----]
                    #       [----******]
                    #        s2       e2
                    #
                    if s1 <= e2:
                        bMerged = (s1, e2, c + 1)
                        blocksMergedOngoing.append(bMerged)
                    #
                    # ongoing:
                    #
                    #            s1        e1
                    #           [------****]
                    #       [----------]
                    #        s2       e2
                    #
                    bMerged = (max(e2 + 1, s1), e1, c)
                    blocksMergedOngoing.append(bMerged)
            # add the last non-overlapping block
            #
            # ongoing:
            #
            #       s1        e1
            #      [----------]
            #           [------****]
            #            s2       e2
            #
            if max(e1 + 1, s2) <= e2:
                bMerged = (max(e1 + 1, s2), e2, 1)
                blocksMergedOngoing.append(bMerged)

        # Add the remaining blocks
        for b in blocksMergedOngoing:
            bMerged = (b[0], b[1], b[2])
            blocksMergedFinalized.append(bMerged)

        if inplace:
            self.blocks = blocksMergedFinalized
        else:
            return BlockList(blocksMergedFinalized)

    @staticmethod
    def parseBed(bedPath, saveFourthColumnAsNumeric):
        """
        Parse a bed file
        :param bedPath: The path to a bed file
        :param saveFourthColumnAsInt: If each track has a 4th column parse it and save it as the 3rd entry of each block in
                                      BlockList
        :return: A dictionary that contains contig names as keys and one BlockList as each value. Each BlockList will
                 contain all the blocks for one contig
        """


        # create a dictionary whose values are BlockLists
        blockListPerContig = defaultdict(BlockList)

        # parse tracks and save them in the dictionary
        with open(bedPath, "r") as f:
            for line in f:
                cols = line.strip().split()
                if len(cols) == 4 and saveFourthColumnAsNumeric == True:
                    c = float(cols[3])
                else:
                    c = 0
                start = int(cols[1])
                end = int(cols[2])
                contig = cols[0]
                blockListPerContig[contig].append((start, end, c))
        # sort blocks per contig
        for contig, blockList in blockListPerContig.items():
            blockList.sort()
        return blockListPerContig


def parseFasta(fastaPath):
    """
    Parse fasta file and save the sequences in a dictionary
    """
    contigSequences = {}
    fasta_sequences = SeqIO.parse(open(fastaPath),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        contigSequences[name] = sequence
    return  contigSequences


class Alignment:
    """
        A class for saving alignment data
        The constructor receives a line from a file with PAF format
        (More about PAF format https://github.com/lh3/miniasm/blob/master/PAF.md)
        and parse the info about the alignment such as the coordinates and the cigar operations.
    """

    def __init__(self, paf_line):

        cols = paf_line.strip().split()
        self.contigName = cols[0]
        self.contigLength = int(cols[1])
        self.contigStart = int(cols[2]) # 0-based closed
        self.contigEnd = int(cols[3]) # 0-based open
        self.orientation = cols[4]
        self.chromName = cols[5]
        self.chromLength = int(cols[6])
        self.chromStart = int(cols[7]) # 0-based closed
        self.chromEnd = int(cols[8]) # 0-based open
        self.isPrimary = False
        if "tp:A:P" in paf_line:
            self.isPrimary = True
        # The cigar string starts after "cg:Z:"
        afterCg = paf_line.strip().split("cg:Z:")[1]
        cigarString = afterCg.split()[0]
        self.cigarList = getCigarList(cigarString)
        if "NM:i:" in paf_line:
            # The edit distance starts after "NM:i:"
            afterNM = paf_line.strip().split("NM:i:")[1]
            editString = afterNM.split()[0]
            self.editDistance = int(editString)
        else:
            self.editDistance = None


def reverseInterval(interval, contigLength):
    """
        Returns the reversed coordinates of the given interval 
        (In other words, what would be the coordinates if we 
        started from the end of the contig).
    """
    #print(interval[0], interval[1])
    assert (interval[0] <= interval[1])
    return (contigLength - interval[1] + 1, contigLength - interval[0] + 1)

def reverseBlocks(blocks, contigLength):
    """
        Return the reversed order and coordinates of the given blocks
        It is assumed that the blocks are given sorted.
    """
    reverseBlocks = []
    for i in range(len(blocks)-1, -1, -1):
        reversedInterval = reverseInterval(blocks[i], contigLength)
        info = blocks[i][2]
        reverseBlocks.append((reversedInterval[0], reversedInterval[1], info))
    return reverseBlocks

def convertIndelsInCigar(cigarList):
    """
        Whenever we want to have the cigar operations with respect to the query sequences instead of the reference
        this function can be used for replacing the insertions by deletions and vice versa.
    """
    newCigarList = []
    for cigarOp, cigarSize in cigarList:
        if cigarOp == 'D':
            newCigarList.append(('I', cigarSize))
        elif cigarOp == 'I':
            newCigarList.append(('D', cigarSize))
        else: # 'M'/'='/'X'
            newCigarList.append((cigarOp, cigarSize))
    return newCigarList



class Projection:
    """
        A class for saving the data related to the projection of a list of blocks
    """

    def __init__(self, orientation, contigLength, isCigarConverted, isCigarReversed):
        self.orientation = orientation
        self.contigLength = contigLength
        self.isCigarConverted = isCigarConverted
        self.isCigarReversed = isCigarReversed
        self.projectionBlocks = []
        self.projectableBlocks = []
        self.projectionCigarList = []

    def update(self,
               projectionStartPos, projectionEndPos,
               projectableStartPos, projectableEndPos,
               info, diff, projectionCigar):
        r = None if diff == None else diff/ (projectionEndPos - projectionStartPos + 1) * 100
        self.projectionBlocks.append((projectionStartPos, projectionEndPos, info, r)) # there is no valid projection for this block
        if self.orientation == '+':
            self.projectableBlocks.append((projectableStartPos, projectableEndPos, info, r))
        else:
            reversedInterval = reverseInterval((projectableStartPos, projectableEndPos), self.contigLength)
            self.projectableBlocks.append((reversedInterval[0], reversedInterval[1], info, r))
        if self.isCigarConverted:
            convertedCigar = convertIndelsInCigar(projectionCigar)
            if self.isCigarReversed:
                self.projectionCigarList.append(convertedCigar[::-1])
            else:
                self.projectionCigarList.append(convertedCigar)
        else:
            self.projectionCigarList.append(projectionCigar)



def findProjectionsInternal(mode, cigarList, forwardBlocks,
                    chromLength, chromStart, chromEnd,
                    contigLength, contigStart, contigEnd,
                    orientation, includeEndingIndel, includePostIndel, isCigarConverted, isCigarReversed):
    """
        Returns:
            * projectableBlocks: A list of tuples which contains the blocks that could be projected 
                                    (onto the reference in 'asm2ref' mode or onto the assembly contig in 'ref2asm' mode)
            * projectionBlocks:  A list of tuples which contains the projections of the projectable blocks 
        Arguments:
            * mode: (str) can be either 'ref2asm' or 'asm2ref'. The function has been written initially for the mode 'asm2ref'
                    The algorithm for the mode 'ref2asm' is pretty similar and we can use the same function for the second
                    one by just swapping and adjusting some of the given arguments. (Read the first comment for more information)
                        'ref2asm':
                                    In this mode the blocks are given in the coordinates of the reference and 
                                    the output will be the projections of those blocks onto the assembly
                        'asm2ref':
                                    In this mode the blocks are given in the coordinates of the assembly and
                                    the output will be the projections of those blocks onto the reference
            * cigarList: A list of tuples which contains the cigar operations and their lengths 
                        [(op1, length1) , (op2, length2), ...]
            * forwardBlocks: A list of tuples which contains the blocks in the coordinates of the contig. Each block is a tuple of (start, end, info). start and end should be 1-based closed, info is optional and can be empty
            * chromLength: (Integer) The total length of the reference chromosome
            * chromStart: (Integer) Where the alignment begins in the coordinates of the reference (should be 1-based and closed)
            * contigLength: (Integer) The total length of the contig
            * contigStart: (Integer) Where the alignment begins in the coordinates of the query contig (should be 1-based and closed)
            * contigEnd: (Integer) Where the alignment ends in the coordinates of the query contig (1-based and closed)
            * orientation: (str) {'+', '-'}
    """
    # If the mode is 'ref2asm' we can simply call the same function but by changing the arguments in a clever manner!:
    #   1. The cigar operations should become with respect to the assembly contig not the reference
    #      to reach that aim we should call "convertIndelsInCigar"
    #   2. If the orientation is negative we should reverse the order of the cigar operations
    #   3. After this conversion we can assume that our reference is the assembly contig so the 
    #      the length, start and end positions of the reference and the assembly contig should be swapped
    #   4. Finally all we need is just to call the function again but in the mode of 'asm2ref'
    if mode == 'ref2asm':
        convertedCigarList = convertIndelsInCigar(cigarList)
        if orientation == '+':
            return findProjectionsInternal('asm2ref', convertedCigarList, forwardBlocks,
                                           contigLength, contigStart, contigEnd,
                                           chromLength, chromStart, chromEnd, 
                                           orientation, includeEndingIndel, includePostIndel, True, False)
        else:
            return findProjectionsInternal('asm2ref', convertedCigarList[::-1], forwardBlocks,
                                           contigLength, contigStart, contigEnd,
                                           chromLength, chromStart, chromEnd,
                                           orientation, includeEndingIndel, includePostIndel, True, True)
    projection = Projection(orientation, contigLength, isCigarConverted, isCigarReversed)
    blockIdx = 0
    nextOpStartRef = chromStart
    nextOpStartContig = None
    currOpStartRef = None
    currOpStartContig = None
    # Return if the blocks has no overlap with the alignment
    if (forwardBlocks[-1][1] < contigStart) or (contigEnd < forwardBlocks[0][0]):
        return projection.projectableBlocks, projection.projectionBlocks, projection.projectionCigarList
    # The cigar starts from the end of the contig if the alignment orientation is negative,
    # so the blocks coordinates and their order should be reversed in that case.
    # (Note that the blocks will be reversed back after the projections are all found) 
    if orientation == '+':
        blocks = forwardBlocks
        nextOpStartContig = contigStart
    else:
        blocks = reverseBlocks(forwardBlocks, contigLength)
        nextOpStartContig = contigLength - contigEnd + 1
    # find the first block not completely clipped and set the blockIdx accordingly
    while (blockIdx < len(blocks)) and (blocks[blockIdx][1] < nextOpStartContig):
        blockIdx += 1
    # If a block starts from before where the whole alignment starts, the projection start point 
    # is where the first operation occurs (The first operation is always M)
    if blocks[blockIdx][0] < nextOpStartContig:
        projectionStartPos = nextOpStartRef
        projectableStartPos = nextOpStartContig

    diff = 0
    projectionCigar = []
    blockEndedOnTheEdgeOfOperation = False
    preCigarOp = None
    preCigarSize = None
    # iterate over cigar elements and find the projections
    for cigarOp, cigarSize in cigarList:
        #print(cigarOp,cigarSize)
        ####################################
        #### Case 1: Mismatch or Match #####
        ####################################
        if cigarOp == 'M' or cigarOp == 'X' or cigarOp == '=':
            currOpStartContig = nextOpStartContig
            currOpStartRef = nextOpStartRef
            nextOpStartContig = currOpStartContig + cigarSize
            nextOpStartRef = currOpStartRef + cigarSize
            # save the previous block if it ended on the (right-most) edge of the previous operation
            if blockEndedOnTheEdgeOfOperation:
                projection.update(projectionStartPos, projectionEndPos,
                                  projectableStartPos, projectableEndPos,
                                  blocks[blockIdx][2], diff, projectionCigar)
                projectionCigar = []
                diff = 0
                blockIdx += 1
                if len(blocks) <= blockIdx: break
                # reset blockEndedOnTheEdgeOfOperation
                blockEndedOnTheEdgeOfOperation = False
            # if the whole operation was within the block
            ###
            # REF:  AAAAAAAAAAAAAAAAAAA
            # ASM:  AAAAAAA^^^^AAAAAAAA
            #              ||||
            # BLK:      [*********]
            ###
            if blocks[blockIdx][0] < currOpStartContig and nextOpStartContig <= blocks[blockIdx][1]:
                # append the whole operation
                projectionCigar.append((cigarOp, cigarSize))
                if cigarOp == 'X': diff += cigarSize
            # When the while loop ends the blockIdx points to the block whose end position is 
            # not within the current operation (each operation can be assumed as an interval)
            # There exists three scenarios for the start position of the block with blockIdx:
            #
            # (1) the start position is before the current operation (The current operation is *M)
            #   operations : -----------[            *M          ][        I       ][   M    ]-----
            #   blocks :     -------[             Block                ]---------------------------
            #
            # (2) the start position is within the current operation
            #   operations : -----------[            *M          ][       I       ][      M   ]-----
            #   blocks :     ------------------[           Block      ]-----------------------------
            #
            # (3) the start position is after the current operation
            #   operations : -----------[            *M          ][        I      ][      M   ]-----
            #   blocks :     ----------------------------------------[       Block      ]-----------
            #
            while (blockIdx < len(blocks)) and (blocks[blockIdx][1] < nextOpStartContig):
                # if the whole block is within the operation
                # update the start positions
                ###
                # REF:  AAAAAAAAAAAAAAAAAAA
                # ASM:  AAAAA^^^^^^^^^^AAAAAA
                #              ||||||
                # BLK:         [****]
                ###
                if currOpStartContig <= blocks[blockIdx][0]:
                    projectionStartPos = currOpStartRef + blocks[blockIdx][0] - currOpStartContig
                    projectableStartPos = blocks[blockIdx][0]
                    if currOpStartContig == blocks[blockIdx][0] and \
                        includePostIndel and \
                        isCigarReversed and \
                        preCigarOp == 'D':
                        projectionStartPos = currOpStartRef - preCigarSize
                        projectionCigar.insert(0, (preCigarOp, preCigarSize))
                # otherwise only the end position in withtin the operation
                # In either case the end positions should be updated
                ###
                # REF:  AAAAAAAAAAAAAAAAAAAAA
                # ASM:  AAAAA^^^^^^^^^^AAAAAA
                #            ||||||
                # BLK:   [   *****]
                ###
                projectionEndPos = currOpStartRef + blocks[blockIdx][1] - currOpStartContig
                projectableEndPos = blocks[blockIdx][1]
                # find the size of the overlap between the current operation and the block
                overlapOpSize = blocks[blockIdx][1] - max(blocks[blockIdx][0], currOpStartContig) + 1
                # append the overlapped operation
                projectionCigar.append((cigarOp, overlapOpSize))
                if cigarOp == 'X': diff += overlapOpSize
                # if the current block is ending exactly at the end of this operation
                # and we should keep the post-ending indel
                # (to be more specific post-ending deletion in 'asm2ref' mode)
                if blocks[blockIdx][1] == nextOpStartContig - 1 and includePostIndel and not isCigarReversed:
                    blockEndedOnTheEdgeOfOperation = True
                    break
                else: # otherwise we can easily finish the projection process of this block
                    projection.update(projectionStartPos, projectionEndPos,
                                      projectableStartPos, projectableEndPos,
                                      blocks[blockIdx][2], diff, projectionCigar)
                projectionCigar = []
                diff = 0
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            # In case of the 2nd scenario mentioned above, the projection start position should be updated
            ###
            # REF:  AAAAAAAAAAAAAAAAAAAAAA
            # ASM:  AAAAA^^^^^^^^^^AAAAAAA
            #                 |||||
            # BLK:            [****    ]
            ###
            if (currOpStartContig <= blocks[blockIdx][0]) and (blocks[blockIdx][0] < nextOpStartContig) and (nextOpStartContig <= blocks[blockIdx][1]):
                # find the overlap size between the current operation and the block
                overlapOpSize = (nextOpStartContig - 1) - blocks[blockIdx][0] + 1
                # append the overlapped operation
                projectionCigar.append((cigarOp, overlapOpSize))
                if cigarOp == 'X': diff += overlapOpSize
                projectionStartPos = currOpStartRef + blocks[blockIdx][0] - currOpStartContig
                projectableStartPos = blocks[blockIdx][0]
                if currOpStartContig == blocks[blockIdx][0] and \
                        includePostIndel and \
                        isCigarReversed and \
                        preCigarOp == 'D':
                    projectionStartPos = currOpStartRef - preCigarSize
                    projectionCigar.insert(0, (preCigarOp, preCigarSize))
        ####################################
        ####### Case 2: Insertion ##########
        ####################################
        elif cigarOp == 'I':
            currOpStartContig = nextOpStartContig
            nextOpStartContig = currOpStartContig + cigarSize
            # For insertion there is no need to shift the reference
            currOpStartRef = nextOpStartRef

            # save the previous block if it ended on the (right-most) edge of the previous operation
            if blockEndedOnTheEdgeOfOperation:
                projection.update(projectionStartPos, projectionEndPos,
                                  projectableStartPos, projectableEndPos,
                                  blocks[blockIdx][2], diff, projectionCigar)
                projectionCigar = []
                diff = 0
                blockIdx += 1
                if len(blocks) <= blockIdx: break
                # reset blockEndedOnTheEdgeOfOperation
                blockEndedOnTheEdgeOfOperation = False

            # if whole operation is within block with the last base not ending
            ###
            # REF:  AAAAAAA----AAAAAAAA
            # ASM:  AAAAAAAAAAAAAAAAAAA
            #              ||||
            # BLK:     [   ****    ]
            ###
            if blocks[blockIdx][0] < currOpStartContig and nextOpStartContig <= blocks[blockIdx][1]:
                projectionCigar.append((cigarOp, cigarSize))
                diff += cigarSize
            # The beginning and endings of the blocks are excluded from the projection if they are insertions
            # So the blocks are trimmed to make their projections start and end in match/mismatch
            while (blockIdx < len(blocks)) and (blocks[blockIdx][1] < nextOpStartContig):
                # If a block is completely within an insertion then there is no valid projection for it
                ###
                # REF:  AAAAAAA---------AAAAAAAA
                # ASM:  AAAAAAAAAAAAAAAAAAAAAAAA
                #                  |||||
                # BLK:             [***]
                ###
                if currOpStartContig <= blocks[blockIdx][0]:
                    projectionCigar.append((cigarOp, cigarSize))
                    projectionStartPos = None
                    projectionEndPos = None
                    if currOpStartContig == blocks[blockIdx][0] and \
                            includePostIndel and \
                            isCigarReversed and \
                            preCigarOp == 'D':
                        projectionStartPos = currOpStartRef - preCigarSize
                        projectionEndPos = currOpStartRef - 1
                        projectionCigar.insert(0, (preCigarOp, preCigarSize))
                    projectableStartPos = blocks[blockIdx][0]
                    projectableEndPos = blocks[blockIdx][1]
                    diff = None # no projection so divergence not defined
                    # if the block ended exactly at the end of this operation
                    # then we may have deletion for the next operation
                    # so we should set the blockEndedOnTheEdge flag to true
                    # then check this flag in the next operation
                    if blocks[blockIdx][1] == nextOpStartContig - 1 and includePostIndel and not isCigarReversed:
                        blockEndedOnTheEdgeOfOperation = True
                        break # break here to go to the next cigar and check if it's deletion
                    else:
                        if includeEndingIndel:
                            projection.update(projectionStartPos, projectionEndPos,
                                              projectableStartPos, projectableEndPos,
                                              blocks[blockIdx][2], diff, projectionCigar)
                        projectionCigar = []
                        # go to the next block
                        blockIdx += 1
                        diff = 0
                        continue
                # There is a block that ends in an insertion
                ###
                # REF:  AAAAAAA---------AAAAAAAA
                # ASM:  AAAAAAAAAAAAAAAAAAAAAAAA
                #              |||||
                # BLK:     [   ****]
                ###
                # (Note that the ending inserted part was not projected in v0.3.2)
                # The end position of the projection block will be one base before where the insertion starts
                # Note that one base before the insertion is absolutely an M
                projectionEndPos = currOpStartRef - 1
                projectableEndPos = blocks[blockIdx][1] if includeEndingIndel else currOpStartContig - 1
                # append the overlapped operation if overlapSize was positive
                # (or if equivalently includeEndingIndel was true)
                overlapOpSize = projectableEndPos - currOpStartContig + 1
                if overlapOpSize > 0: 
                    projectionCigar.append((cigarOp, overlapOpSize))
                    diff += overlapOpSize
                projection.update(projectionStartPos, projectionEndPos,
                                  projectableStartPos, projectableEndPos,
                                  blocks[blockIdx][2], diff, projectionCigar)
                projectionCigar = []
                diff = 0
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            # If the last block starts with insertion but is not completely an insertion,
            # the initial inserted part is not projectable onto the reference so we start 
            # from the next operation (which should be an M)
            ###
            # REF:  AAAAAAA---------AAAAAAAA
            # ASM:  AAAAAAAAAAAAAAAAAAAAAAAA
            #                   ||||
            # BLK:              [***    ]
            ###
            if (currOpStartContig <= blocks[blockIdx][0]) and (blocks[blockIdx][0] < nextOpStartContig) and (nextOpStartContig <= blocks[blockIdx][1]) :
                if currOpStartContig == blocks[blockIdx][0] and \
                        includePostIndel and \
                        isCigarReversed and \
                        preCigarOp == 'D':
                    projectionStartPos = currOpStartRef - preCigarSize
                    projectionCigar.insert(0, (preCigarOp, preCigarSize))
                else:
                    projectionStartPos = currOpStartRef
                projectableStartPos = blocks[blockIdx][0] if includeEndingIndel else nextOpStartContig
                # append the overlapped operation if overlapSize was positive
                # (or if equivalently includeEndingIndel was true)
                overlapOpSize = (nextOpStartContig - 1) - projectableStartPos + 1
                if overlapOpSize > 0: 
                    projectionCigar.append((cigarOp, overlapOpSize))
                    diff += overlapOpSize

        ####################################
        ####### Case 2: Deletion ###########
        ####################################
        elif cigarOp == 'D':
            # if deletion is completely within the block
            ###
            # REF:  AAAAAAAAAAAAAAAAA
            # ASM:  AAAAAA----AAAAAAA
            #             ||||
            # BLK:   [    ****   ]
            ###
            if blocks[blockIdx][0] < nextOpStartContig and nextOpStartContig <= blocks[blockIdx][1]:
                # append the whole deletion
                projectionCigar.append((cigarOp, cigarSize))
                diff += cigarSize
            # if the block ended exactly one base before the deletion
            ###
            # REF:  AAAAAAAAAAAAAAAAA
            # ASM:  AAAAAA----AAAAAAA
            #             ||||
            # BLK:   [   ]****
            ###
            elif blockEndedOnTheEdgeOfOperation:
                # append the whole deletion
                projectionCigar.append((cigarOp, cigarSize))
                diff += cigarSize
                projectionEndPos += cigarSize
                projection.update(projectionStartPos, projectionEndPos,
                                  projectableStartPos, projectableEndPos,
                                  blocks[blockIdx][2], diff, projectionCigar)
                blockEndedOnTheEdgeOfOperation = False
                projectionCigar = []
                diff = 0
                blockIdx += 1
                if blockIdx >= len(blocks): break
            nextOpStartRef += cigarSize
        preCigarOp = cigarOp
        preCigarSize = cigarSize

    # Note that nextOpStart(Contig | Ref) are not pointing to any cigar operation at this moment
    # since iterating over the cigar elements has finished. Those variables are now pointing to one base
    # after the end of the last cigar operation.
    # Handle the last block that has started but not finished by the end of the last cigar operation
    if (blockIdx < len(blocks)) and (blocks[blockIdx][0] < nextOpStartContig):
        projectionEndPos = nextOpStartRef - 1
        projectableEndPos = nextOpStartContig - 1
        projection.update(projectionStartPos, projectionEndPos,
                          projectableStartPos, projectableEndPos,
                          blocks[blockIdx][2], diff, projectionCigar)

    return projection.projectableBlocks, projection.projectionBlocks, projection.projectionCigarList

def findProjections(mode, cigarList, forwardBlocks,
                    chromLength, chromStart, chromEnd,
                    contigLength, contigStart, contigEnd,
                    orientation, includeEndingIndel, includePostIndel):
    return findProjectionsInternal(mode, cigarList, forwardBlocks,
                    chromLength, chromStart, chromEnd,
                    contigLength, contigStart, contigEnd,
                    orientation, includeEndingIndel, includePostIndel, False, False)


def iterateCS(cs_str):
    pattern = re.compile(CS_PATTERN)
    for m in re.finditer(pattern, cs_str):
        op_0 = m.group()[0]
        op_1 = m.group()[1:]
        if op_0 == ':':
            op_len = int(op_1)
        elif op_0 == '*':
            op_len = int((len(op_1) + 1) / 3 )
        else:
            op_len = len(op_1)
        if op_0 == ':':
            op = 'M'
        elif op_0 == '*':
            op = 'X'
        elif op_0 == '+':
            op = 'I'
        elif op_0 == '-':
            op = 'D'
        yield (op, op_len)



def iterateCigar(alignment):
    """
        A generator functions which recieves an alignment and iterates over its cigar operations
        Returns:
            op: The cigar operation could be "M", "X", "=", "I" or "D"
            opLen: The length of the operation
            refInterval:  A tuple of (start, end) which shows the interval of the operation on reference 
            contigInterval: A tuple of (start, end) which shows the interval of the operation on assembly
            Note that start is 0-based closed and end is 0-based open
    """
    refOpStart = alignment.chromStart
    refOpEnd = alignment.chromStart
    orientation = alignment.orientation
    if alignment.orientation == '+':
        contigOpStart = alignment.contigStart
        contigOpEnd = alignment.contigStart
    else: # orientation '-'
        contigOpStart = alignment.contigEnd
        contigOpEnd = alignment.contigEnd
    for op, opLen in alignment.cigarList:
        if op == 'D':
            refOpStart = refOpEnd
            refOpEnd = refOpStart + opLen
            if orientation == '+':
                contigOpStart = contigOpEnd
            else:
                contigOpEnd = contigOpStart
        elif op == 'I':
            refOpStart = refOpEnd
            if orientation == '+':
                contigOpStart = contigOpEnd
                contigOpEnd = contigOpStart + opLen
            else:
                contigOpEnd = contigOpStart
                contigOpStart = contigOpEnd - opLen
        else: # 'M', 'X' or '='
            refOpStart = refOpEnd
            refOpEnd = refOpStart + opLen
            if orientation == '+':
                contigOpStart = contigOpEnd
                contigOpEnd = contigOpStart + opLen
            else:
                contigOpEnd = contigOpStart
                contigOpStart = contigOpEnd - opLen
        refInterval = (refOpStart, refOpEnd)
        contigInterval = (contigOpStart, contigOpEnd)
        yield op, opLen, refInterval, contigInterval
        
        
def getLongInsertionBlocks(alignment, threshold):
    """
        Arguments:
            alignment: An Alignment object containing the alignment info
            threshold: Any insertion longer than this threshold would be returned
        Returns:
            A list of blocks showing the start and end coordinates of long insertions
            in the corresponding contig
            Each block is a tuple of (contig name, start, end)
            Note that start is 0-based closed and end is 0-based open
    """
    blocks = []
    for op, opLen, refInterval, contigInterval in iterateCigar(alignment):
        if op == 'I' and opLen > threshold:
            blocks.append((alignment.contigName, contigInterval[0], contigInterval[1]))
    return blocks
    
def getLongDeletionBlocks(alignment, threshold):
    """
        Arguments:
            alignment: An Alignment object containing the alignment info
            threshold: Any deletion longer than this threshold would be returned
        Returns:
            A list of blocks showing the start and end coordinates of long deletion
            in the corresponding reference chromosome
            Each block is a tuple of (chrom name, start, end)
            Note that start is 0-based closed and end is 0-based open
    """
    blocks = []
    for op, opLen, refInterval, contigInterval in iterateCigar(alignment):
        if op == 'D' and opLen > threshold:
            blocks.append((alignment.chromName, refInterval[0], refInterval[1]))
    return blocks

def parseAssemblyIntervals(faiPath):
    """
        Arguments:
           faiPath: The path to a fasta index file
        Returns:
            A dictionary which has the names of contigs as keys 
            and each value is a list of intervals initialised to [[0, length of contig]]
            encompassing the whole contig
    """
    assemblyIntervals = {}
    with open(faiPath) as f:
        for line in f:
            cols = line.strip().split()
            contigName = cols[0]
            contigLength = int(cols[1])
            assemblyIntervals[contigName] = [[0,contigLength]]
    return assemblyIntervals



def makeCigarString(cigarList):
    cigarFlattened = ["".join((str(opSize), op)) for op, opSize in cigarList]
    return "".join(cigarFlattened)

def runProjection(alignment, mode, blocks, includeEndingIndel, includePostIndel):
    # Extract the alignment attributes like the contig name, alignment boundaries, orientation and cigar
    chromName = alignment.chromName
    contigName = alignment.contigName
    orientation = alignment.orientation
    if alignment.isPrimary == False:
        return [chromName, contigName, orientation, [], []]
    # rBlocks contains the projections and
    # qBlocks contains the projectable blocks
    if mode == "asm2ref":
        if len(blocks[contigName]) == 0: # Continue if there is no block in the contig
            return [chromName, contigName, orientation, [], []]
            #print(blocks[contigName], contigStart, contigEnd, chrom, chromStart, chromEnd)
        projectableBlocks, projectionBlocks, cigarList = findProjections(mode,
                                                                         alignment.cigarList,
                                                                         blocks[contigName],
                                                                         alignment.chromLength,
                                                                         alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                                         alignment.contigLength,
                                                                         alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                                         alignment.orientation,
                                                                         includeEndingIndel, includePostIndel)
    else:
        if len(blocks[chromName]) == 0: # Continue if there is no block in the chrom
            return [chromName, contigName, orientation, [], []]
        projectableBlocks, projectionBlocks, cigarList = findProjections(mode,
                                                                         alignment.cigarList,
                                                                         blocks[chromName],
                                                                         alignment.chromLength,
                                                                         alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                                         alignment.contigLength,
                                                                         alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                                         alignment.orientation,
                                                                         includeEndingIndel, includePostIndel)
    return [chromName, contigName, orientation, projectableBlocks, projectionBlocks, cigarList]

def runProjectionParallel(alignments, mode, blocks, includeEndingIndel, includePostIndel, threads):
    pool = Pool(threads)
    #print("Started projecting")
    results = pool.starmap(runProjection, [(alignment, mode, blocks,includeEndingIndel, includePostIndel) for alignment in alignments])
    pool.close()
    return results


def mergeBlockListsPerContigWithOverlapCount(blockListsPerContig):
    mergedBlockListsPerContig = defaultdict(BlockList)
    for contigName, blockList in blockListsPerContig.items():
        mergedBlockListsPerContig[contigName] = blockList.mergeWithOverlapCount(inplace=False)
    return mergedBlockListsPerContig

def getBlockListsPerRefContig(alignments):
    blockListsPerRefContig = defaultdict(BlockList)
    for alignment in alignments:
        # make coors 1-based
        blockListsPerRefContig[alignment.chromName].append((alignment.chromStart + 1, alignment.chromEnd))

    for contig in blockListsPerRefContig:
        blockListsPerRefContig[contig].sort()

    return blockListsPerRefContig

def getBlockListsPerQueryContig(alignments):
    blockListsPerQueryContig = defaultdict(BlockList)
    for alignment in alignments:
        # make coors 1-based
        blockListsPerQueryContig[alignment.contigName].append((alignment.contigStart + 1, alignment.contigEnd))

    for contig in blockListsPerQueryContig:
        blockListsPerQueryContig[contig].sort()

    return blockListsPerQueryContig

def getBlockListsWithSingleAlignmentPerRefContig(alignments):
    """
    Takes a list of alignments and returns a dictionary that contains the blocks with single alignment in the reference coordinates.
    The dictionary has ref contig names as keys and one BlockList as each value.
    :param alignments: A list of alignments (no need to be sorted)
    """
    blockListsPerRefContig = getBlockListsPerRefContig(alignments)
    mergedBlockListsPerRefContig = mergeBlockListsPerContigWithOverlapCount(blockListsPerRefContig)

    blockListsWithSingleAlignmentPerRefContig = defaultdict(BlockList)
    for refContig, mergedBlockList in mergedBlockListsPerRefContig.items():
        for start, end, count in mergedBlockList.blocks:
            if count  == 1:
                blockListsWithSingleAlignmentPerRefContig[refContig].append((start, end))
    return blockListsWithSingleAlignmentPerRefContig

def getBlockListsWithSingleAlignmentPerQueryContig(alignments):
    """
    Takes a list of alignments and returns a dictionary that contains the blocks with single alignment in the query coordinates.
    The dictionary has query contig names as keys and one BlockList as each value.
    :param alignments: A list of alignments (no need to be sorted)
    """
    blockListsPerQueryContig = getBlockListsPerQueryContig(alignments)
    mergedBlockListsPerQueryContig = mergeBlockListsPerContigWithOverlapCount(blockListsPerQueryContig)

    blockListsWithSingleAlignmentPerQueryContig = defaultdict(BlockList)
    for refContig, mergedBlockList in mergedBlockListsPerQueryContig.items():
        for start, end, count in mergedBlockList.blocks:
            if count  == 1:
                blockListsWithSingleAlignmentPerQueryContig[refContig].append((start, end))
    return blockListsWithSingleAlignmentPerQueryContig

def subsetAlignmentsToRefBlocks(alignments, blockListsPerRefContig):
    """
    Extracts subsets of the given alignments by projecting the blocks given in ref coordinates
    :param alignments: A list of alignments (no need to be sorted)
    :param blockListsPerRefContig: A dictionary that contains blocks in the reference coordinates.
    The dictionary should have ref contig names as keys and one BlockList as each value.
    :return: The list of subset alignments that does not span any bases outside the given blocks in the ref coordinates
             The output is sorted by ref coordinates
    """
    threads = 8
    blocksPerRefContig = {}
    for ctg, blockList in blockListsPerRefContig.items():
        blocksPerRefContig[ctg] = blockList.blocks
    results = runProjectionParallel(alignments, 'ref2asm', blocksPerRefContig, False, False, threads)

    contigLengths = {}
    for alignment in alignments:
        contigLengths[alignment.chromName] = alignment.chromLength
        contigLengths[alignment.contigName] = alignment.contigLength

    subsetAlignments = []
    for res in results:
        rContig = res[0]
        qContig = res[1]
        orientation = res[2]
        projectableBlocks = res[3]
        projectionBlocks = res[4]
        cigarList = res[5]
        if len(projectionBlocks) == 0 or len(projectableBlocks) == 0:
            continue
        # projection blocks are in query coordinates
        # projectable blocks are in ref coordinates
        for rBlock, qBlock, cigar in zip(projectableBlocks, projectionBlocks, cigarList):
            # make a paf line
            pafLine = f"{qContig}\t{contigLengths[qContig]}\t{qBlock[0]-1}\t{qBlock[1]}"
            pafLine += f"\t{orientation}"
            pafLine += f"\t{rContig}\t{contigLengths[rContig]}\t{rBlock[0]-1}\t{rBlock[1]}"
            pafLine += "\ttp:A:P" # it is assumed that only primary alignments are used
            pafLine += f"\tcg:Z:{makeCigarString(cigar)}"

            # make a new alignment based on the current projection
            alignment = Alignment(pafLine)
            subsetAlignments.append(alignment)
    subsetAlignments.sort(key = lambda x : (x.chromName, x.chromStart, x.chromEnd, x.contigName, x.contigStart, x.contigEnd))
    return  subsetAlignments

def subsetAlignmentsToQueryBlocks(alignments, blockListsPerQueryContig):
    """
    Extracts subsets of the given alignments by projecting the blocks given in the query coordinates
    :param alignments: A list of alignments (no need to be sorted)
    :param blockListsPerRefContig: A dictionary that contains blocks in the query coordinates.
    The dictionary should have query contig names as keys and one BlockList as each value.
    :return: The list of subset alignments that does not span any bases outside the given blocks in the query coordinates
             The output is sorted by query coordinates
    """
    threads = 8
    blocksPerQueryContig = {}
    for ctg, blockList in blockListsPerQueryContig.items():
        blocksPerQueryContig[ctg] = blockList.blocks
    results = runProjectionParallel(alignments, 'asm2ref', blocksPerQueryContig, False, False, threads)

    contigLengths = {}
    for alignment in alignments:
        contigLengths[alignment.chromName] = alignment.chromLength
        contigLengths[alignment.contigName] = alignment.contigLength

    subsetAlignments = []
    for res in results:
        rContig = res[0]
        qContig = res[1]
        orientation = res[2]
        projectableBlocks = res[3]
        projectionBlocks = res[4]
        cigarList = res[5]
        if len(projectionBlocks) == 0 or len(projectableBlocks) == 0:
            continue
        # projectable blocks are in query coordinates
        # projection blocks are in ref coordinates
        for qBlock, rBlock, cigar in zip(projectableBlocks, projectionBlocks, cigarList):
            # make a paf line
            pafLine = f"{qContig}\t{contigLengths[qContig]}\t{qBlock[0]-1}\t{qBlock[1]}"
            pafLine += f"\t{orientation}"
            pafLine += f"\t{rContig}\t{contigLengths[rContig]}\t{rBlock[0]-1}\t{rBlock[1]}"
            pafLine += "\ttp:A:P" # it is assumed that only primary alignments are used
            pafLine += f"\tcg:Z:{makeCigarString(cigar)}"

            # make a new alignment based on the current projection
            alignment = Alignment(pafLine)
            subsetAlignments.append(alignment)
    subsetAlignments.sort(key = lambda x : (x.chromName, x.chromStart, x.chromEnd, x.contigName, x.contigStart, x.contigEnd))
    return  subsetAlignments
