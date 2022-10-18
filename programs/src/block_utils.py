import sys
import argparse
from collections import defaultdict
import re

CS_PATTERN = r'(:([0-9]+))|(([+-])([a-z]+)|([\\*]([a-z]+))+)'

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

def findProjections(mode, cigarList, forwardBlocks, 
                    chromLength, chromStart, chromEnd, 
                    contigLength, contigStart, contigEnd, 
                    orientation):
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
    # If the mode is 'asm2ref' we can simply call the same function but by changing the arguments in a clever manner!:
    #   1. The cigar operations should become with respect to the assembly contig not the reference
    #      to reach that aim we should call "convertIndelsInCigar"
    #   2. If the orientation is negative we should reverse the order of the cigar operations
    #   3. After this conversion we can assume that our reference is the assembly contig so the 
    #      the length, start and end positions of the reference and the assembly contig should be swapped
    #   4. Finally all we need is just to call the function again but in the mode of 'ref2asm'
    if mode == 'ref2asm':
        convertedCigarList = convertIndelsInCigar(cigarList)
        if orientation == '+':
            return findProjections('asm2ref', convertedCigarList, forwardBlocks,
                                   contigLength, contigStart, contigEnd,
                                   chromLength, chromStart, chromEnd, 
                                   orientation)
        else:
            return findProjections('asm2ref', convertedCigarList[::-1], forwardBlocks,
                                   contigLength, contigStart, contigEnd,
                                   chromLength, chromStart, chromEnd,
                                   orientation)
    blockIdx = 0
    projectableBlocks = []
    projectionBlocks = []
    nextOpStartRef = chromStart
    nextOpStartContig = None
    currOpStartRef = None
    currOpStartContig = None
    # Return if the blocks has no overlap with the alignment
    if (forwardBlocks[-1][1] < contigStart) or (contigEnd < forwardBlocks[0][0]):
        return projectionBlocks, projectableBlocks
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
    #print(cigarList)
    # iterate over cigar elements and find the projections
    for cigarOp, cigarSize in cigarList:
        #print(cigarOp,cigarSize)
        # Case 1: Mismatch or Match
        if cigarOp == 'M' or cigarOp == 'X' or cigarOp == '=':
            currOpStartContig = nextOpStartContig
            currOpStartRef = nextOpStartRef
            nextOpStartContig = currOpStartContig + cigarSize
            nextOpStartRef = currOpStartRef + cigarSize
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
                if currOpStartContig <= blocks[blockIdx][0]:
                    projectionStartPos = currOpStartRef + blocks[blockIdx][0] - currOpStartContig
                    projectableStartPos = blocks[blockIdx][0]
                projectionEndPos = currOpStartRef + blocks[blockIdx][1] - currOpStartContig
                projectableEndPos = blocks[blockIdx][1]
                info = blocks[blockIdx][2]
                projectionBlocks.append((projectionStartPos, projectionEndPos, info))
                if orientation == '+':
                    projectableBlocks.append((projectableStartPos, projectableEndPos, info))
                else:
                    reversedInterval = reverseInterval((projectableStartPos, projectableEndPos), contigLength)
                    projectableBlocks.append((reversedInterval[0], reversedInterval[1], info))
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            # In case of the 2nd scenario, the projection start position should be updated
            if (currOpStartContig <= blocks[blockIdx][0]) and (blocks[blockIdx][0] < nextOpStartContig):
                projectionStartPos = currOpStartRef + blocks[blockIdx][0] - currOpStartContig
                projectableStartPos = blocks[blockIdx][0]
        # Case 2: Insertion
        elif cigarOp == 'I':
            currOpStartContig = nextOpStartContig
            nextOpStartContig = currOpStartContig + cigarSize
            # For insertion there is no need to shift the reference
            currOpStartRef = nextOpStartRef
            # The beginning and endings of the blocks are excluded from the projection if they are insertions
            # So I trim the blocks to make their projections start and end in match/mismatch
            while (blockIdx < len(blocks)) and (blocks[blockIdx][1] < nextOpStartContig):
                # If a block is completely within an insertion then there is no valid projection for it
                if currOpStartContig <= blocks[blockIdx][0]:
                    blockIdx += 1
                    continue
                # If there is a block that ends in an insertion that inserted part will not be projected
                # So the end position of that projection will be one base before where the insertion starts
                # Note that one base before the insertion is absolutely an M
                projectionEndPos = currOpStartRef - 1
                projectableEndPos = currOpStartContig - 1
                info = blocks[blockIdx][2]
                projectionBlocks.append((projectionStartPos, projectionEndPos, info))
                if orientation == '+':
                    projectableBlocks.append((projectableStartPos, projectableEndPos, info))
                else:
                    reversedInterval = reverseInterval((projectableStartPos, projectableEndPos), contigLength)
                    projectableBlocks.append((reversedInterval[0], reversedInterval[1], info))
                blockIdx += 1
            if blockIdx >= len(blocks):
                break
            # If the last block starts with insertion but is not completely an insertion,
            # the initial inserted part is not projectable onto the reference so we start 
            # from the next operation (which should be an M)
            if (currOpStartContig <= blocks[blockIdx][0]) and (blocks[blockIdx][0] < nextOpStartContig):
                projectionStartPos = currOpStartRef
                projectableStartPos = nextOpStartContig
        # Case 3: Deletion
        elif cigarOp == 'D':
            nextOpStartRef += cigarSize
    # Note that nextOpStart(Contig | Ref) are not pointing to any cigar operation at this moment
    # since iterating over the cigar elements has finished. Those variables are now pointing to one base
    # after the end of the last cigar operation.
    # Handle the last block that has started but not finished by the end of the last cigar operation
    if (blockIdx < len(blocks)) and (blocks[blockIdx][0] < nextOpStartContig):
        projectionEndPos = nextOpStartRef - 1
        projectableEndPos = nextOpStartContig - 1
        info = blocks[blockIdx][2]
        projectionBlocks.append((projectionStartPos, projectionEndPos, info))
        if orientation == '+':
            projectableBlocks.append((projectableStartPos, projectableEndPos, info))
        else:
            reversedInterval = reverseInterval((projectableStartPos, projectableEndPos), contigLength)
            projectableBlocks.append((reversedInterval[0], reversedInterval[1], info))
    return projectableBlocks, projectionBlocks
                

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



def intersectInterval(interval_1, interval_2):
    """
        Receives two intervals; each interval is a tuple of (start, end) and
        returns their intersection
        If there is no overlap returns None
    """
    s1 = interval_1[0]
    e1 = interval_1[1]
    s2 = interval_2[0]
    e2 = interval_2[1]
    
    s = max(s1,s2)
    e = min(e1,e2)
    if e < s:
        return None
    return (s, e)
 

def subtractInterval(intervals, b):
    """
        Arguments:
            intervals: a sorted list of intervals; each interval is a tuple of (start, end)
            b: a single interval that will be subtracted from the list of intervals
            Note that start is 0-based closed and end is 0-based open
        Returns:
            a list of new intervals in which the interval, b, is absent
    """
    newIntervals = []
    for a in intervals:
        
        # Could be either of the two cases below
        # case 1: (with overlap) [ a - b ] is what should remain
        #
        # [      a      ]
        # [a - b][        b      ]
        #
        # or case 2: (with no overlap)
        #
        # [      a       ]
        # [    a - b     ]    [       b       ]
        if a[0] < b[0]:
            newInterval = (a[0], min(a[1], b[0]))
            newIntervals.append(newInterval)
            
        # Could be either of the two cases below
        # case 1: (with overlap)
        #
        #            [      a      ]
        # [        b      ][ a - b ]
        #
        # or case 2: (with no overlap)
        #
        #                     [      a       ]
        # [       b       ]   [     a - b    ]
        if b[1] < a[1]:
            newInterval = (max(a[0], b[1]), a[1])
            newIntervals.append(newInterval)
    return newIntervals


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
