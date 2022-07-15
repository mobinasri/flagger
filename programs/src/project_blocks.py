import sys
import argparse
from collections import defaultdict
import re

def getCigarList(cigarString):
    """
        Returns a list of tuples based on the cigar string. 
        Each tuple has two elements; the first one is showing 
        the operation which could be one of M, I or D and the 
        second element is an integer which shows the length of 
        the associated operation.
    """
    cigarOps = re.compile("[0-9]+").split(cigarString)[1:]
    # H (hard clipping) or S (Soft clipping) are not included since in the PAF format 
    # the start and end positions are shifted instead of adding H or S to the cigar string
    cigarSizes = [int(size) for size in re.compile("M|I|D|X|=").split(cigarString)[:-1]]
    cigarList = [ (op, size) for op, size in zip(cigarOps, cigarSizes)]
    return cigarList

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
        reverseBlocks.append(reverseInterval(blocks[i], contigLength))
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
            * forwardBlocks: A list of tuples which contains the blocks in the coordinates of the contig. Each block is a tuple of (start, end).
            * chromLength: (Integer) The total length of the reference chromosome
            * chromStart: (Integer) Where the alignment begins in the coordinates of the reference (1-based and inclusively)
            * contigLength: (Integer) The total length of the contig
            * contigStart: (Integer) Where the alignment begins in the coordinates of the query contig (1-based and inclusively)
            * contigEnd: (Integer) Where the alignment ends in the coordinates of the query contig (1-based and inclusively)
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
    # Return if the blocks has not overlap with the alignment
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
                projectionBlocks.append((projectionStartPos, projectionEndPos))
                if orientation == '+':
                    projectableBlocks.append((projectableStartPos, projectableEndPos))
                else:
                    projectableBlocks.append(reverseInterval((projectableStartPos, projectableEndPos), contigLength))
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
                projectionBlocks.append((projectionStartPos, projectionEndPos))
                if orientation == '+':
                    projectableBlocks.append((projectableStartPos, projectableEndPos))
                else:
                    projectableBlocks.append(reverseInterval((projectableStartPos, projectableEndPos), contigLength))
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
        projectionBlocks.append((projectionStartPos, projectionEndPos))
        if orientation == '+':
            projectableBlocks.append((projectableStartPos, projectableEndPos))
        else:
            projectableBlocks.append(reverseInterval((projectableStartPos, projectableEndPos), contigLength))
    return projectableBlocks, projectionBlocks
                


def main():
    parser = argparse.ArgumentParser(description='Given the alignments find the projection of a set of assembly blocks onto the reference (\'asm2ref\') or vice versa (\'ref2asm\')')
    parser.add_argument('--mode', type=str, default='asm2ref',
                    help='(Str) Default=\'asm2ref\' It can be either {\'ref2asm\' or \'asm2ref\'}. \
                           1.\'asm2ref\': \
                                    In this mode the blocks are given in the coordinates of the assembly and \
                                    the output will be the projections of those blocks onto the reference. \
                           2.\'ref2asm\': \
                                    In this mode the blocks are given in the coordinates of the reference and \
                                    the output will be the projections of those blocks onto the assembly.')
    parser.add_argument('--paf', type=str,
                    help='(PAF format) The alignments of the assembly to the reference. It should include the cigar format.')
    parser.add_argument('--blocks', type=str,
                    help='(BED format) The desired blocks in the assembly (or in the reference if the mode is \'ref2asm\'). It should be sorted and with no overlaps.')
    parser.add_argument('--outputProjectable', type=str,
                    help='(BED format) A path for saving the query blocks that could be projected. The projection of each line is available in the same line of the projection output')
    parser.add_argument('--outputProjection', type=str,
                    help='(BED format) A path for saving the projections of the query blocks.Note that the lines may not be sorted and may have overlaps because of its correspondence with the projected bed file. It is recommended to run bedtools sort (and merge) on this output')
    
    # Fetch the arguments
    args = parser.parse_args()
    mode = args.mode
    pafPath = args.paf
    blocksPath = args.blocks
    outputProjectable = args.outputProjectable
    outputProjection = args.outputProjection

    # Read the desired blocks. Their start and end positions are converted into the 1-based format
    blocks = defaultdict(list)
    with open(blocksPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            contigName = attrbs[0]
            # start is 0-based in bed format, it gets converted to 1-based here
            start = int(attrbs[1]) + 1
            end = int(attrbs[2])
            blocks[contigName].append((start, end))

    # Read the alignments one by one and for each of them find the projections by calling findProjections
    with open(pafPath,"r") as fPaf, open(outputProjection, "w") as fRef, open(outputProjectable, "w") as fQuery:
        for line in fPaf:
            # Extract the alignment attributes like the contig name, alignment boundaries, orientation and cigar 
            attrbs = line.strip().split()
            contigName = attrbs[0]
            contigLength = int(attrbs[1])
            contigStart = int(attrbs[2]) + 1 # 1-based
            contigEnd = int(attrbs[3])
            orientation = attrbs[4]
            chrom = attrbs[5]
            chromLength = int(attrbs[6])
            chromStart = int(attrbs[7]) + 1 # 1-based
            chromEnd = int(attrbs[8])
            # The cigar string starts after "cg:Z:"
            afterCg = line.strip().split("cg:Z:")[1]
            cigarString = afterCg.split()[0]
            cigarList = getCigarList(cigarString)
            #print(orientation,chromStart,chromEnd,cigarString)
            # rBlocks contains the projections and 
            # qBlocks contains the projectable blocks
            if mode == "asm2ref":
                if len(blocks[contigName]) == 0: # Continue if there is no block in the contig
                    continue
                #print(blocks[contigName], contigStart, contigEnd, chrom, chromStart, chromEnd)
                qBlocks, rBlocks = findProjections(mode, cigarList, blocks[contigName],
                                                    chromLength, chromStart, chromEnd, 
                                                    contigLength, contigStart, contigEnd,
                                                    orientation)
            else:
                if len(blocks[chrom]) == 0: # Continue if there is no block in the chrom
                    continue
                qBlocks, rBlocks = findProjections(mode, cigarList, blocks[chrom],
                                                    chromLength, chromStart, chromEnd,
                                                    contigLength, contigStart, contigEnd,
                                                    orientation)
            if mode == "asm2ref":
                for rBlock in rBlocks:
                    fRef.write("{}\t{}\t{}\n".format(chrom, rBlock[0] - 1, rBlock[1]))
                for qBlock in qBlocks:
                    fQuery.write("{}\t{}\t{}\n".format(contigName, qBlock[0] - 1, qBlock[1]))
            else: # mode = "ref2asm"
                for rBlock in rBlocks:
                    fRef.write("{}\t{}\t{}\n".format(contigName, rBlock[0] - 1, rBlock[1]))
                for qBlock in qBlocks:
                    fQuery.write("{}\t{}\t{}\n".format(chrom, qBlock[0] - 1, qBlock[1]))
main()
