from block_utils import *
from collections import defaultdict
from copy import deepcopy
import random
import numpy as np
import re
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class HomologyBlock:
    """
        A class for saving info of one block from the original assembly
        and also the name of the new contig it belongs to in the new 
        (falsified) assembly. This class also contains the relative order 
        of the block w.r.t to the other blocks in the newly generated contig.
    """

    def __init__(self, origCtg, origStart, origEnd, origStrand, newCtg, orderIndex):
        self.origCtg = origCtg # name of the original contig from the given assembly
        self.origStart = origStart # the 1-based start location of this block
        self.origEnd = origEnd # the 1-based start location of this block
        self.origStrand = origStrand # ['+' or '-']: if '-' the original block should be rev-complemented
        self.newCtg = newCtg # the name of the new contig where this block is localized in
        self.orderIndex = orderIndex # relative order of the block w.r.t to the other blocks in the new contig
        self.annotationBlockLists = defaultdict(BlockList)
        self.misAssemblyBlockList = BlockList()

        # attributes related to sampling locations
        # they will be updated with "setSamplingLengthAttributes" or "copySamplingLengthAttributes"
        self.annotationStartBlockListsForSampling = defaultdict(BlockList)
        self.annotationStartBlockLengthsForSampling = defaultdict(list)
        self.annotationStartTotalLengthsForSampling = defaultdict(int)

        # attributes for sampling misjoin locations
        self.annotationBlockListsForSamplingMisjoin = defaultdict(BlockList)
        self.annotationBlockLengthsForSamplingMisjoin = defaultdict(list)
        self.annotationTotalLengthsForSamplingMisjoin = defaultdict(int)

        self.misAssemblyLength = 0
        self.minOverlapRatioWithEachAnnotation = 0.5
        self.minMarginLength = 1000
        self.containsMisAssembly = False

    def getSequence(self, origCtgSequences, singleBaseErrorRate):
        origCtgSeq = origCtgSequences[self.origCtg]
        # note that start and end are 1-based closed
        blockSeq = origCtgSeq[self.origStart - 1: self.origEnd]

        # some parts of the block may be "Err"
        # split block into "Err" and non-"Err" parts since
        # "Err" parts are going to be changed and contain random
        # single-base errors
        blockSeqParts = []
        preE = 0
        for s, e, comp in self.misAssemblyBlockList.blocks:
            if comp == "Err":
                # add the previous sequence with no "Err" misassembly
                if preE < s - 1:
                    blockSeqParts.append(blockSeq[preE:s-1])
                # generate erroneous sequence and add this to the list
                erroneousSeq = induceSingleBaseErrors(blockSeq[s-1:e], singleBaseErrorRate)
                blockSeqParts.append(erroneousSeq)
                preE = e
        # add the last block with no "Err" misassembly
        if preE < len(blockSeq):
            blockSeqParts.append(blockSeq[preE:])

        finalBlockSeq = "".join(blockSeqParts)
        if self.origStrand == '-':
            finalBlockSeq = reverseComplement(finalBlockSeq)
        return finalBlockSeq

    def reverseAllBlockLists(self):
        """
            Make the coordinates of all the block lists in this homology block w.r.t to the end of whole block
        """
        blockLength = self.origEnd - self.origStart + 1
        for annotation in self.annotationBlockLists:
            self.annotationBlockLists[annotation].reverse(blockLength, inplace=True)
            # reverse blocks and lengths for sampling "Msj"
            self.annotationBlockListsForSamplingMisjoin[annotation].reverse(blockLength, inplace=True)
            self.annotationBlockLengthsForSamplingMisjoin[annotation].reverse()
            # reverse blocks and lengths for sampling "Sw", "Err", "Dup", "Col"
            self.annotationStartBlockListsForSampling[annotation].reverse(blockLength, inplace=True)
            self.annotationStartBlockLengthsForSampling[annotation].reverse()
        self.misAssemblyBlockList.reverse(blockLength, inplace=True)

    def clearAnnotationBlocksForSampling(self):
        for annotation in self.annotationBlockLists:
            self.annotationStartBlockListsForSampling[annotation] = BlockList([])
            self.annotationStartBlockLengthsForSampling[annotation] = []
            self.annotationStartTotalLengthsForSampling[annotation] = 0
            self.annotationBlockListsForSamplingMisjoin[annotation] = BlockList([])
            self.annotationBlockLengthsForSamplingMisjoin[annotation] = []
            self.annotationTotalLengthsForSamplingMisjoin[annotation] = 0

    def addMisAssemblyBlockList(self, blockList):
        """
        :param blockList: A BlockList containing the coordinates of the misassembled part
                          (Note that coordinates should be 1-based and relative to
                           the start position of this block rather than the original contig.
                           The "shift" method of "BlockList" can be used for shifting
                           coordinates prior to calling this method)
        """
        self.misAssemblyBlockList.extend(blockList)

    def addAnnotationBlockList(self, name, blockList):
        """
        :param name: The name of the annotation to add
        :param blockList: A BlockList containing the coordinates of the annotation
                          (Note that coordinates should be 1-based and relative to
                           the start position of this block rather than the original contig.
                           The "shift" method of "BlockList" can be used for shifting
                           coordinates prior to adding blockList)
        """
        self.annotationBlockLists[name] = blockList.copy()

    def setSamplingLengthAttributes(self, misAssemblyLength, minOverlapRatioWithEachAnnotation, minMarginLength):
        """
        :param misAssemblyLength: The length of misassemblies that are going to be generated after calling this method
        :param minOverlapRatioWithEachAnnotation: Minimum overlap ratio each misassembly should have with the desired annotation
        :param minMarginLength: The minimum margin length from both ends of the whole block
        """
        self.misAssemblyLength = misAssemblyLength
        self.minOverlapRatioWithEachAnnotation = minOverlapRatioWithEachAnnotation
        self.minMarginLength = minMarginLength

    def copySamplingLengthAttributes(self, otherBlock):
        """
        copy length attributes related to sampling process from the other block
        """
        self.setSamplingLengthAttributes(otherBlock.misAssemblyLength,
                                         otherBlock.minOverlapRatioWithEachAnnotation,
                                         otherBlock.minMarginLength)

    def updateAnnotationBlocksForSampling(self):
        """
        This method performs two main truncations and one extension on the annotation blocks to find where
        the start locations of misassemblies are going to be sampled from:
            1. Truncate each continuous annotation block from the right side to exclude start locations for which
               the misassembly overlap will be less than "minOverlapRatioWithEachAnnotation"
            2. Extend each continuous annotation block from the left side as far as the overlap limitation allows
            3. Truncate each continuous annotation block from both sides so that it does not have
               any overlap with the margins of the whole block.
               The misassembly length should be added to the right margin to make sure that
               the whole misassembled segment is fully within the block

        The attributes that end with "Misjoin" are for sampling the locations where misjoins can happen
        The truncation (1.) and extension (2.) are not applied on annotation blocks to find misjoin
        locations since there is no overlap restriction for misjoins (each misjoin is a single-point misassembly)
        """

        assert(0.5 <= self.minOverlapRatioWithEachAnnotation)
        #The length to be excluded from sampling in each continuous annotation
        #block (from the right side)
        lengthToTruncateFromRight = int(self.misAssemblyLength * self.minOverlapRatioWithEachAnnotation)

        # The length to be extended from the left side
        lengthToExtendFromLeft = self.misAssemblyLength - lengthToTruncateFromRight - 1
        # The margin of the whole block to be excluded from sampling
        # (from left side)
        wholeBlockMarginFromLeft = self.minMarginLength
        # The margin of the whole block to be excluded from sampling
        # (from right side)
        wholeBlockMarginFromRight = self.misAssemblyLength + self.minMarginLength - 1

        wholeBlockWithoutMargins = BlockList([(1, self.origEnd - self.origStart + 1)]).truncateFromLeft(wholeBlockMarginFromLeft, inplace=False)
        wholeBlockWithoutMargins.truncateFromRight(wholeBlockMarginFromRight, inplace=True)

        # margins for misjoin are more inclusive on the right side
        wholeBlockWithoutMarginsForMisjoin = BlockList([(1, self.origEnd - self.origStart + 1)]).truncateFromLeft(wholeBlockMarginFromLeft, inplace=False)
        wholeBlockWithoutMarginsForMisjoin.truncateFromRight(self.minMarginLength, inplace=True)

        for name in self.annotationBlockLists:
            ###############################################################
            ### Construct sampling blocks for the start locations of  #####
            ###       misassemblies other than misjoins               #####
            ###############################################################

            self.annotationStartBlockListsForSampling[name] = self.annotationBlockLists[name].copy()

            # remove blocks shorter than minimum overlap size
            self.annotationStartBlockListsForSampling[name].removeBlocksShorterThan(lengthToTruncateFromRight, inplace=True)

            # truncate the right side of the blocks
            self.annotationStartBlockListsForSampling[name].truncateFromRight(lengthToTruncateFromRight, inplace=True)

            # extend blocks from the left side
            self.annotationStartBlockListsForSampling[name].extendFromLeft(lengthToExtendFromLeft,
                                                                           leftMostPosition= 1,
                                                                           inplace=True)


            # truncate the blocks to make sure they do not have overlap with the margins
            self.annotationStartBlockListsForSampling[name].intersect(wholeBlockWithoutMargins, inplace=True)
            self.annotationStartBlockLengthsForSampling[name] = [i[1] - i[0] + 1 for i in self.annotationStartBlockListsForSampling[name].blocks]
            self.annotationStartTotalLengthsForSampling[name] = sum(self.annotationStartBlockLengthsForSampling[name])

            ###############################################################
            ####    Construct sampling blocks for the locations of    #####
            ####         misassemblies other than misjoins            #####
            ###############################################################

            # Since misjoin is a single point event
            # there is no need to extend/truncate the annotation blocks
            # from left or right to meet the overlap limitation
            # Additionally there is no need
            # to filter short annotation blocks
            self.annotationBlockListsForSamplingMisjoin[name] = self.annotationBlockLists[name].copy()

            # truncate the blocks to make sure they do not have overlap with the margins
            self.annotationBlockListsForSamplingMisjoin[name].intersect(wholeBlockWithoutMarginsForMisjoin, inplace=True)
            self.annotationBlockLengthsForSamplingMisjoin[name] = [i[1] - i[0] + 1 for i in self.annotationBlockListsForSamplingMisjoin[name].blocks]
            self.annotationTotalLengthsForSamplingMisjoin[name] = sum(self.annotationBlockLengthsForSamplingMisjoin[name])


    def sampleMisjoinLocation(self, name):
        selectedInterval = random.choices(self.annotationBlockListsForSamplingMisjoin[name].blocks,
                                               weights=self.annotationBlockLengthsForSamplingMisjoin[name],
                                               k=1)[0]
        selectedLocation = random.randint(selectedInterval[0], selectedInterval[1])

        return selectedLocation

    def sampleMisAssemblyInterval(self, name, misAssemblyLength):
        # For getting a random interval for any misassembly
        # other than misjoin (which is a single-point event),
        # this method can be called
        # self.misAssemblyLength was set the last time
        # "updateAnnotationBlocksForSampling" was invoked
        # If the previous value does match the current one
        # return None since it is essential to update the
        # annotation start locations for sampling based on
        # the correct misassembly length (This is mainly for
        # the overlap criteria)
        if self.misAssemblyLength != misAssemblyLength:
            return None
        selectedStartInterval = random.choices(self.annotationStartBlockListsForSampling[name].blocks,
                                               weights=self.annotationStartBlockLengthsForSampling[name],
                                               k=1)[0]
        selectedStartPos = random.randint(selectedStartInterval[0], selectedStartInterval[1])
        selectedEndPos = selectedStartPos + misAssemblyLength - 1

        return selectedStartPos, selectedEndPos

    def extractAnnotationsFromParentBlock(self, parentBlock, start, end):
        """
        This function is useful for moving annotation coordinates from the parent block to the
        child one (self), which is one part of the parent block

        :param parentBlock: The block which was split and then one part of it is the current block
        :param start: 1-based location of the parent block which is the first base in this block
        :param end: 1-based location of the parent block which is the last base in this block
        """
        for name, blockList in parentBlock.annotationBlockLists.items():
            subsetBlockList = blockList.intersect(BlockList([(start, end)]), inplace=False)
            subsetBlockList.shift( -(start - 1), minCoordinate = 1, maxCoordinate = end - start + 1, inplace = True)
            self.addAnnotationBlockList(name, subsetBlockList)

    def extractMisAssembliesFromParentBlock(self, parentBlock, start, end):
        """
        This function is useful for moving annotation coordinates from the parent block to the
        child one (self), which is one part of the parent block

        :param parentBlock: The block which was split and then one part of it is the current block
        :param start: 1-based location in the parent block which is the first base in this block
        :param end: 1-based location in the parent block which is the last base in this block
        """
        subsetBlockList = parentBlock.misAssemblyBlockList.intersect(BlockList([(start, end)]), inplace=False)
        subsetBlockList.shift( -(start - 1), minCoordinate = 1, maxCoordinate = end - start + 1, inplace = True)
        self.addMisAssemblyBlockList(subsetBlockList)



class HomologyRelation:
    """
        A class for saving two homologous blocks and the alignment between them
    """

    def __init__(self, block: HomologyBlock, homologousBlock: HomologyBlock, cigarList: list, orientation: str):
        self.block = block
        self.homologousBlock = homologousBlock
        if cigarList != None and len(cigarList) > 0:
            blockLen = block.origEnd - block.origStart + 1
            homologousBlockLen = homologousBlock.origEnd - homologousBlock.origStart + 1
            alignment = Alignment(f"query_ctg\t{homologousBlockLen}\t0\t{homologousBlockLen}\t{orientation}\tref_ctg\t{blockLen}\t0\t{blockLen}\t0\t0\t60\tcg:Z:{makeCigarString(cigarList)}\ttp:A:P")
            self.alignment = alignment
        else:
            self.alignment = None

    @staticmethod
    def createRef2QueryRelationFromAlignment(alignment: Alignment, newCtgSuffix: str):
        """

        :param alignment: An alignment to be used for creating the homology blocks and the homology relation
        :param newCtgSuffix: A suffix to be added to the origCtg to get the new contig name
        :return: A homology relation created based on the given alignment. This homology block
                 contains the reference/hap1 block as the "block" attribute and the query/hap2 block
                 as "homologousBlock" attribute, so it is representing the relation from reference to query
        """
        block = HomologyBlock(alignment.chromName,
                              alignment.chromStart + 1,
                              alignment.chromEnd,
                              '+',
                              f'{alignment.chromName}{newCtgSuffix}',
                              None)
        homologousBlock = HomologyBlock(alignment.contigName,
                                        alignment.contigStart + 1,
                                        alignment.contigEnd,
                                        '+',
                                        f'{alignment.contigName}{newCtgSuffix}',
                                        None)
        homologyRelation = HomologyRelation(block, homologousBlock, alignment.cigarList, alignment.orientation)
        return  homologyRelation

    @staticmethod
    def createVoidRelationFromInterval(origCtg: str, origCtgStart: int, origCtgEnd: int, newCtgSuffix: str):
        """
        :param origCtg: The name of the original contig this block is taken from
        :param origCtgStart: The 1-based (closed) start position of this block in the original contig
        :param origCtgEnd: The 1-base (closed) end position of this block in the original contig
        :param newCtgSuffix: The suffix to be added to the origCtg to get the new contig name
        :return: a void homology relation, which contains no homologous block and no alignment
        """
        block = HomologyBlock(origCtg,
                              origCtgStart,
                              origCtgEnd,
                              '+',
                              f'{origCtg}{newCtgSuffix}',
                              None)
        voidRelation = HomologyRelation(block, None, None, None)
        return  voidRelation

    def splitIntoThreeParts(self, start, end):
        """
        Split the relation into three smaller relations based on the given coordinates
        IMPORTANT NOTE:
        This method only works correctly when the strands of both "block" and "homologousBlock" are positive


        :param start: The start location of the middle part (1-based and in the block coordinates).
                      "start" cannot be 1 (Part 1 cannot be empty)
        :param end: The end location of the middle part (1-based and in the block coordinates).
                    "end" cannot be the last location of the block (Part 3 cannot be empty)
        :return: A list of three homology relations
        """
        assert(self.block.origStrand == '+')
        assert(self.homologousBlock.origStrand == '+')
        assert(1 < start)
        assert(end < self.alignment.chromLength)

        forwardBlocks = [(1, start - 1, ""), (start, end, ""), (end + 1, self.alignment.chromLength, "")]
        includeEndingIndel = True
        includePostIndel = True
        projectableBlocks, projectionBlocks, cigarLists = \
            findProjections('ref2asm', self.alignment.cigarList, forwardBlocks,
                            self.alignment.chromLength, self.alignment.chromStart + 1, self.alignment.chromEnd,
                            self.alignment.contigLength, self.alignment.contigStart + 1, self.alignment.contigEnd,
                            self.alignment.orientation, includeEndingIndel, includePostIndel)

        assert(len(projectableBlocks) == 3)
        assert(len(projectionBlocks) == 3)
    
        
        projectionsOrigCoor = []
        projectionsRelCoor = []
        for i in range(3):
            if projectionBlocks[i][0] == None or projectionBlocks[i][1] == None:
                return None
            projectionsOrigCoor.append([projectableBlocks[i][0] + self.block.origStart - 1,
                                        projectableBlocks[i][1] + self.block.origStart - 1,
                                        projectionBlocks[i][0] + self.homologousBlock.origStart - 1,
                                        projectionBlocks[i][1] + self.homologousBlock.origStart - 1,
                                        cigarLists[i]])
            projectionsRelCoor.append([projectableBlocks[i][0],
                                       projectableBlocks[i][1],
                                       projectionBlocks[i][0],
                                       projectionBlocks[i][1],
                                       cigarLists[i]])
        # sort projections by start position of the ref haplotype
        projectionsOrigCoor.sort(key = lambda x : x[0])
        projectionsRelCoor.sort(key = lambda x : x[0])

        rBlock = self.block
        qBlock = self.homologousBlock

        # create one homology block per projection
        # the middle projection will be used for switching

        # ref blocks

        rBlockPart1 =  HomologyBlock(rBlock.origCtg,
                                     projectionsOrigCoor[0][0],
                                     projectionsOrigCoor[0][1],
                                     '+',
                                     rBlock.newCtg,
                                     rBlock.orderIndex)
        rBlockPart1.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[0][0], projectionsRelCoor[0][1])
        rBlockPart1.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[0][0], projectionsRelCoor[0][1])
        rBlockPart1.copySamplingLengthAttributes(rBlock)
        rBlockPart1.updateAnnotationBlocksForSampling()

        rBlockPart2 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[1][0],
                                    projectionsOrigCoor[1][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 1)
        rBlockPart2.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.copySamplingLengthAttributes(rBlock)
        rBlockPart2.updateAnnotationBlocksForSampling()

        rBlockPart3 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[2][0],
                                    projectionsOrigCoor[2][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 2)
        rBlockPart3.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])
        rBlockPart3.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])
        rBlockPart3.copySamplingLengthAttributes(rBlock)
        rBlockPart3.updateAnnotationBlocksForSampling()

        # query blocks
        qOrderIndexPart1 = qBlock.orderIndex if self.alignment.orientation  == '+' else qBlock.orderIndex + 2
        qOrderIndexPart2 = qBlock.orderIndex + 1 if self.alignment.orientation  == '+' else qBlock.orderIndex + 1
        qOrderIndexPart3 = qBlock.orderIndex + 2 if self.alignment.orientation  == '+' else qBlock.orderIndex

        qBlockPart1 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[0][2],
                                    projectionsOrigCoor[0][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart1)
        qBlockPart1.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])
        qBlockPart1.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])
        qBlockPart1.clearAnnotationBlocksForSampling()
        #qBlockPart1.copySamplingLengthAttributes(qBlock)
        #qBlockPart1.updateAnnotationBlocksForSampling()

        qBlockPart2 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[1][2],
                                    projectionsOrigCoor[1][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart2)
        qBlockPart2.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.clearAnnotationBlocksForSampling()
        #qBlockPart2.copySamplingLengthAttributes(qBlock)
        #qBlockPart2.updateAnnotationBlocksForSampling()

        qBlockPart3 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[2][2],
                                    projectionsOrigCoor[2][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart3)
        qBlockPart3.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[2][2], projectionsRelCoor[2][3])
        qBlockPart3.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[2][2], projectionsRelCoor[2][3])
        qBlockPart3.clearAnnotationBlocksForSampling()
        #qBlockPart3.copySamplingLengthAttributes(qBlock)
        #qBlockPart3.updateAnnotationBlocksForSampling()

        relationPart1 = HomologyRelation(rBlockPart1,
                                         qBlockPart1,
                                         projectionsOrigCoor[0][4],
                                         self.alignment.orientation)
        relationPart2 = HomologyRelation(rBlockPart2,
                                         qBlockPart2,
                                         projectionsOrigCoor[1][4],
                                         self.alignment.orientation)
        relationPart3 = HomologyRelation(rBlockPart3,
                                         qBlockPart3,
                                         projectionsOrigCoor[2][4],
                                         self.alignment.orientation)

        homologyRelations = [relationPart1, relationPart2, relationPart3]

        return  homologyRelations

    def splitIntoTwoParts(self, loc):
        """
        Split the relation into two smaller relations based on the given location
        IMPORTANT NOTE:
        This method only works correctly when the strands of both "block" and "homologousBlock" are positive

        :param loc: The location where the relation should be split; [,loc], [loc+1,]
        :return: A list of two homology relations
        """
        assert(self.block.origStrand == '+')
        assert(self.homologousBlock.origStrand == '+')
        assert(1 < loc)
        assert(loc < self.alignment.chromLength)

        forwardBlocks = [(1, loc, ""), (loc+1, self.alignment.chromLength, "")]
        includeEndingIndel = True
        includePostIndel = True
        projectableBlocks, projectionBlocks, cigarLists = \
            findProjections('ref2asm', self.alignment.cigarList, forwardBlocks,
                            self.alignment.chromLength, self.alignment.chromStart + 1, self.alignment.chromEnd,
                            self.alignment.contigLength, self.alignment.contigStart + 1, self.alignment.contigEnd,
                            self.alignment.orientation, includeEndingIndel, includePostIndel)

        assert(len(projectableBlocks) == 2)
        assert(len(projectionBlocks) == 2)


        projectionsOrigCoor = []
        projectionsRelCoor = []
        for i in range(2):
            if projectionBlocks[i][0] == None or projectionBlocks[i][1] == None:
                return None
            projectionsOrigCoor.append([projectableBlocks[i][0] + self.block.origStart - 1,
                                        projectableBlocks[i][1] + self.block.origStart - 1,
                                        projectionBlocks[i][0] + self.homologousBlock.origStart - 1,
                                        projectionBlocks[i][1] + self.homologousBlock.origStart - 1,
                                        cigarLists[i]])
            projectionsRelCoor.append([projectableBlocks[i][0],
                                       projectableBlocks[i][1],
                                       projectionBlocks[i][0],
                                       projectionBlocks[i][1],
                                       cigarLists[i]])
        # sort projections by start position of the ref haplotype
        projectionsOrigCoor.sort(key = lambda x : x[0])
        projectionsRelCoor.sort(key = lambda x : x[0])

        rBlock = self.block
        qBlock = self.homologousBlock

        # create one homology block per projection
        # the middle projection will be used for switching

        # ref blocks
        rBlockPart1 =  HomologyBlock(rBlock.origCtg,
                                     projectionsOrigCoor[0][0],
                                     projectionsOrigCoor[0][1],
                                     '+',
                                     rBlock.newCtg,
                                     rBlock.orderIndex)
        rBlockPart1.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[0][0], projectionsRelCoor[0][1])
        rBlockPart1.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[0][0], projectionsRelCoor[0][1])
        rBlockPart1.copySamplingLengthAttributes(rBlock)
        rBlockPart1.updateAnnotationBlocksForSampling()

        rBlockPart2 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[1][0],
                                    projectionsOrigCoor[1][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 1)
        rBlockPart2.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.copySamplingLengthAttributes(rBlock)
        rBlockPart2.updateAnnotationBlocksForSampling()

        # query blocks
        qOrderIndexPart1 = qBlock.orderIndex if self.alignment.orientation  == '+' else qBlock.orderIndex + 1
        qOrderIndexPart2 = qBlock.orderIndex + 1 if self.alignment.orientation  == '+' else qBlock.orderIndex

        qBlockPart1 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[0][2],
                                    projectionsOrigCoor[0][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart1)
        qBlockPart1.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])
        qBlockPart1.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])
        qBlockPart1.clearAnnotationBlocksForSampling()
        #qBlockPart1.copySamplingLengthAttributes(qBlock)
        #qBlockPart1.updateAnnotationBlocksForSampling()

        qBlockPart2 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[1][2],
                                    projectionsOrigCoor[1][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart2)
        qBlockPart2.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.clearAnnotationBlocksForSampling()
        #qBlockPart2.copySamplingLengthAttributes(qBlock)
        #qBlockPart2.updateAnnotationBlocksForSampling()


        relationPart1 = HomologyRelation(rBlockPart1,
                                         qBlockPart1,
                                         projectionsOrigCoor[0][4], #4th element is cigar
                                         self.alignment.orientation)
        relationPart2 = HomologyRelation(rBlockPart2,
                                         qBlockPart2,
                                         projectionsOrigCoor[1][4], #4th element is cigar
                                         self.alignment.orientation)

        homologyRelations = [relationPart1, relationPart2]

        return  homologyRelations

class HomologyRelationChains:
    def __init__(self, alignments, origContigLengths, origCtgListRefOnly, newCtgSuffix):
        """
            A class for saving chains of homology relations for each new contig

            origCtgListRefOnly: The names of the original contigs from ref/hap1 only
            For other params read the documentation for "createAllInclusiveRelationChainsFromAlignments"
        """
        self.relationChains = HomologyRelationChains.createAllInclusiveRelationChainsFromAlignments(alignments,
                                                                                                    origContigLengths,
                                                                                                    newCtgSuffix)
        # a dictionary with annotation names as keys and list of
        # weights as values; one weight per newCtg
        # it will be filled by calling "updateAnnotationBlocksForSampling"
        self.origRefContigNames = origCtgListRefOnly.copy()
        self.newCtgAnnotationWeightsForSampling = defaultdict(list)
        self.newCtgAnnotationWeightsForSamplingMisjoin = defaultdict(list)
        self.newCtgListForSampling = [c + newCtgSuffix for c in origCtgListRefOnly]
        self.newCtgToIndexForSampling = {c: i for i, c in enumerate(self.newCtgListForSampling)}
        self.misAssemblyLength = 0

    @staticmethod
    def createAllInclusiveRelationChainsFromAlignments(alignments: list, contigLengths: dict, newCtgSuffix: str) -> defaultdict:
        """

        :param alignments: A list of alignments to be used for creating homology blocks. It is recommended that alignments
                            to be free of any overlap either on the reference/hap1 or query/hap2 coordinates.
                            These functions can be used for filtering overlapping alignments:
                                - "getBlockListsWithSingleAlignmentPerRefContig"
                                - "getBlockListsWithSingleAlignmentPerQueryContig"
                                - "subsetAlignmentsToRefBlocks"
                                - "subsetAlignmentsToQueryBlocks"
        :param contigLengths: A dictionary of contig lengths with original contig names as keys and
                              original contig lengths as values
        :param newCtgSuffix: A suffix to be added to the origCtg to get the new contig name
        :return: A dictionary of relations chains with new contig names as keys and relation lists as values
        """

        # a dictionary with new contig names as keys and
        # the relation/alignment free intervals saved as
        # in BlockLists as values
        relationFreeIntervals = {}
        for origCtg, origCtgLen in contigLengths.items():
            #BlockList receives 1-based coordinates
            print(origCtg + newCtgSuffix)
            relationFreeIntervals[origCtg + newCtgSuffix] = BlockList([(1, origCtgLen)])

        relationChains = defaultdict(list)

        # create two-way relations; ref2query and query2ref and add them to the dictionary
        for alignment in alignments:
            # the homology blocks (non-void ones) are shared between ref2queryRelation and query2refRelation
            # to be able to have access to the updated index of the homologous block
            ref2queryRelation = HomologyRelation.createRef2QueryRelationFromAlignment(alignment, newCtgSuffix)
            query2refRelation = HomologyRelation(ref2queryRelation.homologousBlock,
                                                 ref2queryRelation.block,
                                                 None,
                                                 None)
            relationChains[ref2queryRelation.block.newCtg].append(ref2queryRelation)
            relationChains[query2refRelation.block.newCtg].append(query2refRelation)

            # update relation-free intervals
            # to create void relations from them
            # after the alignments are all iterated
            relationFreeIntervals[ref2queryRelation.block.newCtg].subtract(BlockList([(alignment.chromStart + 1, alignment.chromEnd)]), inplace=True)
            relationFreeIntervals[query2refRelation.block.newCtg].subtract(BlockList([(alignment.contigStart + 1, alignment.contigEnd)]), inplace=True)

        # create void relations for the intervals without alignments
        for origCtg in contigLengths:
            for interval in relationFreeIntervals[origCtg + newCtgSuffix].blocks:
                voidRelation = HomologyRelation.createVoidRelationFromInterval(origCtg, interval[0], interval[1], newCtgSuffix)
                relationChains[voidRelation.block.newCtg].append(voidRelation)

        # sort the relations for each contig based on start coordinates
        for newCtg in relationChains:
            relationChains[newCtg].sort(key = lambda x : x.block.origStart)
            for i, relation in enumerate(relationChains[newCtg]):
                relation.block.orderIndex = i
        return relationChains

    def fillAnnotationBlockListsFromOriginalContigs(self, annotationBlockListsPerOrigContig, origContigLengths, newCtgSuffix):
        """
        This method takes annotations in the original contig coordinates and extract the annotations related to the "block"
        attribute of each homology relation

        :param annotationBlockListsPerOrigContig: A nested dictionary that contains annotations in the original
                                                  contig coordinates. It has origCtg names as keys and a dictionary
                                                  of annotations per value like below:
                                                  {"contig1": {"annot1" : BlockList([(1,10), (20,30)]),
                                                               "annot2": BlockList([(11, 19)])}}
        :param origContigLengths: A dictionary of contig lengths with original contig names as keys and
                                  original contig lengths as values
        :param newCtgSuffix: A suffix to be added to the origCtg to get the new contig name
        """
        for origCtg, annotationBlockLists in annotationBlockListsPerOrigContig.items():
            newCtgName = origCtg + newCtgSuffix
            # create a homology block for the whole original contig,
            # including the annotations
            if origCtg not in origContigLengths:
                print(f"{origCtg} is present in annotation but not in the fasta file so its annotation is skipped!")
                continue
            wholeOrigContigBlock = HomologyBlock(origCtg, 1, origContigLengths[origCtg], '+', newCtgName, 0)
            for name, blockList in annotationBlockLists.items():
                wholeOrigContigBlock.addAnnotationBlockList(name, blockList)
            # the created homology block will then be used for extracting the
            # annotations related to each relation.block
            for relation in self.relationChains[newCtgName]:
                relation.block.extractAnnotationsFromParentBlock(wholeOrigContigBlock,
                                                                 relation.block.origStart,
                                                                 relation.block.origEnd)

    def induceSwitchMisAssembly(self, newCtg, orderIndex, switchStart, switchEnd, switchEffectWindowLength):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # to make sure margin length is greater than or equal to
        # the size of the blocks created beside the edges of the switch errors
        assert(switchEffectWindowLength <= relationToSplit.block.minMarginLength)

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(switchStart, switchEnd)

        if ref2querySplitRelations == None:
            return False

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)


        # the blocks that have to be swapped
        relationPart2 = ref2querySplitRelations[1]
        rBlockPart2 = relationPart2.block
        qBlockPart2 = relationPart2.homologousBlock

        # swap blocks, order indices and new contig names
        rBlockPart2.orderIndex, qBlockPart2.orderIndex = qBlockPart2.orderIndex,  rBlockPart2.orderIndex
        rBlockPart2.newCtg, qBlockPart2.newCtg = qBlockPart2.newCtg,  rBlockPart2.newCtg
        relationPart2.block, relationPart2.homologousBlock = relationPart2.homologousBlock, relationPart2.block

        # convert cigar operations (DEL to INS and INS to DEL) since the blocks are switched
        relationPart2.alignment.cigarList = convertIndelsInCigar(relationPart2.alignment.cigarList)

        # update the strand orientation of the switched blocks
        # and reverse the coordinates
        # if the alignment's orientation was negative
        if relationPart2.alignment.orientation == '-':
            rBlockPart2.reverseAllBlockLists()
            qBlockPart2.reverseAllBlockLists()
            rBlockPart2.origStrand = '-'
            qBlockPart2.origStrand = '-'


        ### Adding misassembly blocks expected to be represented in read alignments ####

        # add misassembly with "Sw" label to both ends of the middle reference block
        rBlockPart2Length = rBlockPart2.origEnd - rBlockPart2.origStart + 1
        if rBlockPart2Length <= 2.5 * switchEffectWindowLength:
            rBlockPart2.addMisAssemblyBlockList(BlockList([(1, rBlockPart2Length, "Sw")]))
        else:
            rBlockPart2.addMisAssemblyBlockList(BlockList([(1, switchEffectWindowLength, "Sw"),
                                                           (rBlockPart2Length - switchEffectWindowLength + 1,
                                                            rBlockPart2Length,
                                                            "Sw")]))
        rBlockPart2.containsMisAssembly= True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationBlocksForSampling()

        # add misassembly with "Sw" label to both ends of the middle query block
        qBlockPart2Length = qBlockPart2.origEnd - qBlockPart2.origStart + 1
        if qBlockPart2Length <= 2.5 * switchEffectWindowLength:
            qBlockPart2.addMisAssemblyBlockList(BlockList([(1, qBlockPart2Length, "Sw")]))
        else:
            qBlockPart2.addMisAssemblyBlockList(BlockList([(1, switchEffectWindowLength, "Sw"),
                                                           (qBlockPart2Length - switchEffectWindowLength + 1,
                                                            qBlockPart2Length,
                                                            "Sw")]))
        qBlockPart2.containsMisAssembly= True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        qBlockPart2.clearAnnotationBlocksForSampling()


        # add misassembly with the "Sw" label to the end part of the left reference block
        rBlockPart1 = ref2querySplitRelations[0].block
        rBlockPart1Length = rBlockPart1.origEnd - rBlockPart1.origStart + 1
        if rBlockPart1Length <= switchEffectWindowLength:
            rBlockPart1.addMisAssemblyBlockList(BlockList([(1, rBlockPart1Length, "Sw")]))
        else:
            rBlockPart1.addMisAssemblyBlockList(BlockList([(rBlockPart1Length - switchEffectWindowLength + 1,
                                                            rBlockPart1Length,
                                                            "Sw")]))


        # add misassembly with the "Sw" label to the beginning part of the right reference block
        rBlockPart3 = ref2querySplitRelations[2].block
        rBlockPart3Length = rBlockPart3.origEnd - rBlockPart3.origStart + 1
        if rBlockPart3Length <= switchEffectWindowLength:
            rBlockPart3.addMisAssemblyBlockList(BlockList([(1, rBlockPart3Length, "Sw")]))
        else:
            rBlockPart3.addMisAssemblyBlockList(BlockList([(1, switchEffectWindowLength, "Sw")]))


        # add misassembly with the "Sw" label to the left query block
        # it will be added to the end if orientation was positive
        # it will be added to the beginning if orientation was negative
        qBlockPart1 = ref2querySplitRelations[0].homologousBlock
        qBlockPart1Length = qBlockPart1.origEnd - qBlockPart1.origStart + 1
        if qBlockPart1Length <= switchEffectWindowLength:
            qBlockPart1.addMisAssemblyBlockList(BlockList([(1, qBlockPart1Length, "Sw")]))
        else:
            if relationToSplit.alignment.orientation == '+':
                qBlockPart1.addMisAssemblyBlockList(BlockList([(qBlockPart1Length - switchEffectWindowLength + 1,
                                                                qBlockPart1Length,
                                                                "Sw")]))
            else: # '-'
                qBlockPart1.addMisAssemblyBlockList(BlockList([(1, switchEffectWindowLength, "Sw")]))


        # add misassembly with the "Sw" label to the right query block
        # it will be added to the beginning if orientation was positive
        # it will be added to the end if orientation was negative
        qBlockPart3 = ref2querySplitRelations[2].homologousBlock
        qBlockPart3Length = qBlockPart3.origEnd - qBlockPart3.origStart + 1
        if qBlockPart3Length <= switchEffectWindowLength:
            qBlockPart3.addMisAssemblyBlockList(BlockList([(1, qBlockPart3Length, "Sw")]))
        else:
            if relationToSplit.alignment.orientation == '+':
                qBlockPart3.addMisAssemblyBlockList(BlockList([(1, switchEffectWindowLength, "Sw")]))
            else: # '-'
                qBlockPart3.addMisAssemblyBlockList(BlockList([(qBlockPart3Length - switchEffectWindowLength + 1,
                                                                qBlockPart3Length,
                                                                "Sw")]))



        # create the equivalent list of relations from query to ref
        # these relations will show the same connections between blocks
        # but in the other way around
        query2refSplitRelations = []
        if relationToSplit.alignment.orientation == '+':
            for relation in ref2querySplitRelations:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)
        else:
            for relation in ref2querySplitRelations[::-1]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)

        # insert split relations to relation chain of the "newCtg"
        for relation in ref2querySplitRelations:
            self.relationChains[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        # insert split relations to relation chain of the "otherHapNewCtg"
        for relation in query2refSplitRelations:
            self.relationChains[otherHapNewCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[otherHapNewCtg][otherHapOrderIndex + 3:]:
            relation.block.orderIndex += 2

        
        self.updateNewCtgAnnotationWeightsForSampling(newCtg, relationToSplit, ref2querySplitRelations)

        return True

    def breakIntoTwoContigs(self, newCtg, lastOrderIndexOnLeft):
        """
        split the relation chain of the new contig into two parts
        remove new contig from the chains table and insert two smaller contigs/chains
        (Reminder: each contig is represented as a sorted chain of relations)

        The split contigs are named according to the below text:

        if the contig has not been split before
        add "_p1"/"_p2" to the end of the name; one suffix for each part
        if it was split before, so it already has "_p1"/"_p2" then
        just add an extension ".1"/".2" (1 for the left part and 2 for the right part)
        therefore the final name of split contigs can be like examples below:

        ${ORIGINAL_NAME}_f_p2 or ${ORIGINAL_NAME}_f_p1 which means that it was split only once

        or like ${ORIGINAL_NAME}_f_p1.1.2 which means that it was split three times "_p1" -> ".1" -> ".2"
        for the first time it was on the left side ("_p1"),
        for the second time on the left side (".1")
        and for the third time it was on the right side (".2")

        :param newCtg: The name of the new contig (parent chain) to break
        :param lastOrderIndexOnLeft: The relations of the new contig will be split into two parts:
                                        - part 1: From index 0 to lastOrderIndexOnLeft inclusively
                                        - part 2: From index lastOrderIndexOnLeft + 1 till the end inclusively
        """
        z = re.findall("(?<=\.p)[0-9.]+$", newCtg)
        if len(z) == 1:
            newCtgLeft = newCtg + "_1"
            newCtgRight = newCtg + "_2"
        else:
            newCtgLeft = newCtg + ".p_1"
            newCtgRight = newCtg + ".p_2"

        assert (lastOrderIndexOnLeft < len(self.relationChains[newCtg]))
        # insert all relations up to lastOrderIndexOnLeft to "newCtgLeft"
        # their indices don't need to be changed
        for relation in self.relationChains[newCtg][:lastOrderIndexOnLeft + 1]:
            relation.block.newCtg = newCtgLeft
            self.relationChains[newCtgLeft].append(relation)

        # insert all relations after otherHapOrderIndex to "otherHapNewCtgRight"
        for relation in self.relationChains[newCtg][lastOrderIndexOnLeft + 1:]:
            relation.block.newCtg = newCtgRight
            self.relationChains[newCtgRight].append(relation)

        # set the indices of all the blocks in the new contig on the right side
        for i, relation in enumerate(self.relationChains[newCtgRight]):
            relation.block.orderIndex = i

        # delete the parent contig
        del self.relationChains[newCtg]


    def addMisAssemblyToBeginning(self, block, misAssemblyType, misAssemblyLength):
        blockLength = block.origEnd - block.origStart + 1
        if blockLength <= misAssemblyLength:
            block.addMisAssemblyBlockList(BlockList([(1, blockLength, misAssemblyType)]))
        else:
            block.addMisAssemblyBlockList(BlockList([(1, misAssemblyLength, misAssemblyType)]))

    def addMisAssemblyToEnd(self, block, misAssemblyType, misAssemblyLength):
        blockLength = block.origEnd - block.origStart + 1
        if blockLength <= misAssemblyLength:
            block.addMisAssemblyBlockList(BlockList([(1, blockLength, misAssemblyType)]))
        else:
            block.addMisAssemblyBlockList(BlockList([(blockLength - misAssemblyLength + 1, blockLength, misAssemblyType)]))


    def induceMisjoinMisAssembly(self, newCtg_1, orderIndex_1, loc_1, newCtg_2, orderIndex_2, loc_2, misjoinEffectWindowLength):

        # In contrast to other misassembly generating methods
        # this methods needs two relation chains since
        # a misjoin is an event that joins two segments of the genome
        # that are far from each other. For the ease of implementation
        # this method takes two distinct chains, splits and joins them in the given
        # locations. Misjoin could happen in single chain but that is not implemented
        # yet
        # The variables in this method either have "_1" or "_2" suffix
        # Variables that end with "_1" are related to the first relation
        # chain and the ones with "_2" are related to the second chain.
        relationToSplit_1 = self.relationChains[newCtg_1][orderIndex_1]
        relationToSplit_2 = self.relationChains[newCtg_2][orderIndex_2]

        # get the order index and the name of the new contigs
        # of the homologous blocks
        otherHapOrderIndex_1 = relationToSplit_1.homologousBlock.orderIndex
        otherHapNewCtg_1 = relationToSplit_1.homologousBlock.newCtg

        otherHapOrderIndex_2 = relationToSplit_2.homologousBlock.orderIndex
        otherHapNewCtg_2 = relationToSplit_2.homologousBlock.newCtg

        # split the homology relation into two parts
        ref2querySplitRelations_1 = relationToSplit_1.splitIntoTwoParts(loc_1)
        ref2querySplitRelations_2 = relationToSplit_2.splitIntoTwoParts(loc_2)

        if ref2querySplitRelations_1 == None:
            return False
        if ref2querySplitRelations_2 == None:
            return False

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg_1].pop(orderIndex_1)
        self.relationChains[otherHapNewCtg_1].pop(otherHapOrderIndex_1)

        self.relationChains[newCtg_2].pop(orderIndex_2)
        self.relationChains[otherHapNewCtg_2].pop(otherHapOrderIndex_2)

        relationPart1_1 = ref2querySplitRelations_1[0]
        relationPart2_1 = ref2querySplitRelations_1[1]

        relationPart1_2 = ref2querySplitRelations_2[0]
        relationPart2_2 = ref2querySplitRelations_2[1]

        # add misassembly blocks to the homology blocks in the first chain with the suffix of "_1"
        rBlockPart1_1 = relationPart1_1.block
        rBlockPart2_1 = relationPart2_1.block
        self.addMisAssemblyToEnd(rBlockPart1_1, "Msj", misjoinEffectWindowLength)
        self.addMisAssemblyToBeginning(rBlockPart2_1, "Msj", misjoinEffectWindowLength)

        qBlockPart1_1 = relationPart1_1.homologousBlock
        qBlockPart2_1 = relationPart2_1.homologousBlock
        self.addMisAssemblyToEnd(qBlockPart1_1, "Msj", misjoinEffectWindowLength)
        self.addMisAssemblyToBeginning(qBlockPart2_1, "Msj", misjoinEffectWindowLength)

        # add misassembly blocks to the homology blocks in the first chain with the suffix of "_2"
        rBlockPart1_2 = relationPart1_2.block
        rBlockPart2_2 = relationPart2_2.block
        self.addMisAssemblyToEnd(rBlockPart1_2, "Msj", misjoinEffectWindowLength)
        self.addMisAssemblyToBeginning(rBlockPart2_2, "Msj", misjoinEffectWindowLength)

        qBlockPart1_2 = relationPart1_2.homologousBlock
        qBlockPart2_2 = relationPart2_2.homologousBlock
        self.addMisAssemblyToEnd(qBlockPart1_2, "Msj", misjoinEffectWindowLength)
        self.addMisAssemblyToBeginning(qBlockPart2_2, "Msj", misjoinEffectWindowLength)

        # create the equivalent list of relations from query to ref
        # these relations will show the same connections between blocks
        # but in the other way around
        query2refSplitRelations_1 = []
        query2refSplitRelations_2 = []
        for relation in ref2querySplitRelations_1:
            query2refRelation_1 = HomologyRelation(relation.homologousBlock,
                                                 relation.block,
                                                 None,
                                                 None)
            query2refSplitRelations_1.append(query2refRelation_1)

        for relation in ref2querySplitRelations_2:
            query2refRelation_2 = HomologyRelation(relation.homologousBlock,
                                                   relation.block,
                                                   None,
                                                   None)
            query2refSplitRelations_2.append(query2refRelation_2)

        ## Adjust query2ref relations
        # insert split relations to relation chain of the "otherHapNewCtg_1"
        for relation in query2refSplitRelations_1:
            self.relationChains[otherHapNewCtg_1].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by one
        for relation in self.relationChains[otherHapNewCtg_1][otherHapOrderIndex_1 + 2:]:
            relation.block.orderIndex += 1

        # insert split relations to relation chain of the "otherHapNewCtg_2"
        for relation in query2refSplitRelations_2:
            self.relationChains[otherHapNewCtg_2].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by one
        for relation in self.relationChains[otherHapNewCtg_2][otherHapOrderIndex_2 + 2:]:
            relation.block.orderIndex += 1

        # create two misjoins by swapping the links between ref2query relations
        # the old relation is already popped out above
        relationsPrePart1_1 = self.relationChains[newCtg_1][:orderIndex_1]
        relationsPostPart2_1 = self.relationChains[newCtg_1][orderIndex_1:]


        relationsPrePart1_2 = self.relationChains[newCtg_2][:orderIndex_2]
        relationsPostPart2_2 = self.relationChains[newCtg_2][orderIndex_2:]

        newRelations_1 = relationsPrePart1_1
        newRelations_1.extend([relationPart1_1, relationPart2_2])
        newRelations_1.extend(relationsPostPart2_2)

        newRelations_2 = relationsPrePart1_2
        newRelations_2.extend([relationPart1_2, relationPart2_1])
        newRelations_2.extend(relationsPostPart2_1)

        for i, relation in enumerate(newRelations_1):
            relation.block.orderIndex = i

        for i, relation in enumerate(newRelations_2):
            relation.block.orderIndex = i

        # rename new contigs
        newCtgRenamed_1 = newCtg_1 + ".Msj_" + newCtg_2
        newCtgRenamed_2 = newCtg_2 + ".Msj_" + newCtg_1

        # remove old chains
        # add new chains
        del self.relationChains[newCtg_1]
        self.relationChains[newCtgRenamed_1] = newRelations_1

        del self.relationChains[newCtg_2]
        self.relationChains[newCtgRenamed_2] = newRelations_2

        # rename newCtg attributes
        for relation in self.relationChains[newCtgRenamed_1]:
            relation.block.newCtg = newCtgRenamed_1

        for relation in self.relationChains[newCtgRenamed_2]:
            relation.block.newCtg = newCtgRenamed_2

        # get indices of previous contigs
        # to set sampling weights to zero
        newCtgIndex_1 = self.newCtgToIndexForSampling[newCtg_1]
        newCtgIndex_2 = self.newCtgToIndexForSampling[newCtg_2]

        # set the sampling weights of the previous chains to zero
        # since they no longer exist
        for annot, weightList in self.newCtgAnnotationWeightsForSampling.items():
            weightList[newCtgIndex_1] = 0
            weightList[newCtgIndex_2] = 0
        # for misjoins
        for annot, weightList in self.newCtgAnnotationWeightsForSamplingMisjoin.items():
            weightList[newCtgIndex_1] = 0
            weightList[newCtgIndex_2] = 0

        # add the names of the new misjoined contigs to the contig list for sampling
        self.newCtgListForSampling.extend([newCtgRenamed_1, newCtgRenamed_2])
        self.newCtgToIndexForSampling[newCtgRenamed_1] = len(self.newCtgListForSampling) - 2
        self.newCtgToIndexForSampling[newCtgRenamed_2] = len(self.newCtgListForSampling) - 1

        # add a new contig weight corresponding to newCtgRenamed_1 for sampling misassembly intervals
        for annot, weightList in self.newCtgAnnotationWeightsForSampling.items():
            weight = 0
            for relation in self.relationChains[newCtgRenamed_1]:
                weight += relation.block.annotationStartTotalLengthsForSampling[annot]
            weightList.append(weight)

        # add a new contig weight corresponding to newCtgRenamed_2 for sampling misassembly intervals
        for annot, weightList in self.newCtgAnnotationWeightsForSampling.items():
            weight = 0
            for relation in self.relationChains[newCtgRenamed_2]:
                weight += relation.block.annotationStartTotalLengthsForSampling[annot]
            weightList.append(weight)

        # add a new contig weight corresponding to newCtgRenamed_1 for sampling misjoin locations
        for annot, weightList in self.newCtgAnnotationWeightsForSamplingMisjoin.items():
            weight = 0
            for relation in self.relationChains[newCtgRenamed_1]:
                weight += relation.block.annotationTotalLengthsForSamplingMisjoin[annot]
            weightList.append(weight)

        # add a new contig weight corresponding to newCtgRenamed_2 for sampling misjoin locations
        for annot, weightList in self.newCtgAnnotationWeightsForSamplingMisjoin.items():
            weight = 0
            for relation in self.relationChains[newCtgRenamed_2]:
                weight += relation.block.annotationTotalLengthsForSamplingMisjoin[annot]
            weightList.append(weight)

        return True


    def induceCollapseMisAssembly(self, newCtg, orderIndex, collapseStart, collapseEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(collapseStart, collapseEnd)

        if ref2querySplitRelations == None:
            return False

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)


        # the homologous block, which is qBlockPart2, should be removed
        # from the middle relation. The alignment should also be removed
        # from this relation
        relationPart2 = ref2querySplitRelations[1]
        relationPart2.homologousBlock = None
        relationPart2.alignment = None

        rBlockPart2 = relationPart2.block
        rBlockPart2.addMisAssemblyBlockList(BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1, "Col_Del")]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationBlocksForSampling()

        # since qBlockPart2 is going to be ignored, orderIndex of qBlockPart1/3 should be adjusted
        qBlockPart1 = ref2querySplitRelations[0].homologousBlock
        qBlockPart3 = ref2querySplitRelations[2].homologousBlock
        if relationToSplit.alignment.orientation == '+':
            qBlockPart3.orderIndex -= 1
        else:
            qBlockPart1.orderIndex -= 1

        # create the equivalent list of relations from query to ref
        # these relations will show the same connections between blocks
        # but in the other way around
        query2refSplitRelations = []
        if relationToSplit.alignment.orientation == '+':
            for relation in [ref2querySplitRelations[0], ref2querySplitRelations[2]]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)
        else:
            for relation in [ref2querySplitRelations[2], ref2querySplitRelations[0]]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)

        ## ref to query ##

        # insert split relations to relation chain of the "newCtg"
        for relation in ref2querySplitRelations:
            self.relationChains[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        ## query to ref ##

        # insert split relations to relation chain of the "otherHapNewCtg"
        for relation in query2refSplitRelations:
            self.relationChains[otherHapNewCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by one
        for relation in self.relationChains[otherHapNewCtg][otherHapOrderIndex + 2:]:
            relation.block.orderIndex += 1

        # break the query contigs into two parts
        # one part ends right before the collapsed block
        # the other part starts right after the collapsed block
        self.breakIntoTwoContigs(newCtg=otherHapNewCtg,
                                 lastOrderIndexOnLeft=otherHapOrderIndex)

        self.updateNewCtgAnnotationWeightsForSampling(newCtg, relationToSplit, ref2querySplitRelations)

        return True

    def induceDuplicationMisAssembly(self, newCtg, orderIndex, duplicationStart, duplicationEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(duplicationStart, duplicationEnd)

        if ref2querySplitRelations == None:
            return False

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)


        # the block that has to be duplicated
        rBlockPart2 = ref2querySplitRelations[1].block
        rBlockPart2.addMisAssemblyBlockList(BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1, "Dup")]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationBlocksForSampling()

        # falsely duplicated block, this is the duplication of the middle part of rBlock
        newDupCtg = rBlockPart2.origCtg + f".Dup_{rBlockPart2.origStart}_{rBlockPart2.origEnd}"
        rBlockPart2Dup = HomologyBlock(rBlockPart2.origCtg,
                                       rBlockPart2.origStart,
                                       rBlockPart2.origEnd,
                                       '+',
                                       newDupCtg,
                                       0)
        rBlockPart2Dup.annotationBlockLists = deepcopy(rBlockPart2.annotationBlockLists)
        rBlockPart2Dup.addMisAssemblyBlockList(BlockList([(1, rBlockPart2Dup.origEnd - rBlockPart2Dup.origStart + 1, "Dup")]))
        rBlockPart2Dup.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2Dup.clearAnnotationBlocksForSampling()

        self.relationChains[newDupCtg] = [HomologyRelation(rBlockPart2Dup,
                                                          None,
                                                          None,
                                                          None)]

        # create the equivalent list of relations from query to ref
        # these relations will show the same connections between blocks
        # but in the other way around
        query2refSplitRelations = []
        if relationToSplit.alignment.orientation == '+':
            for relation in ref2querySplitRelations:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)
        else:
            for relation in ref2querySplitRelations[::-1]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)

        # insert split relations to relation chain of the "newCtg"
        for relation in ref2querySplitRelations:
            self.relationChains[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        # insert split relations to relation chain of the "otherHapNewCtg"
        for relation in query2refSplitRelations:
            self.relationChains[otherHapNewCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[otherHapNewCtg][otherHapOrderIndex + 3:]:
            relation.block.orderIndex += 2

        self.updateNewCtgAnnotationWeightsForSampling(newCtg, relationToSplit, ref2querySplitRelations)

        return True


    def induceBaseErrorMisAssembly(self, newCtg, orderIndex, errorStart, errorEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(errorStart, errorEnd)

        if ref2querySplitRelations == None:
            return False

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)


        # the ref block that has to be contaminated with base errors
        rBlockPart2 = ref2querySplitRelations[1].block
        rBlockPart2.addMisAssemblyBlockList(BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1, "Err")]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationBlocksForSampling()

        # the query block that is supposed to be collapsed since the
        # ref haplotype is highly erroneous
        qBlockPart2 = ref2querySplitRelations[1].homologousBlock
        qBlockPart2.addMisAssemblyBlockList(BlockList([(1, qBlockPart2.origEnd - qBlockPart2.origStart + 1, "Col_Err")]))
        qBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        qBlockPart2.clearAnnotationBlocksForSampling()

        # the erroneous block will not have the previous alignment
        ref2querySplitRelations[1].alignment = None

        # create the equivalent list of relations from query to ref
        # these relations will show the same connections between blocks
        # but in the other way around
        query2refSplitRelations = []
        if relationToSplit.alignment.orientation == '+':
            for relation in ref2querySplitRelations:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)
        else:
            for relation in ref2querySplitRelations[::-1]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refSplitRelations.append(query2refRelation)

        # insert split relations to relation chain of the "newCtg"
        for relation in ref2querySplitRelations:
            self.relationChains[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        # insert split relations to relation chain of the "otherHapNewCtg"
        for relation in query2refSplitRelations:
            self.relationChains[otherHapNewCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in self.relationChains[otherHapNewCtg][otherHapOrderIndex + 3:]:
            relation.block.orderIndex += 2

        self.updateNewCtgAnnotationWeightsForSampling(newCtg, relationToSplit, ref2querySplitRelations)

        return True

    def updateNewCtgAnnotationWeightsForSampling(self, newCtg, parentRelation, childRelations):
        
        assert(newCtg in self.newCtgToIndexForSampling)
        newCtgIndex = self.newCtgToIndexForSampling[newCtg]

        # subtract the old weights related to the parent block
        parentRefBlock = parentRelation.block
        for annot, total in parentRefBlock.annotationStartTotalLengthsForSampling.items():
            # if "updateAnnotationBlocksForSampling" is not called
            # then the weight list would be empty
            if 0 < len(self.newCtgAnnotationWeightsForSampling[annot]):
                self.newCtgAnnotationWeightsForSampling[annot][newCtgIndex] -= total

        # add the new weights related to the child blocks
        for relation in childRelations:
            for annot, total in relation.block.annotationStartTotalLengthsForSampling.items():
                # if "updateAnnotationBlocksForSampling" is not called
                # then the weight list would be empty
                if 0 < len(self.newCtgAnnotationWeightsForSampling[annot]):
                    self.newCtgAnnotationWeightsForSampling[annot][newCtgIndex] += total

    def getListOfSamplingLengths(self, newCtg, annotation):
        relations = self.relationChains[newCtg]
        lengths = []
        for relation in relations:
            lengths.append(relation.block.annotationStartTotalLengthsForSampling[annotation])
        return  lengths

    def getTotalSamplingLength(self, newCtg, annotation):
        return sum(self.getListOfSamplingLengths(newCtg, annotation))


    def updateAnnotationBlocksForSampling(self, annotations, misAssemblyLength, minOverlapRatioWithEachAnnotation, minMarginLength):
        """
        This function has to be called each time the misassembly length is changed
        """
        self.misAssemblyLength = misAssemblyLength
        for newCtg, relations in self.relationChains.items():
            for relation in relations:
                # Sampling will happen only in the blocks containing
                # a single 1-to-1 alignment with no previously created misassembly
                # And also only in the reference/hap1 blocks
                # (TODO:Maybe enabling misassembly creation from hap2 later)
                if relation.homologousBlock is not None and \
                        relation.block.containsMisAssembly is False and \
                        relation.alignment is not None:

                    relation.block.setSamplingLengthAttributes(misAssemblyLength,
                                                               minOverlapRatioWithEachAnnotation,
                                                               minMarginLength)
                    relation.block.updateAnnotationBlocksForSampling()
                else:
                    relation.block.clearAnnotationBlocksForSampling()

        # for sampling intervals of "Sw", "Err", "Dup", and "Col"
        for annotation in annotations:
            newCtgWeights = []
            for newCtg in self.newCtgListForSampling:
                newCtgWeights.append(self.getTotalSamplingLength(newCtg, annotation))
            self.newCtgAnnotationWeightsForSampling[annotation] =  newCtgWeights

        # for sampling intervals of "Msj"
        for annotation in annotations:
            newCtgWeights = []
            for newCtg in self.newCtgListForSampling:
                newCtgWeights.append(self.getTotalSamplingMisjoinLength(newCtg, annotation))
            self.newCtgAnnotationWeightsForSamplingMisjoin[annotation] =  newCtgWeights

    def getWeightedRandomNewCtgForSampling(self, annotation):
        """
        Given the annotation name select one new contig randomly by taking the
        total length of sampling regions as the sampling weight for each new contig

        :param annotation: The annotation name
        :return: The randomly selected new contig
        """
        #print(annotation, self.newCtgAnnotationWeightsForSampling)
        if sum(self.newCtgAnnotationWeightsForSampling[annotation]) == 0:
            return None
        selectedNewCtg = random.choices(self.newCtgListForSampling,
                                        weights=self.newCtgAnnotationWeightsForSampling[annotation],
                                        k=1)[0]
        return selectedNewCtg

    def getWeightedRandomOrderIndexForSampling(self, newCtg, annotation):
        """
        Given the annotation and the new contig name select one relation index randomly by taking the
        total length of sampling regions as the sampling weight for each relation

        :param annotation: The annotation name
        :param newCtg: The name of the new contig
        :return: The randomly selected order index
        """
        weights = self.getListOfSamplingLengths(newCtg, annotation)
        orderIndex = random.choices(np.arange(len(self.relationChains[newCtg])),
                                    weights=weights,
                                    k=1)[0]
        return orderIndex

    def getRandomIntervalFromRelationBlock(self, newCtg, annotation, orderIndex, misAssemblyLength):
        """
        :param newCtg: The name of the new contig
        :param annotation: The annotation name
        :param orderIndex:  The order index of the selected relation/block in its relation chain
        :param misAssemblyLength: Misassembly length
        :return: The randomly selected start and end
        """
        block = self.relationChains[newCtg][orderIndex].block
        start, end = block.sampleMisAssemblyInterval(annotation, misAssemblyLength)
        return start, end

    def getRandomMisAssemblyInterval(self, annotation, misAssemblyLength):
        """
        This method can be called to obtain a random location for inducing a misassembly
        overlapping with the annotation of interest

        Note that the amount of overlap was determined while calling "updateAnnotationBlocksForSampling"

        :param annotation: The annotation name
        :param misAssemblyLength: Misassembly length
        :return: a randomly selected new contig name, index of relation, start and end coordinates
        """
        newCtg = self.getWeightedRandomNewCtgForSampling(annotation)
        if newCtg == None:
            return None, None, None, None
        orderIndex = self.getWeightedRandomOrderIndexForSampling(newCtg, annotation)
        start, end = self.getRandomIntervalFromRelationBlock(newCtg, annotation, orderIndex, misAssemblyLength)
        return newCtg, orderIndex, start, end


    # functions for sampling misjoin locations
    def getListOfSamplingMisjoinLengths(self, newCtg, annotation):
        relations = self.relationChains[newCtg]
        lengths = []
        for relation in relations:
            lengths.append(relation.block.annotationTotalLengthsForSamplingMisjoin[annotation])
        return  lengths

    def getTotalSamplingMisjoinLength(self, newCtg, annotation):
        return sum(self.getListOfSamplingMisjoinLengths(newCtg, annotation))

    def getWeightedRandomNewCtgForSamplingMisjoin(self, annotation, newCtgListToExclude):
        """
        Given the annotation name select one new contig randomly by taking the
        total length of sampling (for misjoin) regions as the sampling weight for each new contig

        :param annotation: The annotation name
        :param newCtgListToExclude: a list of new contigs to exclude from the sampling process
        :return: The randomly selected new contig
        """
        #print(annotation, self.newCtgAnnotationWeightsForSampling)
        newCtgIndicesToExlucde = [self.newCtgToIndexForSampling[c] for c in newCtgListToExclude]
        weights = np.array(self.newCtgAnnotationWeightsForSamplingMisjoin[annotation])
        newCtgList = np.array(self.newCtgListForSampling)
        # exclude contigs
        weights = np.delete(weights, newCtgIndicesToExlucde)
        newCtgList = np.delete(newCtgList, newCtgIndicesToExlucde)
        if sum(weights) == 0:
            return None
        selectedNewCtg = random.choices(newCtgList,
                                        weights=weights,
                                        k=1)[0]
        return selectedNewCtg

    def getWeightedRandomOrderIndexForSamplingMisjoin(self, newCtg, annotation):
        """
        Given the annotation and the new contig name select one relation index randomly by taking the
        total length of sampling (for misjoin) regions as the sampling weight for each relation

        :param annotation: The annotation name
        :param newCtg: The name of the new contig
        :return: The randomly selected order index
        """
        weights = self.getListOfSamplingMisjoinLengths(newCtg, annotation)
        orderIndex = random.choices(np.arange(len(self.relationChains[newCtg])),
                                    weights=weights,
                                    k=1)[0]
        return orderIndex

    def getRandomLocationFromRelationBlock(self, newCtg, annotation, orderIndex):
        """
        :param newCtg: The name of the new contig
        :param annotation: The annotation name
        :param orderIndex:  The order index of the selected relation/block in its relation chain
        :return: The randomly selected location
        """
        block = self.relationChains[newCtg][orderIndex].block
        loc = block.sampleMisjoinLocation(annotation)
        return loc

    def getRandomMisjoinLocation(self, annotation, newCtgListToExclude):
        """
        This method can be called to obtain a random location for inducing a misjoin
        overlapping with the annotation of interest

        :param annotation: The annotation name
        :param newCtgListToExclude: a list of new contigs to exclude from the sampling process
        :return: a randomly selected new contig name, index of relation, and location
        """
        newCtg = self.getWeightedRandomNewCtgForSamplingMisjoin(annotation, newCtgListToExclude)
        if newCtg == None:
            return None, None, None
        orderIndex = self.getWeightedRandomOrderIndexForSamplingMisjoin(newCtg, annotation)
        loc = self.getRandomLocationFromRelationBlock(newCtg, annotation, orderIndex)
        return newCtg, orderIndex, loc



    def getNewCtgSequence(self, newCtg, origCtgSequences, singleBaseErrorRate):
        newCtgSeqList = []
        for relation in self.relationChains[newCtg]:
            blockSeq = relation.block.getSequence(origCtgSequences, singleBaseErrorRate)
            newCtgSeqList.append(blockSeq)
        return  "".join(newCtgSeqList)

    def yieldNewCtgSequences(self, origCtgSequences, singleBaseErrorRate):
        newCtgList = sorted(list(self.relationChains.keys()))
        for newCtg in newCtgList:
            yield newCtg, self.getNewCtgSequence(newCtg, origCtgSequences, singleBaseErrorRate)

    def getTotalCountOfLongerBlocks(self, annotation, minBlockSize, onlyRefInHomology=False):
        """
        It counts the number of blocks longer than minimum block size for one annotation and
        returns the total count
        :param annotation: annotation name
        :param minBlockSize: Minimum block size
        :param onlyRefInHomology: If True it will consider only ref blocks with 1-to-1 mappings
        :return: totalCount
        """
        totalCount = 0
        for newCtg, relations in self.relationChains.items():
            for relation in relations:
                # if there is no homology for this block, skip it
                if (relation.homologousBlock is None and
                        onlyRefInHomology):
                    continue
                # if there is a homology but the original contig of this block is not from reference, skip it
                if (relation.homologousBlock is not None and
                        onlyRefInHomology and
                        relation.block.origCtg not in self.origRefContigNames):
                    continue
                for block in relation.block.annotationBlockLists[annotation].blocks:
                    if minBlockSize < (block[1] - block[0]):
                        totalCount += 1

        return totalCount

    def getTotalCountOfLongerBlocksForAllAnnotations(self, annotations, minBlockSizes, onlyRefInHomology=False):
        """
        Run getTotalCountOfLongerBlocks for each given annotation and minimum block size
        :param annotations: A list of annotation names
        :param minBlockSizes: A list of minimum block sizes
        :param onlyRefInHomology: If True it will consider only ref blocks in 1-to-1 mappings
        :return: A dictionary with annotation names as keys and count lists as values. If the input blockSizes
                 is None each value in the output dictionary will be an integer instead of list. It can be used for
                 cases when we don't want to filter annotation blocks.
        """
        if minBlockSizes is None:
            totalCounts = defaultdict(int)
        else:
            totalCounts = defaultdict(list)
        for annotation in annotations:
            if minBlockSizes is None:
                totalCounts[annotation] = self.getTotalCountOfLongerBlocks(annotation, 0, onlyRefInHomology)
            else:
                for minBlockSize in minBlockSizes:
                    totalCounts[annotation].append(self.getTotalCountOfLongerBlocks(annotation, minBlockSize, onlyRefInHomology))
        return totalCounts

    def getTotalLengthOfLongerBlocks(self, annotation, minBlockSize, onlyRefInHomology=False):
        """
        It computes the total length of blocks longer than minimum block size for one annotation and
        returns the total length
        :param annotation: annotation name
        :param minBlockSize: Minimum block size
        :param onlyRefInHomology: If True it will consider only ref blocks with 1-to-1 mappings
        :return: totalLength
        """
        totalLength = 0
        for newCtg, relations in self.relationChains.items():
            for relation in relations:
                # if there is no homology for this block, skip it
                if (relation.homologousBlock is None and
                        onlyRefInHomology):
                    continue
                # if there is a homology but the original contig of this block is not from reference, skip it
                if (relation.homologousBlock is not None and
                        onlyRefInHomology and
                        relation.block.origCtg not in self.origRefContigNames):
                    continue
                for block in relation.block.annotationBlockLists[annotation].blocks:
                    if minBlockSize < (block[1] - block[0]):
                        totalLength += block[1] - block[0]
        return totalLength

    def getTotalLengthOfLongerBlocksForAllAnnotations(self, annotations, minBlockSizes, onlyRefInHomology=False):
        """
        Run getTotalLengthOfLongerBlocks for each given annotation and minimum block size
        :param annotations: A list of annotation names
        :param minBlockSizes: A list of block sizes (it can be None)
        :param onlyRefInHomology: If True it will consider only ref blocks in 1-to-1 mappings
        :return: A dictionary with annotation names as keys and total length lists as values. If minBlockSizes
                 is None each value in the output dictionary will be an integer instead of list. It can be used for cases
                 when we don't want to filter annotation blocks.
        """
        if minBlockSizes is None:
            totalLengths = defaultdict(int)
        else:
            totalLengths = defaultdict(list)
        for annotation in annotations:
            if minBlockSizes is None:
                totalLengths[annotation] = self.getTotalLengthOfLongerBlocks(annotation, 0, onlyRefInHomology)
            else:
                for minBlockSize in minBlockSizes:
                    totalLengths[annotation].append(self.getTotalLengthOfLongerBlocks(annotation, minBlockSize, onlyRefInHomology))
        return totalLengths

    def getLowerBoundOnNumberOfMisassemblies(self, annotation, misAssemblySize, marginSize):
        """
        It computes the minimum number of misassemblies of a specific length that can be created in the given annotation
        :param annotation:  annotation name
        :param misAssemblySize: misassembly size
        :return: lowerBound
        """
        lowerBound = 0
        for newCtg, relations in self.relationChains.items():
            for relation in relations:
                # if there is no homology for this block, skip it
                if relation.homologousBlock is None:
                    continue
                # if there is a homology but the original contig of this block is not from reference, skip it
                if relation.homologousBlock is not None and relation.block.origCtg not in self.origRefContigNames:
                    continue
                for block in relation.block.annotationBlockLists[annotation].blocks:
                    if misAssemblySize < (block[1] - block[0]):
                        # X = block size
                        # L = effective misassembly size
                        ## n = [log2( X/L + 1)] - 1
                        # how many times we can split the whole block until we get a block smaller than the
                        # misassembly size
                        numberOfSplits = np.floor(np.log2((block[1] - block[0]) / misAssemblySize + 1)) - 1
                        numberOfSplits = 0 if numberOfSplits < 0 else numberOfSplits
                        lowerBound += np.pow(2, numberOfSplits) + 1
        return lowerBound


    def writeNewContigsToFasta(self, origCtgSequences, fastaPath, singleBaseErrorRate):
        handle = open(fastaPath, "w")
        writer = FastaWriter(handle)
        for newCtg, newSeq in self.yieldNewCtgSequences(origCtgSequences, singleBaseErrorRate):
            record = SeqRecord(Seq(newSeq),
                               id = newCtg,
                               description = "")
            writer.write_record(record)
        handle.close()


    def getNewCtgCoordinatesOfAnnotation(self, newCtg, annotation):
        """
        :param newCtg: The name of the new contig
        :param annotation: The name of the annotation
        :return: a BlockList that contains the coordinates of "annotation" in "newCtg"
        """
        newCtgAnnotationBlockList = BlockList()
        newStart = 1
        for relation in self.relationChains[newCtg]:
            blockLen = relation.block.origEnd - relation.block.origStart + 1
            annotationBlockList = relation.block.annotationBlockLists[annotation].copy()
            # shifting coors to make them to be w.r.t. the start of the whole new contig
            annotationBlockList.shift(newStart - 1, newStart, newStart + blockLen - 1, inplace=True)
            # add coors of the annotation of this block to the block list for the whole new contig
            newCtgAnnotationBlockList.extend(annotationBlockList)
            # update the start coordinate for the next block
            newStart += blockLen
        return newCtgAnnotationBlockList


    def getNewCtgCoordinatesOfMisAssemblies(self, newCtg):
        """
        :param newCtg: The name of the new contig
        :return: a BlockList that contains the coordinates of all misassemblies in "newCtg"
        """
        newCtgMisAssembliesBlockList = BlockList()
        newStart = 1
        for relation in self.relationChains[newCtg]:
            blockLen = relation.block.origEnd - relation.block.origStart + 1
            misAssemblyBlockList = relation.block.misAssemblyBlockList.copy()
            # shift coors to make them to be w.r.t. the start of the whole new contig
            misAssemblyBlockList.shift(newStart - 1, newStart, newStart + blockLen - 1, inplace=True)
            # add coors of the annotation of this block to the block list for the whole new contig
            newCtgMisAssembliesBlockList.extend(misAssemblyBlockList)
            # update the start coordinate for the next block
            newStart += blockLen
        return newCtgMisAssembliesBlockList

    def writeAnnotationCoordinatesToBed(self, annotation, bedPath):
        """
        Write the coordinates of an annotation in the whole (falsified) assembly
        to a bed file
        :param annotation: The name of the annotation
        :param bedPath: The path to the output bed file
        """
        newCtgList = sorted(list(self.relationChains.keys()))
        with open(bedPath, "w") as f:
            for newCtg in newCtgList:
                annotationBlockList = self.getNewCtgCoordinatesOfAnnotation(newCtg, annotation)
                for block in annotationBlockList.blocks:
                    # start should be 0-based in bed
                    f.write(f"{newCtg}\t{block[0] - 1}\t{block[1]}\n")

    def writeMisAssemblyCoordinatesToBed(self, bedPath):
        """
        Write the coordinates of mis-assemblies in the whole (falsified) assembly
        to a bed file
        :param bedPath: The path to the output bed file
        """
        newCtgList = sorted(list(self.relationChains.keys()))
        with open(bedPath, "w") as f:
            for newCtg in newCtgList:
                misAssemblyBlockList = self.getNewCtgCoordinatesOfMisAssemblies(newCtg)
                for block in misAssemblyBlockList.blocks:
                    # start should be 0-based in bed
                    f.write(f"{newCtg}\t{block[0] - 1}\t{block[1]}\t{block[2]}\n")
















