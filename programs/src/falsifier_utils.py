from block_utils import *
from collections import defaultdict
from copy import deepcopy
import random
import numpy as np
import re

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
        self.misAssemblyBlockLists = defaultdict(BlockList)

        # attributes related to sampling locations
        # they will be updated with "setSamplingLengthAttributes" or "copySamplingLengthAttributes"
        self.annotationStartBlockListsForSampling = defaultdict(BlockList)
        self.annotationStartBlockLengthsForSampling = defaultdict(list)
        self.annotationStartTotalLengthsForSampling = {}
        self.misAssemblyLength = 0
        self.minOverlapRatioWithEachAnnotation = 0.5
        self.minMarginLength = 10000
        self.containsMisAssembly = False

    def getSequence(self, origCtgSequences):
        origCtgSeq = origCtgSequences[self.origCtg]
        # note that start and end are 1-based closed
        blockSeq = origCtgSeq[self.origStart - 1: self.origEnd]
        if self.origStrand == '-':
            blockSeq = reverseComplement(blockSeq)
        return blockSeq

    def reverseAllBlockLists(self):
        """
            Make the coordinates of all the block lists in this homology block w.r.t to the end of whole block
        """
        blockLength = self.origEnd - self.origStart + 1
        for annotation in self.annotationBlockLists:
            self.annotationBlockLists[annotation].reverse(blockLength, inplace=True)
            self.annotationStartBlockListsForSampling[annotation].reverse(blockLength, inplace=True)
            self.annotationStartBlockLengthsForSampling[annotation].reverse()
        for name in self.misAssemblyBlockLists:
            self.misAssemblyBlockLists[name].reverse(blockLength, inplace=True)

    def clearAnnotationStartBlocksForSampling(self):
        for annotation in self.annotationBlockLists:
            self.annotationStartBlockListsForSampling[annotation] = BlockList([])
            self.annotationStartBlockLengthsForSampling[annotation] = []
            self.annotationStartTotalLengthsForSampling[annotation] = 0

    def addMisAssemblyBlockList(self, name, blockList):
        """

        :param name: The name of the misassembly to add [either "Err", "Dup" or "Col"]
        :param blockList: A BlockList containing the coordinates of the misassembled part
                          (Note that coordinates should be 1-based and relative to
                           the start position of this block rather than the original contig.
                           The "shift" method of "BlockList" can be used for shifting
                           coordinates prior to adding blockList)
        """
        self.misAssemblyBlockLists[name] = blockList.copy()

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

    def updateAnnotationStartLocationsForSampling(self):
        """
        This method performs two main truncations to narrow down the annotation blocks where the start locations of
        misassemblies are going to be sampled from:
            - Truncate each continuous annotation block from the right side as long as "lengthToTruncateFromEnd"
            - Truncate each continuous annotation block from both sides as far as it won't have
               overlap with the margins of the whole block.
               The length of the right margin should be greater than the misassembly length to make sure that
               the whole misassembled segment is fully within the block
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
        wholeBlockMarginFromRight = self.misAssemblyLength + self.minMarginLength

        wholeBlockWithoutMargins = BlockList([(1, self.origEnd - self.origStart + 1)]).truncateFromLeft(wholeBlockMarginFromLeft, inplace=False)
        wholeBlockWithoutMargins.truncateFromRight(wholeBlockMarginFromRight - 1, inplace=True)

        for name in self.annotationBlockLists:
            self.annotationStartBlockListsForSampling[name] = self.annotationBlockLists[name].copy()

            # remove blocks shorter than minimum overlap 
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


    def sampleMisAssemblyInterval(self, name, misAssemblyLength):
        # self.misAssemblyLength was set the last time
        # "updateAnnotationStartLocationsForSampling" was invoked
        # If the previous value does match the current one
        # return None since it is essential to update the
        # annotation start locations for sampling based on
        # the correct misassembly length
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
        :param start: 1-based location of the parent block which is the first base in this block
        :param end: 1-based location of the parent block which is the last base in this block
        """
        for name, blockList in parentBlock.misAssemblyBlockLists.items():
            subsetBlockList = blockList.intersect(BlockList([(start, end)]), inplace=False)
            subsetBlockList.shift( -(start - 1), minCoordinate = 1, maxCoordinate = end - start + 1, inplace = True)
            self.addMisAssemblyBlockList(name, subsetBlockList)



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
        rBlockPart1.updateAnnotationStartLocationsForSampling()

        rBlockPart2 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[1][0],
                                    projectionsOrigCoor[1][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 1)
        rBlockPart2.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])
        rBlockPart2.copySamplingLengthAttributes(rBlock)
        rBlockPart2.updateAnnotationStartLocationsForSampling()

        rBlockPart3 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[2][0],
                                    projectionsOrigCoor[2][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 2)
        rBlockPart3.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])
        rBlockPart3.extractMisAssembliesFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])
        rBlockPart3.copySamplingLengthAttributes(rBlock)
        rBlockPart3.updateAnnotationStartLocationsForSampling()

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
        qBlockPart1.clearAnnotationStartBlocksForSampling()
        #qBlockPart1.copySamplingLengthAttributes(qBlock)
        #qBlockPart1.updateAnnotationStartLocationsForSampling()

        qBlockPart2 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[1][2],
                                    projectionsOrigCoor[1][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart2)
        qBlockPart2.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])
        qBlockPart2.clearAnnotationStartBlocksForSampling()
        #qBlockPart2.copySamplingLengthAttributes(qBlock)
        #qBlockPart2.updateAnnotationStartLocationsForSampling()

        qBlockPart3 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[2][2],
                                    projectionsOrigCoor[2][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart3)
        qBlockPart3.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[2][2], projectionsRelCoor[2][3])
        qBlockPart3.extractMisAssembliesFromParentBlock(qBlock, projectionsRelCoor[2][2], projectionsRelCoor[2][3])
        qBlockPart3.clearAnnotationStartBlocksForSampling()
        #qBlockPart3.copySamplingLengthAttributes(qBlock)
        #qBlockPart3.updateAnnotationStartLocationsForSampling()

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
        # it will be filled by calling "updateAnnotationStartLocationsForSampling"
        self.newCtgAnnotationWeightsForSampling = defaultdict(list)
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

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(switchStart, switchEnd)

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

        # add misassembly with "Msj" label to both ends of the middle reference block
        rBlockPart2Length = rBlockPart2.origEnd - rBlockPart2.origStart + 1
        if rBlockPart2Length <= 2.5 * switchEffectWindowLength:
            rBlockPart2.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, rBlockPart2Length)]))
        else:
            rBlockPart2.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, switchEffectWindowLength),
                                                           (rBlockPart2Length - switchEffectWindowLength + 1, rBlockPart2Length)]))
        rBlockPart2.containsMisAssembly= True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationStartBlocksForSampling()

        # add misassembly with "Msj" label to both ends of the middle query block
        qBlockPart2Length = qBlockPart2.origEnd - qBlockPart2.origStart + 1
        if qBlockPart2Length <= 2.5 * switchEffectWindowLength:
            qBlockPart2.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, qBlockPart2Length)]))
        else:
            qBlockPart2.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, switchEffectWindowLength),
                                                           (qBlockPart2Length - switchEffectWindowLength + 1, qBlockPart2Length)]))
        qBlockPart2.containsMisAssembly= True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        qBlockPart2.clearAnnotationStartBlocksForSampling()


        # add misassembly with the "Msj" label to the end part of the left reference block
        rBlockPart1 = ref2querySplitRelations[0].block
        rBlockPart1Length = rBlockPart1.origEnd - rBlockPart1.origStart + 1
        if rBlockPart1Length <= switchEffectWindowLength:
            rBlockPart1.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, rBlockPart1Length)]))
        else:
            rBlockPart1.addMisAssemblyBlockList("Msj",
                                                BlockList([(rBlockPart1Length - switchEffectWindowLength + 1, rBlockPart1Length)]))


        # add misassembly with the "Msj" label to the beginning part of the right reference block
        rBlockPart3 = ref2querySplitRelations[2].block
        rBlockPart3Length = rBlockPart3.origEnd - rBlockPart3.origStart + 1
        if rBlockPart3Length <= switchEffectWindowLength:
            rBlockPart3.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, rBlockPart3Length)]))
        else:
            rBlockPart3.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, switchEffectWindowLength)]))


        # add misassembly with the "Msj" label to the left query block
        # it will be added to the end if orientation was positive
        # it will be added to the beginning if orientation was negative
        qBlockPart1 = ref2querySplitRelations[0].homologousBlock
        qBlockPart1Length = qBlockPart1.origEnd - qBlockPart1.origStart + 1
        if qBlockPart1Length <= switchEffectWindowLength:
            qBlockPart1.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, qBlockPart1Length)]))
        else:
            if relationToSplit.alignment.orientation == '+':
                qBlockPart1.addMisAssemblyBlockList("Msj",
                                                    BlockList([(qBlockPart1Length - switchEffectWindowLength + 1, qBlockPart1Length)]))
            else: # '-'
                qBlockPart1.addMisAssemblyBlockList("Msj",
                                                    BlockList([(1, switchEffectWindowLength)]))


        # add misassembly with the "Msj" label to the right query block
        # it will be added to the beginning if orientation was positive
        # it will be added to the end if orientation was negative
        qBlockPart3 = ref2querySplitRelations[2].homologousBlock
        qBlockPart3Length = qBlockPart3.origEnd - qBlockPart3.origStart + 1
        if qBlockPart3Length <= switchEffectWindowLength:
            qBlockPart3.addMisAssemblyBlockList("Msj",
                                                BlockList([(1, qBlockPart3Length)]))
        else:
            if relationToSplit.alignment.orientation == '+':
                qBlockPart3.addMisAssemblyBlockList("Msj",
                                                    BlockList([(1, switchEffectWindowLength)]))
            else: # '-'
                qBlockPart3.addMisAssemblyBlockList("Msj",
                                                    BlockList([(qBlockPart3Length - switchEffectWindowLength + 1, qBlockPart3Length)]))



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
        z = re.findall("(?<=_f_p)[0-9.]+$", newCtg)
        if len(z) == 1:
            newCtgLeft = newCtg + ".1"
            newCtgRight = newCtg + ".2"
        else:
            newCtgLeft = newCtg + "_p1"
            newCtgRight = newCtg + "_p2"

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

    def induceCollapseMisAssembly(self, newCtg, orderIndex, collapseStart, collapseEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(collapseStart, collapseEnd)

        # the homologous block, which is qBlockPart2, should be removed
        # from the middle relation. The alignment should also be removed
        # from this relation
        relationPart2 = ref2querySplitRelations[1]
        relationPart2.homologousBlock = None
        relationPart2.alignment = None

        rBlockPart2 = relationPart2.block
        rBlockPart2.addMisAssemblyBlockList("Col",
                                            BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1)]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationStartBlocksForSampling()

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

    def induceDuplicationMisAssembly(self, newCtg, orderIndex, duplicationStart, duplicationEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(duplicationStart, duplicationEnd)

        # the block that has to be duplicated
        rBlockPart2 = ref2querySplitRelations[1].block
        rBlockPart2.addMisAssemblyBlockList( "Dup",
                                                BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1)]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationStartBlocksForSampling()

        # falsely duplicated block, this is the duplication of the middle part of rBlock
        newDupCtg = rBlockPart2.origCtg + f"_Dup_{rBlockPart2.origStart}_{rBlockPart2.origEnd}"
        rBlockPart2Dup = HomologyBlock(rBlockPart2.origCtg,
                                       rBlockPart2.origStart,
                                       rBlockPart2.origEnd,
                                       '+',
                                       newDupCtg,
                                       0)
        rBlockPart2Dup.annotationBlockLists = deepcopy(rBlockPart2.annotationBlockLists)
        rBlockPart2Dup.addMisAssemblyBlockList( "Dup",
                                                BlockList([(1, rBlockPart2Dup.origEnd - rBlockPart2Dup.origStart + 1)]))
        rBlockPart2Dup.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2Dup.clearAnnotationStartBlocksForSampling()

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


    def induceBaseErrorMisAssembly(self, newCtg, orderIndex, errorStart, errorEnd):

        relationToSplit = self.relationChains[newCtg][orderIndex]

        # get the order index and the name of the new contig
        # for the homologous block
        otherHapOrderIndex = relationToSplit.homologousBlock.orderIndex
        otherHapNewCtg = relationToSplit.homologousBlock.newCtg

        # remove previous relation
        # both from ref2query and from query2ref
        self.relationChains[newCtg].pop(orderIndex)
        self.relationChains[otherHapNewCtg].pop(otherHapOrderIndex)

        # split the homology relation into three parts
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(errorStart, errorEnd)

        # the ref block that has to be contaminated with base errors
        rBlockPart2 = ref2querySplitRelations[1].block
        rBlockPart2.addMisAssemblyBlockList( "Err",
                                             BlockList([(1, rBlockPart2.origEnd - rBlockPart2.origStart + 1)]))
        rBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        rBlockPart2.clearAnnotationStartBlocksForSampling()

        # the query block that is supposed to be collapsed since the
        # ref haplotype is highly erroneous
        qBlockPart2 = ref2querySplitRelations[1].homologousBlock
        qBlockPart2.addMisAssemblyBlockList( "Col",
                                             BlockList([(1, qBlockPart2.origEnd - qBlockPart2.origStart + 1)]))
        qBlockPart2.containsMisAssembly = True
        # blocks with misassembly cannot be used for creating
        # another misassmbley later
        qBlockPart2.clearAnnotationStartBlocksForSampling()

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

    def updateNewCtgAnnotationWeightsForSampling(self, newCtg, parentRelation, childRelations):
        
        assert(newCtg in self.newCtgToIndexForSampling)
        newCtgIndex = self.newCtgToIndexForSampling[newCtg]

        # subtract the old weights related to the parent block
        parentRefBlock = parentRelation.block
        for annot, total in parentRefBlock.annotationStartTotalLengthsForSampling.items():
            # if "updateAnnotationStartLocationsForSampling" is not called
            # then the weight list would be empty
            if 0 < len(self.newCtgAnnotationWeightsForSampling[annot]):
                self.newCtgAnnotationWeightsForSampling[annot][newCtgIndex] -= total

        # add the new weights related to the child blocks
        for relation in childRelations:
            for annot, total in relation.block.annotationStartTotalLengthsForSampling.items():
                # if "updateAnnotationStartLocationsForSampling" is not called
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


    def updateAnnotationStartLocationsForSampling(self, annotations, misAssemblyLength, minOverlapRatioWithEachAnnotation, minMarginLength):
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
                    relation.block.updateAnnotationStartLocationsForSampling()
                else:
                    relation.block.clearAnnotationStartBlocksForSampling()

        for annotation in annotations:
            newCtgWeights = []
            for newCtg in self.newCtgListForSampling:
                newCtgWeights.append(self.getTotalSamplingLength(newCtg, annotation))
            self.newCtgAnnotationWeightsForSampling[annotation] =  newCtgWeights

    def getWeightedRandomNewCtgForSampling(self, annotation):
        """
        Given the annotation name select one new contig randomly by taking the
        total length of sampling regions as the sampling weight for each new contig

        :param annotation: The annotation name
        :return: The randomly selected new contig
        """
        selectedNewCtg = random.choices(self.newCtgList,
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

        Note that the amount of overlap was determined while calling "updateAnnotationStartLocationsForSampling"

        :param annotation: The annotation name
        :param misAssemblyLength: Misassembly length
        :return: a randomly selected new contig name, index of relation, start and end coordinates
        """
        newCtg = self.getWeightedRandomNewCtgForSampling(annotation)
        orderIndex = self.getWeightedRandomOrderIndexForSampling(newCtg, annotation)
        start, end = self.getRandomIntervalFromRelationBlock(newCtg, annotation, orderIndex, misAssemblyLength)
        return newCtg, orderIndex, start, end

    def getNewCtgSequence(self, newCtg, origCtgSequences):
        newCtgSeqList = []
        for relation in self.relationChains[newCtg]:
            blockSeq = relation.block.getSequence(origCtgSequences)
            newCtgSeqList.append(blockSeq)
        return  "".join(newCtgSeqList)

    def yeildNewCtgSequences(self, origCtgSequences):
        for newCtg in self.relationChains:
            yield newCtg, self.getNewCtgSequence(newCtg, origCtgSequences)





