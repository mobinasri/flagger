from block_utils import *
from collections import defaultdict
from copy import deepcopy

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
        self.annotationBlockListsToBeSampled = defaultdict(BlockList)
        self.misAssemblyBlockLists = defaultdict(BlockList)

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

    def updateOneAnnotationBlockListToBeSampled(self, name, lengthToTruncateFromEnd, wholeBlockMargin):
        """

        :param name: The name of the annotation to update
        :param lengthToTruncateFromEnd: The length of each annotation block to
                                        be excluded from sampling (from the right side)
        :param wholeBlockMargin: The margin of the whole block to be excluded from sampling
                                 (from both side of the whole block)
        """
        annotationBlockListToBeSampled = self.annotationBlockLists[name].copy()
        # truncate the right side of the blocks
        annotationBlockListToBeSampled.truncateFromEnd(lengthToTruncateFromEnd, inplace=True)

        # truncate the blocks to make sure they are far enough
        # from the edges of the whole block
        wholeBlockWithoutMargin = BlockList([(1, self.origEnd - self.origStart + 1)]).truncateFromBothSides(wholeBlockMargin, inplace=False)
        annotationBlockListToBeSampled.intersect(wholeBlockWithoutMargin, inplace=True)
        self.annotationBlockListsToBeSampled[name] = annotationBlockListToBeSampled

    def updateAllAnnotationBlockListsToBeSampled(self, lengthToTruncateFromEnd, wholeBlockMargin):
        """
        Run "updateOneAnnotationBlockListToBeSampled" for all existing annotations
        """
        for name in self.annotationBlockLists:
            self.updateOneAnnotationBlockListToBeSampled(name, lengthToTruncateFromEnd, wholeBlockMargin)
    
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

        rBlockPart2 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[1][0],
                                    projectionsOrigCoor[1][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 1)
        rBlockPart2.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])

        rBlockPart3 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[2][0],
                                    projectionsOrigCoor[2][1],
                                    '+',
                                    rBlock.newCtg,
                                    rBlock.orderIndex + 2)
        rBlockPart3.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])

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

        qBlockPart2 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[1][2],
                                    projectionsOrigCoor[1][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart2)
        qBlockPart2.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])

        qBlockPart3 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[2][2],
                                    projectionsOrigCoor[2][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart3)
        qBlockPart3.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[2][2], projectionsRelCoor[2][3])

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
    def __init__(self, alignments, origContigLengths, rightAnnotationMarginLength, blockMarginLength, newCtgSuffix):
        """
            A class for saving chains of homology relations for each contig

            :param rightAnnotationMarginLength: The length of the margin from the right side of each continuous annotation
            block in the reference/hap1 block of each homology relation
            These margins are NOT ALLOWED to be the start coordinates of misassemblies and they must be set based on
            how many bases of each misassembly is desired to have overlap with each annotation

            :param blockMarginLength: The length of the margins from the both sides of the reference/hap1 block of
            each homology relation. It MUST BE GREATER THAN the length of the desired misassembly
            The main purpose of this these margins is to make sure that misassemblies do not exceed the interval of each
            homology relation and the whole altered segment is happening within the corresponding relation.
            Additionally, the misassemblies will not be created near to the edges of alignments where they could exist
            unreliable homology relation.

            Note than blockMarginLength and rightAnnotationMarginLength should be modified once the length of desired
            misassemblies is changed
            For other params read the documentation for "createAllInclusiveRelationChainsFromAlignments"
        """
        self.relationChains = HomologyRelationChains.createAllInclusiveRelationChainsFromAlignments(alignments,
                                                                                                    origContigLengths,
                                                                                                    newCtgSuffix)
        self.rightAnnotationMarginLength = rightAnnotationMarginLength
        self.blockMarginLength = blockMarginLength

    def updateMargins(self, rightAnnotationMarginLength, blockMarginLength):
        self.rightAnnotationMarginLength = rightAnnotationMarginLength
        self.blockMarginLength = blockMarginLength

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

    def induceSwitchMisAssembly(self, newCtg, orderIndex, switchStart, switchEnd):

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
        ref2querySplitRelations = relationToSplit.splitIntoThreeParts(switchStart, switchEnd)

        # the blocks that have to be swapped
        relationPart2 = ref2querySplitRelations[1]
        rBlockPart2 = relationPart2.block
        qBlockPart2 = relationPart2.homologousBlock

        # swap blocks, order indices and new contig names
        rBlockPart2.orderIndex, qBlockPart2.orderIndex = qBlockPart2.orderIndex,  rBlockPart2.orderIndex
        rBlockPart2.newCtg, qBlockPart2.newCtg = qBlockPart2.newCtg,  rBlockPart2.newCtg
        relationPart2.block, relationPart2.homologousBlock = relationPart2.homologousBlock, relationPart2.block

        # convert cigar if the alignment orientation is negative
        if relationPart2.alignment.orientation == '-':
            relationPart2.alignment.cigarList = convertIndelsInCigar(relationPart2.alignment.cigarList)
            rBlockPart2.origStrand = '-'
            qBlockPart2.origStrand = '-'

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


        # falsely duplicated block, this is the duplication of the middle part of the rBlock
        newDupCtg = rBlockPart2.origCtg + f"_Dup_{rBlockPart2.origStart}_{rBlockPart2.origEnd}"
        rBlockPart2Dup = HomologyBlock(rBlockPart2.origCtg,
                                       rBlockPart2.origStart,
                                       rBlockPart2.origEnd,
                                       '+',
                                       newDupCtg,
                                       0)
        rBlockPart2Dup.annotationBlockLists = deepcopy(rBlockPart2.annotationBlockLists)
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


