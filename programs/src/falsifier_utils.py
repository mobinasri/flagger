from block_utils import *
from collections import defaultdict


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

    # switchStart and switchEnd must be 1-based and closed
    def induceSwitchError(self, switchStart, switchEnd):

        # This function only works when the strand of the whole block is positive
        # In other words it is not possible to create a switch within a previously
        # created haplotype switch
        assert(self.block.origStrand == '+')
        assert(self.homologousBlock.origStrand == '+')
        forwardBlocks = [(1, switchStart - 1, ""), (switchStart, switchEnd, ""), (switchEnd + 1, self.alignment.chromLength, "")]
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

        rBlockPart2 = HomologyBlock(qBlock.origCtg,
                                       projectionsOrigCoor[1][2],
                                       projectionsOrigCoor[1][3],
                                       self.alignment.orientation,
                                       rBlock.newCtg,
                                       rBlock.orderIndex + 1)
        rBlockPart2.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[1][2], projectionsRelCoor[1][3])

        rBlockPart3 = HomologyBlock(rBlock.origCtg,
                                      projectionsOrigCoor[2][0],
                                      projectionsOrigCoor[2][1],
                                      '+',
                                      rBlock.newCtg,
                                      rBlock.orderIndex + 2)
        rBlockPart3.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[2][0], projectionsRelCoor[2][1])

        # query blocks
        qOrderIndexPart1 = qBlock.orderIndex if self.alignment.orientation  == '+' else qBlock.orderIndex + 2
        qOrderIndexPart2 = qBlock.orderIndex + 1
        qOrderIndexPart3 = qBlock.orderIndex + 2 if self.alignment.orientation  == '+' else qBlock.orderIndex

        qBlockPart1 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[0][2],
                                    projectionsOrigCoor[0][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart1)
        qBlockPart1.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])

        qBlockPart2 = HomologyBlock(rBlock.origCtg,
                                    projectionsOrigCoor[1][0],
                                    projectionsOrigCoor[1][1],
                                    self.alignment.orientation,
                                    qBlock.newCtg,
                                    qOrderIndexPart2)
        qBlockPart2.extractAnnotationsFromParentBlock(rBlock, projectionsRelCoor[1][0], projectionsRelCoor[1][1])

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
                                                 projectionsOrigCoor[1][4] if self.alignment.orientation == '+' else convertIndelsInCigar(projectionsOrigCoor[1][4]),
                                                 self.alignment.orientation)
        relationPart3 = HomologyRelation(rBlockPart3,
                                                 qBlockPart3,
                                                 projectionsOrigCoor[2][4],
                                                 self.alignment.orientation)

        homologyRelations = [relationPart1, relationPart2, relationPart3]

        return  homologyRelations

    def induceCollapseError(self, collapseStart, collapseEnd):

        # This function only works when the strand of the whole block is positive
        # In other words it is not possible to create a switch within a previously
        # created haplotype switch
        assert(self.block.origStrand == '+')
        assert(self.homologousBlock.origStrand == '+')
        forwardBlocks = [(1, collapseStart - 1, ""), (collapseStart, collapseEnd, ""), (collapseEnd + 1, self.alignment.chromLength, "")]
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
        qOrderIndexPart1 = qBlock.orderIndex if self.alignment.orientation  == '+' else qBlock.orderIndex + 1
        qOrderIndexPart3 = qBlock.orderIndex + 1 if self.alignment.orientation  == '+' else qBlock.orderIndex

        qBlockPart1 = HomologyBlock(qBlock.origCtg,
                                    projectionsOrigCoor[0][2],
                                    projectionsOrigCoor[0][3],
                                    '+',
                                    qBlock.newCtg,
                                    qOrderIndexPart1)
        qBlockPart1.extractAnnotationsFromParentBlock(qBlock, projectionsRelCoor[0][2], projectionsRelCoor[0][3])

        # skip qBlockPart2 since we want to make rBlockPart1 completely collapsed

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
                                         None,
                                         None,
                                         None)
        relationPart3 = HomologyRelation(rBlockPart3,
                                         qBlockPart3,
                                         projectionsOrigCoor[2][4],
                                         self.alignment.orientation)

        homologyRelations = [relationPart1, relationPart2, relationPart3]

        return  homologyRelations

    @staticmethod
    def createRef2QueryRelationFromAlignment(alignment: Alignment, newCtgSuffix: str):
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

    # given start should be 1-based
    # given end should be 1-based
    @staticmethod
    def createVoidRelationFromInterval(ctgName: str, ctgStart: int, ctgEnd: int, newCtgSuffix: str):
        block = HomologyBlock(ctgName,
                              ctgStart,
                              ctgEnd,
                              '+',
                              f'{ctgName}{newCtgSuffix}',
                              None)
        voidRelation = HomologyRelation(block, None, None, None)
        return  voidRelation

    @staticmethod
    def createAllInclusiveRelationDictFromAlignments(alignments: list, contigLengths: dict, newCtgSuffix: str) -> list:
        relationFreeIntervals = {}
        for ctgName, ctgLen in contigLengths.items():
            relationFreeIntervals[ctgName + newCtgSuffix] = BlockList([(1, ctgLen)])
        relationDict = defaultdict(list)

        # create two-way relations; ref2query and query2ref and add them to the dictionary
        for alignment in alignments:
            # the homology blocks (non-void ones) are shared between ref2queryRelation and query2refRelation
            # to be able to have access to the updated index of the homologous block
            ref2queryRelation = HomologyRelation.createRef2QueryRelationFromAlignment(alignment, newCtgSuffix)
            query2refRelation = HomologyRelation(ref2queryRelation.homologousBlock,
                                                 ref2queryRelation.block,
                                                 None,
                                                 None)
            relationDict[ref2queryRelation.block.newCtg].append(ref2queryRelation)
            relationDict[query2refRelation.block.newCtg].append(query2refRelation)

            # update relation-free intervals
            # to create void relations from them
            # after the alignments are all iterated
            relationFreeIntervals[ref2queryRelation.block.newCtg].subtract(BlockList([(alignment.chromStart + 1, alignment.chromEnd)]), inplace=True)
            relationFreeIntervals[query2refRelation.block.newCtg].subtract(BlockList([(alignment.contigStart + 1, alignment.contigEnd)]), inplace=True)

        # create void relations for the intervals without alignments
        for ctgName in contigLengths:
            for interval in relationFreeIntervals[ctgName + newCtgSuffix].blocks:
                voidRelation = HomologyRelation.createVoidRelationFromInterval(ctgName, interval[0], interval[1], newCtgSuffix)
                relationDict[voidRelation.block.newCtg].append(voidRelation)

        # sort the relations for each contig based on start coordinates
        for ctgName in relationDict:
            relationDict[ctgName].sort(key = lambda x : x.block.origStart)
            for i, relation in enumerate(relationDict[ctgName]):
                relation.block.orderIndex = i
        return relationDict

    @staticmethod
    def fillAnnotationBlockListsFromOriginalContigs(relationsDict, annotationBlockListsPerOrigContig, contigLengths, newCtgSuffix):
        for ctgName, annotationBlockLists in annotationBlockListsPerOrigContig.items():
            newCtgName = ctgName + newCtgSuffix
            # create a homology block for the whole original contig,
            # including the annotations
            wholeOrigContigBlock = HomologyBlock(ctgName, 1, contigLengths[ctgName], '+', newCtgName, 0)
            for name, blockList in annotationBlockLists.items():
                wholeOrigContigBlock.addAnnotationBlockList(name, blockList)
            # the created homology block will then be used for extracting the
            # annotations related to each relation.block
            for relation in relationsDict[newCtgName]:
                relation.block.extractAnnotationsFromParentBlock(wholeOrigContigBlock,
                                                                 relation.block.origStart,
                                                                 relation.block.origEnd)

    @staticmethod
    def induceSwitchErrorAndUpdateRelationsInNewContig(relationsDict, newCtg, orderIndex, switchStart, switchEnd):
        relationsOneCtg = relationsDict[newCtg]
        switchingRelation = relationsOneCtg[orderIndex]
        newCtgOtherHap = switchingRelation.homologousBlock.newCtg
        orderIndexOtherHap = switchingRelation.homologousBlock.orderIndex

        # remove previous relation
        # both from ref2query and from query2ref
        relationsDict[newCtg].pop(orderIndex)
        relationsDict[newCtgOtherHap].pop(orderIndexOtherHap)

        # split relation into three parts and switch the middle part
        ref2queryRelations = switchingRelation.induceSwitchError(switchStart, switchEnd)
        query2refRelations = []

        if switchingRelation.alignment.orientation == '+':
            for relation in ref2queryRelations:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refRelations.append(query2refRelation)
        else:
            for relation in ref2queryRelations[::-1]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refRelations.append(query2refRelation)

        for relation in ref2queryRelations:
            relationsDict[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in relationsDict[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        for relation in query2refRelations:
            relationsDict[newCtgOtherHap].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in relationsDict[newCtgOtherHap][orderIndex + 3:]:
            relation.block.orderIndex += 2

    @staticmethod
    def induceCollapseErrorAndUpdateRelationsInNewContig(relationsDict, newCtg, orderIndex, collapseStart, collapseEnd):
        relationsOneCtg = relationsDict[newCtg]
        collapsingRelation = relationsOneCtg[orderIndex]
        newCtgOtherHap = collapsingRelation.homologousBlock.newCtg
        orderIndexOtherHap = collapsingRelation.homologousBlock.orderIndex

        # remove previous relation
        # both from ref2query and from query2ref
        relationsDict[newCtg].pop(orderIndex)
        relationsDict[newCtgOtherHap].pop(orderIndexOtherHap)

        # split relation into three parts and switch the middle part
        ref2queryRelations = collapsingRelation.induceCollapseError(collapseStart, collapseEnd)
        query2refRelations = []

        if collapsingRelation.alignment.orientation == '+':
            for relation in [ref2queryRelations[0], ref2queryRelations[2]]:
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refRelations.append(query2refRelation)
        else:
            for relation in [ref2queryRelations[2], ref2queryRelations[0]]: # skip the middle relation in collapsing
                query2refRelation = HomologyRelation(relation.homologousBlock,
                                                     relation.block,
                                                     None,
                                                     None)
                query2refRelations.append(query2refRelation)

        for relation in ref2queryRelations:
            relationsDict[newCtg].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in relationsDict[newCtg][orderIndex + 3:]:
            relation.block.orderIndex += 2

        for relation in query2refRelations:
            relationsDict[newCtgOtherHap].insert(relation.block.orderIndex, relation)

        # shift the indices of all the blocks after the last added relation by two
        for relation in relationsDict[newCtgOtherHap][orderIndex + 2:]:
            relation.block.orderIndex += 1




