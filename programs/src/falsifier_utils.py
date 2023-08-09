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
    def createRef2QueryRelationFromAlignment(alignment: Alignment, newCtgSuffix: str) -> HomologyRelation:
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

    # given start should be 0-based
    # given end should be 1-based
    @staticmethod
    def createVoidRelationFromInterval(ctgName: str, ctgStart: int, ctgEnd: int, newCtgSuffix: str) -> HomologyRelation:
        block = HomologyBlock(ctgName,
                              ctgStart + 1,
                              ctgEnd,
                              '+',
                              f'{ctgName}{newCtgSuffix}',
                              None)
        voidRelation = HomologyRelation(block, None, None)
        return  voidRelation

    @staticmethod
    def createAllInclusiveRelationDictFromAlignments(alignments: list, contigLengths: dict, newCtgSuffix: str) -> list:
        relationFreeIntervals = []
        for ctgName, ctgLen in contigLengths.items():
            relationFreeIntervals[ctgName] = [(0, ctgLen)]
        relationDict = defaultdict(list)

        # create two-way relations; ref2query and query2ref and add them to the dictionary
        for alignment in alignments:
            # the homology blocks (non-void ones) are shared between ref2queryRelation and query2refRelation
            # to be able to have access to the updated index of the homologous block
            ref2queryRelation = HomologyRelation.createRef2QueryRelationFromAlignment(alignment, newCtgSuffix)
            query2refRelation = HomologyRelation(ref2queryRelation.homologousBlock,
                                                 ref2queryRelation.block,
                                                 None)
            relationDict[alignment.chromName].append(ref2queryRelation)
            relationDict[alignment.contigName].append(query2refRelation)

            # update relation-free intervals
            # to create void relations from them
            # after the alignments are all iterated
            subtractInterval(relationFreeIntervals[alignment.chromName], (alignment.chromStart, alignment.chromEnd))
            subtractInterval(relationFreeIntervals[alignment.contigName], (alignment.contigStart, alignment.contigEnd))

        # create void relations for the intervals without alignments
        for ctgName, intervals in relationFreeIntervals.items():
            for interval in intervals:
                voidRelation = HomologyRelation.createVoidRelationFromInterval(ctgName, interval[0], interval[1], newCtgSuffix)
                relationDict[ctgName].append(voidRelation)

        # sort the relations for each contig based on start coordinates
        for ctgName in contigLengths:
            relationDict[ctgName].sort(key = lambda x : x.block.origStart)
        return relationDict



