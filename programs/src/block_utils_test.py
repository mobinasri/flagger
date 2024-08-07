from block_utils import *
import unittest
from project_blocks_multi_thread import makeCigarString

class TestProjection(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """Define two alignments one in positive orientation and one in negative"""
        self.alignmentPositive = Alignment("ctg\t350\t100\t150\t+\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignmentNegative = Alignment("ctg\t350\t200\t250\t-\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignmentOverlapping1 = Alignment("ctg2\t10\t0\t10\t-\tctg1\t14\t0\t10\t8\t10\t60\tcg:Z:3=1X2I1=2D1=1X1=\ttp:A:P")
        self.alignmentOverlapping2 = Alignment("ctg3\t5\t0\t5\t+\tctg1\t14\t10\t13\t3\t3\t60\tcg:Z:1=2I1X1=\ttp:A:P")
        self.alignmentOverlapping3 = Alignment("ctg2\t10\t3\t8\t+\tctg1\t14\t4\t12\t5\t8\t60\tcg:Z:4=3D1X\ttp:A:P")
        self.alignmentOverlapping4 = Alignment("ctg2\t10\t5\t10\t+\tctg1\t14\t4\t12\t5\t8\t60\tcg:Z:4=3D1X\ttp:A:P")
        self.alignmentsOverlapping = [self.alignmentOverlapping1, self.alignmentOverlapping2, self.alignmentOverlapping3]
        self.alignmentsOverlapping2 = [self.alignmentOverlapping1, self.alignmentOverlapping2, self.alignmentOverlapping4]
        self.blockList1 = BlockList([(1,10), (14, 20), (21, 30)])
        self.blockList2 = BlockList([(4,6), (8, 18), (25, 40)])
        print(f"Tests:")

    def testIntersectBlockList(self):
        outputIntersect = self.blockList1.intersect(self.blockList2, inplace=False)

        truthIntersect = BlockList([(4, 6, 0), (8, 10, 0), (14, 18, 0), (25, 30, 0)])
        self.assertEqual(len(outputIntersect.blocks) , len(truthIntersect.blocks), "Number of blocks is not correct")

        for i in range(len(truthIntersect.blocks)):
            self.assertEqual(truthIntersect.blocks[i][0], outputIntersect.blocks[i][0], "Start is not correct")
            self.assertEqual(truthIntersect.blocks[i][1], outputIntersect.blocks[i][1], "End is not correct")

    def testSubtractBlockList(self):
        outputSubtract = self.blockList1.subtract(self.blockList2, inplace=False)

        truthSubctract = BlockList([(1, 3, 0), (7, 7, 0), (19, 20, 0), (21, 24, 0)])
        self.assertEqual(len(outputSubtract.blocks) , len(truthSubctract.blocks), "Number of blocks is not correct")

        for i in range(len(truthSubctract.blocks)):
            self.assertEqual(truthSubctract.blocks[i][0], outputSubtract.blocks[i][0], "Start is not correct")
            self.assertEqual(truthSubctract.blocks[i][1], outputSubtract.blocks[i][1], "End is not correct")

    def testPositiveAsm2Ref(self):
        alignment = self.alignmentPositive
        mode = "asm2ref"
        includeEndingIndel = False
        includePostIndel = False
        # block start and end are 1-based
        blocks = [(81, 106, "NA"), (111, 116, "NA"), (124, 126, "NA")]
        projectionBlocks = [(11, 15), (29, 31), (32, 39)]
        projectableBlocks = [(101, 105), (111, 113), (124, 126)]
        cigarListTruth = ["4=1X", "1X2=", "1X1=5D1="]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testPositiveAsm2RefIncludeEndingIndel(self):
        alignment = self.alignmentPositive
        mode = "asm2ref"
        includeEndingIndel = True
        includePostIndel = False
        # block start and end are 1-based
        blocks = [(101, 110, "NA"), (111, 116, "NA"), (117, 126, "NA")]
        projectionBlocks = [(11, 18), (29, 31), (32, 39)]
        projectableBlocks = [(101, 110), (111, 116), (117, 126)]
        cigarListTruth = ["4=1X2I2=1X", "1X2=3I", "7I1X1=5D1="]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testPositiveAsm2RefIncludeEndingAndPostIndel(self):
        alignment = self.alignmentPositive
        mode = "asm2ref"
        includeEndingIndel = True
        includePostIndel = True
        # block start and end are 1-based
        blocks = [(101, 110, "NA"), (111, 116, "NA"), (117, 126, "NA")]
        projectionBlocks = [(11, 28), (29, 31), (32, 39)]
        projectableBlocks = [(101, 110), (111, 116), (117, 126)]
        cigarListTruth = ["4=1X2I2=1X10D", "1X2=3I", "7I1X1=5D1="]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testPositiveRef2Asm(self):
        alignment = self.alignmentPositive
        mode = "ref2asm"
        includeEndingIndel = False
        includePostIndel = False
        # block start and end are 1-based
        blocks = [(11, 16, "NA"), (29, 31, "NA"), (32, 35, "NA")]
        projectionBlocks = [(101, 108), (111, 113), (124,125)]
        projectableBlocks = [(11, 16), (29, 31), (32, 33)]
        cigarListTruth = ["4=1X2I1=", "1X2=", "1X1="]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testPositiveRef2AsmIncludeEndingAndPostIndel(self):
        alignment = self.alignmentPositive
        mode = "ref2asm"
        includeEndingIndel = True
        includePostIndel = True
        # block start and end are 1-based
        blocks = [(11, 20, "NA"), (21, 31, "NA"), (32, 35, "NA")]
        projectionBlocks = [(101, 110), (111, 123), (124,125)]
        projectableBlocks = [(11, 20), (21, 31), (32, 35)]
        cigarListTruth = ["4=1X2I2=1X2D", "8D1X2=10I", "1X1=2D"]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testPositiveRef2AsmIncludeEndingIndel(self):
        alignment = self.alignmentPositive
        mode = "ref2asm"
        includeEndingIndel = True
        includePostIndel = True
        # block start and end are 1-based
        blocks = [(11, 20, "NA"), (21, 31, "NA"), (32, 48, "NA")]
        projectionBlocks = [(101, 110), (111, 123), (124,140)]
        projectableBlocks = [(11, 20), (21, 31), (32, 48)]
        cigarListTruth = ["4=1X2I2=1X2D", "8D1X2=10I", "1X1=5D5=1X4=5I"]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            self.assertEqual(qBlocks[i][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[i][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[i][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[i][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[i], makeCigarString(cigarList[i]), "Incorrect CIGAR string")


    def testNegativeAsm2Ref(self):
        alignment = self.alignmentNegative
        mode = "asm2ref"
        includeEndingIndel = False
        includePostIndel = False
        # block start and end are 1-based
        blocks = [(150, 211, "NA"), (212, 225, "NA"), (240, 241, "NA")]
        projectionBlocks = [(49, 58), (39, 48), (18, 29)]
        projectableBlocks = [(201, 210), (216, 225), (240, 241)]
        cigarListTruth = ["3=1X6=", "5=1X4=", "1X10D1X"]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            j = -1-i
            self.assertEqual(qBlocks[j][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[j][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[j][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[j][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[j], makeCigarString(cigarList[i]), "Incorrect CIGAR string")


    def testNegativeAsm2RefIncludeEndingIndel(self):
        alignment = self.alignmentNegative
        mode = "asm2ref"
        includeEndingIndel = True
        includePostIndel = False
        # block start and end are 1-based
        blocks = [(150, 211, "NA"), (212, 225, "NA"), (240, 241, "NA")]
        projectionBlocks = [(49, 58), (39, 48), (18, 29)]
        projectableBlocks = [(201, 211), (212, 225), (240, 241)]
        cigarListTruth = ["1I3=1X6=", "5=1X4=4I", "1X10D1X"]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            j = -1-i
            self.assertEqual(qBlocks[j][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[j][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[j][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[j][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[j], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testNegativeAsm2RefIncludeEndingAndPostIndel(self):
        alignment = self.alignmentNegative
        mode = "asm2ref"
        includeEndingIndel = True
        includePostIndel = True
        # block start and end are 1-based
        blocks = [(150, 211, "NA"), (226, 227, "NA"), (240, 241, "NA")]
        projectionBlocks = [(49, 58), (32, 38), (18, 29)]
        projectableBlocks = [(201, 211), (226, 227), (240, 241)]
        cigarListTruth = ["1I3=1X6=", "1X1=5D", "1X10D1X"]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(3):
            j = -1-i
            self.assertEqual(qBlocks[j][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[j][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[j][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[j][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[j], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testNegativeRef2AsmIncludeEndingAndPostIndel(self):
        alignment = self.alignmentNegative
        mode = "ref2asm"
        includeEndingIndel = True
        includePostIndel = True
        # block start and end are 1-based
        blocks = [(14, 20, "NA"), (39, 48, "NA"), (49, 58, "NA")]
        projectionBlocks = [(241, 247), (211, 225), (201, 210)]
        projectableBlocks = [(14, 20), (39, 48) , (49, 58)]
        cigarListTruth = ["1=1X2I2=1X2D" , "5=1X4=5I", "3=1X6="]
        qBlocks, rBlocks, cigarList = findProjections(mode,
                                                      alignment.cigarList,
                                                      blocks,
                                                      alignment.chromLength,
                                                      alignment.chromStart + 1, alignment.chromEnd, # make 1-based start
                                                      alignment.contigLength,
                                                      alignment.contigStart + 1, alignment.contigEnd, # make 1-based start
                                                      alignment.orientation,
                                                      includeEndingIndel, includePostIndel)
        for i in range(2):
            j = -1-i
            self.assertEqual(qBlocks[j][0], projectableBlocks[i][0], "Incorrect projectable start position")
            self.assertEqual(qBlocks[j][1], projectableBlocks[i][1], "Incorrect projectable end position")
            self.assertEqual(rBlocks[j][0], projectionBlocks[i][0], "Incorrect projection start position")
            self.assertEqual(rBlocks[j][1], projectionBlocks[i][1], "Incorrect projection end position")
            self.assertEqual(cigarListTruth[j], makeCigarString(cigarList[i]), "Incorrect CIGAR string")

    def testCountAlignmentOverlapsRef(self):
        alignments = self.alignmentsOverlapping
        blockListsPerRefContig = getBlockListsPerRefContig(alignments)
        mergedBlockListsOutput = mergeBlockListsPerContigWithOverlapCount(blockListsPerRefContig)

        mergedBlockListsTruth = {"ctg1": BlockList([ (1, 4, 1), (5, 10, 2), (11, 12, 2), (13, 13, 1)]) }
 
        self.assertListEqual(list(mergedBlockListsTruth.keys()), list(mergedBlockListsOutput.keys()), "Incorrect contig names")
        self.assertEqual(len(mergedBlockListsTruth["ctg1"].blocks), len(mergedBlockListsOutput["ctg1"].blocks), "Incorrect blocks length")

        for i in range(len(mergedBlockListsTruth["ctg1"].blocks)):
            truthBlock = mergedBlockListsTruth["ctg1"].blocks[i]
            outputBlock = mergedBlockListsOutput["ctg1"].blocks[i]
            self.assertEqual(truthBlock[0], outputBlock[0], "Incorrect start position")
            self.assertEqual(truthBlock[1], outputBlock[1], "Incorrect end position")
            

    #@unittest.SkipTest
    def testRefUniqueAlignments(self):
        alignments = self.alignmentsOverlapping
        uniqueBlockListsPerRefContig = getBlockListsWithSingleAlignmentPerRefContig(alignments)
        alignmentsOutput = subsetAlignmentsToRefBlocks(alignments, uniqueBlockListsPerRefContig)

        refUniqueAlignment1 = Alignment("ctg2\t10\t6\t10\t-\tctg1\t14\t0\t4\t4\t4\t60\tcg:Z:3=1X\ttp:A:P")
        refUniqueAlignment2 = Alignment("ctg3\t5\t4\t5\t+\tctg1\t14\t12\t13\t1\t1\t60\tcg:Z:1=\ttp:A:P")
        refUniqueAlignments = [refUniqueAlignment1, refUniqueAlignment2]

        #print(alignmentsOutput[0].cigarList, alignmentsOutput[0].contigStart, alignmentsOutput[0].contigEnd)
        #print(alignmentsOutput[1].cigarList, alignmentsOutput[1].contigStart, alignmentsOutput[1].contigEnd)

        #print(alignmentsOutput[2].cigarList)

        self.assertEqual(len(refUniqueAlignments), len(alignmentsOutput), "Incorrect length of total alignments")
        for i in range(len(alignmentsOutput)):
            a = alignmentsOutput[i]
            t = refUniqueAlignments[i]
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

    def testQueryUniqueAlignments(self):
        alignments = self.alignmentsOverlapping
        uniqueBlockListsPerQueryContig = getBlockListsWithSingleAlignmentPerQueryContig(alignments)
        alignmentsOutput = subsetAlignmentsToQueryBlocks(alignments, uniqueBlockListsPerQueryContig)

        queryUniqueAlignment1 = Alignment("ctg2\t10\t8\t10\t-\tctg1\t14\t0\t2\t2\t2\t60\tcg:Z:2=\ttp:A:P")
        queryUniqueAlignment2 = Alignment("ctg2\t10\t0\t3\t-\tctg1\t14\t7\t10\t1\t1\t60\tcg:Z:1=1X1=\ttp:A:P")
        queryUniqueAlignment3 = Alignment("ctg3\t5\t0\t5\t+\tctg1\t14\t10\t13\t3\t3\t60\tcg:Z:1=2I1X1=\ttp:A:P")
        queryUniqueAlignments = [queryUniqueAlignment1, queryUniqueAlignment2, queryUniqueAlignment3]


        self.assertEqual(len(queryUniqueAlignments), len(alignmentsOutput), "Incorrect length of total alignments")
        for i in range(len(alignmentsOutput)):
            a = alignmentsOutput[i]
            t = queryUniqueAlignments[i]
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")


    def testTwoWayUniqueAlignments(self):
        alignments = self.alignmentsOverlapping
        uniqueBlockListsPerQueryContig = getBlockListsWithSingleAlignmentPerQueryContig(alignments)
        uniqueBlockListsPerRefContig = getBlockListsWithSingleAlignmentPerRefContig(alignments)
        refUniqueAlignments = subsetAlignmentsToRefBlocks(alignments, uniqueBlockListsPerRefContig)
        alignmentsOutput = subsetAlignmentsToQueryBlocks(refUniqueAlignments, uniqueBlockListsPerQueryContig)



        twoWayUniqueAlignment1 = Alignment("ctg2\t10\t8\t10\t-\tctg1\t14\t0\t2\t2\t2\t60\tcg:Z:2=\ttp:A:P")
        twoWayUniqueAlignment2 = Alignment("ctg3\t5\t4\t5\t+\tctg1\t14\t12\t13\t1\t1\t60\tcg:Z:1=\ttp:A:P")
        twoWayUniqueAlignments = [twoWayUniqueAlignment1, twoWayUniqueAlignment2]


        self.assertEqual(len(twoWayUniqueAlignments), len(alignmentsOutput), "Incorrect length of total alignments")
        for i in range(len(alignmentsOutput)):
            a = alignmentsOutput[i]
            t = twoWayUniqueAlignments[i]
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

    def testTwoWayUniqueAlignments_2(self):
        alignments = self.alignmentsOverlapping2
        uniqueBlockListsPerQueryContig = getBlockListsWithSingleAlignmentPerQueryContig(alignments)
        uniqueBlockListsPerRefContig = getBlockListsWithSingleAlignmentPerRefContig(alignments)
        refUniqueAlignments = subsetAlignmentsToRefBlocks(alignments, uniqueBlockListsPerRefContig)
        alignmentsOutput = subsetAlignmentsToQueryBlocks(refUniqueAlignments, uniqueBlockListsPerQueryContig)

        twoWayUniqueAlignment1 = Alignment("ctg3\t5\t4\t5\t+\tctg1\t14\t12\t13\t1\t1\t60\tcg:Z:1=\ttp:A:P")
        twoWayUniqueAlignments = [twoWayUniqueAlignment1]

        self.assertEqual(len(twoWayUniqueAlignments), len(alignmentsOutput), "Incorrect length of total alignments")
        for i in range(len(alignmentsOutput)):
            a = alignmentsOutput[i]
            t = twoWayUniqueAlignments[i]
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")
    
    def testBlockListPerContigSplit_1(self):
        blockListPerContig = {"contig1": BlockList(),
                              "contig2": BlockList([(1,10), (41, 55)]),
                              "contig3": BlockList([(16,50)])}
        outParts = BlockList.split(blockListPerContig=blockListPerContig,
                                   numberOfParts=4)
        part1 = {"contig2": BlockList([(1, 10), (41,45)])}
        part2 = {"contig2": BlockList([(46, 55)]),
                 "contig3": BlockList([(16, 20)])}
        part3 = {"contig3": BlockList([(21, 35)])}
        part4 = {"contig3": BlockList([(36, 50)])}

        truthParts = [part1, part2, part3, part4]

        self.assertEqual(len(truthParts), len(outParts))
        for truthPart, outPart in zip(truthParts, outParts):
            self.assertCountEqual(list(truthPart.keys()), list(outPart.keys()))
            for contig, truthBlockList in truthPart.items():
                self.assertTrue(truthBlockList.isEqual(outPart[contig]))

    def testBlockListPerContigSplit_2(self):
        blockListPerContig = {"contig1": BlockList(),
                              "contig2": BlockList([(1,10), (41, 55)]),
                              "contig3": BlockList([(16,51)]),
                              "contig4": BlockList()}
        outParts = BlockList.split(blockListPerContig=blockListPerContig,
                                   numberOfParts=4)
        part1 = {"contig2": BlockList([(1, 10), (41,46)])}
        part2 = {"contig2": BlockList([(47, 55)]),
                 "contig3": BlockList([(16, 22)])}
        part3 = {"contig3": BlockList([(23, 38)])}
        part4 = {"contig3": BlockList([(39, 51)])}

        truthParts = [part1, part2, part3, part4]

        self.assertEqual(len(truthParts), len(outParts))
        for truthPart, outPart in zip(truthParts, outParts):
            self.assertCountEqual(list(truthPart.keys()), list(outPart.keys()))
            for contig, truthBlockList in truthPart.items():
                self.assertTrue(truthBlockList.isEqual(outPart[contig]))

def main():
    unittest.main(verbosity=2)

if __name__ == '__main__':
    main()

