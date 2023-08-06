from block_utils import *
import unittest
from project_blocks_multi_thread import makeCigarString

class TestProjection(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """Define two alignments one in positive orientation and one in negative"""
        self.alignmentPositive = Alignment("ctg\t350\t100\t150\t+\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignmentNegative = Alignment("ctg\t350\t200\t250\t-\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignmentOverlapping1 = Alignment("ctg1\t14\t0\t10\t-\tctg2\t10\t0\t10\t8\t10\t60\tcg:Z:3=1X2I1=2D1=1X1=\ttp:A:P")
        self.alignmentOverlapping2 = Alignment("ctg1\t14\t10\t13\t+\tctg3\t5\t0\t5\t3\t3\t60\tcg:Z:1=2I1X1=\ttp:A:P")
        self.alignmentOverlapping3 = Alignment("ctg1\t14\t4\t12\t+\tctg4\t5\t0\t5\t3\t5\t60\tcg:Z:4=3D1X\ttp:A:P")
        self.alignmentsOverlapping = [self.alignmentOverlapping1, self.alignmentOverlapping2, self.alignmentOverlapping3]
        print(f"Tests:")

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

<<<<<<< HEAD

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

    def testRefUniqueAlignments(self):
        alignments = self.alignmentsOverlapping
        refUniqueBlocksPerContig = getBlocksWithSingleAlignmentPerRefContig(alignments)
        alignmentsOutput = subsetAlignmentsToRefBlocks(alignments, refUniqueBlocksPerContig)

        refUniqueAlignment1 = Alignment("ctg1\t14\t0\t4\t-\tctg2\t10\t6\t10\t4\t4\t60\tcg:Z:3=1X\ttp:A:P")
        refUniqueAlignment2 = Alignment("ctg1\t14\t12\t13\t+\tctg3\t5\t4\t5\t1\t1\t60\tcg:Z:1=\ttp:A:P")
        refUniqueAlignments = [refUniqueAlignment1, refUniqueAlignment2]

        self.assertCountEqual(refUniqueAlignments, alignmentsOutput, "Incorrect length of total alignments")
        for i in range(len(alignmentsOutput)):
            a = alignmentsOutput[i]
            t = refUniqueAlignments[i]
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

def main():
    unittest.main(verbosity=2)

if __name__ == '__main__':
    main()

