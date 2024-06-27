from block_utils import *
from alignment_utils import *
import unittest
from project_blocks_multi_thread import makeCigarString

class TestProjection(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """Define two alignments one in positive orientation and one in negative"""
        self.alignmentPositive = Alignment("ctg\t350\t100\t150\t+\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignmentNegative = Alignment("ctg\t350\t200\t250\t-\tref\t158\t10\t58\t27\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        print(f"Tests:")

    
    def testSplitIntoAlignmentsWithShortGaps_Positive_shortGap(self):
        alignmentPositiveSplit1 = Alignment("ctg\t350\t100\t110\t+\tref\t158\t10\t18\t0\t0\t0\tcg:Z:4=1X2I2=1X\ttp:A:P")
        alignmentPositiveSplit2 = Alignment("ctg\t350\t110\t113\t+\tref\t158\t28\t31\t0\t0\t0\tcg:Z:1X2=\ttp:A:P")
        alignmentPositiveSplit3 = Alignment("ctg\t350\t123\t150\t+\tref\t158\t31\t58\t0\t0\t0\tcg:Z:1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        splitAlignmentsTruth = [alignmentPositiveSplit1, alignmentPositiveSplit2, alignmentPositiveSplit3]

        splitAlignmentsActual = splitIntoAlignmentsWithShortGaps(self.alignmentPositive, 5)

        for a, t in zip(splitAlignmentsActual, splitAlignmentsTruth):
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

    def testSplitIntoAlignmentsWithShortGaps_Positive_longGap(self):
        splitAlignmentsTruth = [self.alignmentPositive]

        splitAlignmentsActual = splitIntoAlignmentsWithShortGaps(self.alignmentPositive, 100)

        for a, t in zip(splitAlignmentsActual, splitAlignmentsTruth):
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

    def testSplitIntoAlignmentsWithShortGaps_Negative_shortGap(self):
        alignmentPositiveSplit1 = Alignment("ctg\t350\t240\t250\t+\tref\t158\t10\t18\t0\t0\t0\tcg:Z:4=1X2I2=1X\ttp:A:P")
        alignmentPositiveSplit2 = Alignment("ctg\t350\t237\t240\t+\tref\t158\t28\t31\t0\t0\t0\tcg:Z:1X2=\ttp:A:P")
        alignmentPositiveSplit3 = Alignment("ctg\t350\t200\t227\t+\tref\t158\t31\t58\t0\t0\t0\tcg:Z:1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        splitAlignmentsTruth = [alignmentPositiveSplit1, alignmentPositiveSplit2, alignmentPositiveSplit3]

        splitAlignmentsActual = splitIntoAlignmentsWithShortGaps(self.alignmentNegative, 5)

        for a, t in zip(splitAlignmentsActual, splitAlignmentsTruth):
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

    def testSplitIntoAlignmentsWithShortGaps_Negative_longGap(self):
        splitAlignmentsTruth = [self.alignmentNegative]

        splitAlignmentsActual = splitIntoAlignmentsWithShortGaps(self.alignmentNegative, 100)

        for a, t in zip(splitAlignmentsActual, splitAlignmentsTruth):
            self.assertEqual(t.chromLength, a.chromLength, "Incorrect chrom length")
            self.assertEqual(t.chromStart, a.chromStart, "Incorrect chrom start")
            self.assertEqual(t.chromEnd, a.chromEnd, "Incorrect chrom end")
            self.assertEqual(t.contigLength, a.contigLength, "Incorrect contig length")
            self.assertEqual(t.contigStart, a.contigStart, "Incorrect contig start")
            self.assertEqual(t.contigEnd, a.contigEnd, "Incorrect contig end")
            self.assertEqual(t.orientation, a.orientation, "Incorrect orientation")
            self.assertListEqual(t.cigarList, a.cigarList, "Incorrect CIGAR")

def main():
    unittest.main(verbosity=2)

if __name__ == '__main__':
    main()

