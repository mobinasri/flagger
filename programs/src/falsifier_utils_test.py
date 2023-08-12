from block_utils import *
from falsifier_utils import *
import unittest
from project_blocks_multi_thread import makeCigarString

class TestProjection(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """Define two alignments one in positive orientation and one in negative"""
        self.alignment1 = Alignment("ctg2\t30\t8\t18\t+\tctg1\t30\t6\t15\t8\t9\t60\tcg:Z:1=1X2I3=1D3=\ttp:A:P")
        self.alignment2 = Alignment("ctg2\t30\t21\t27\t-\tctg1\t30\t20\t27\t5\t7\tcg:Z:3=2D1X1I1=\ttp:A:P")
        self.alignments = [self.alignment1, self.alignment2]
        self.contigLengths = {"ctg1":30, "ctg2":30, "ctg3":10}
        self.annotationBlockListsPerOrigContig = {"ctg1":{"annot1": BlockList([(4, 10),(14, 18)]),
                                                          "annot2": BlockList([(1, 3), (11, 13), (19, 30)])},
                                                  "ctg2":{"annot1": BlockList([(9, 12), (21, 27)]),
                                                          "annot2": BlockList([(1,8), (13, 20), (28, 30)])},
                                                  "ctg3":{"annot1": BlockList([(1,10)]),
                                                          "annot2": BlockList()}
                                                  }
        print(f"Tests:")

    def testCreatingHomologyRelationsDict(self):
        outputRelations = HomologyRelation.createAllInclusiveRelationDictFromAlignments(self.alignments, self.contigLengths, "_f")

        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 15, '+', "ctg1_f", 1)
        ctg1HomologyBlock3 = HomologyBlock("ctg1", 16, 20, '+', "ctg1_f", 2)
        ctg1HomologyBlock4 = HomologyBlock("ctg1", 21, 27, '+', "ctg1_f", 3)
        ctg1HomologyBlock5 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 4)

        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 8, '+', "ctg2_f", 0)
        ctg2HomologyBlock2 = HomologyBlock("ctg2", 9, 18, '+', "ctg2_f", 1)
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 19, 21, '+', "ctg2_f", 2)
        ctg2HomologyBlock4 = HomologyBlock("ctg2", 22, 27, '+', "ctg2_f", 3)
        ctg2HomologyBlock5 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f", 4)

        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)

        #blockAlignment1 = Alignment(f"query_ctg\t10\t0\t10\t+\tref_ctg\t9\t0\t9\t0\t0\t60\tcg:Z:1=1X2I3=1D3=\ttp:A:P")
        #blockAlignment2 = Alignment(f"query_ctg\t11\t0\t11\t+\tref_ctg\t10\t0\t10\t0\t0\t60\tcg:Z:3=2D1X1I1=\ttp:A:P")

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                               HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                               HomologyRelation(ctg1HomologyBlock3, None, None, None),
                               HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock4, getCigarList("3=2D1X1I1="), '-'),
                               HomologyRelation(ctg1HomologyBlock5, None, None, None)]
        truthRelations["ctg2_f"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                              HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                              HomologyRelation(ctg2HomologyBlock3, None, None, None),
                              HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock4, None, None),
                              HomologyRelation(ctg2HomologyBlock5, None, None, None)]
        truthRelations["ctg3_f"]= [HomologyRelation(ctg3HomologyBlock1, None, None, None)]

        for ctgName in ["ctg1_f", "ctg2_f", "ctg3_f"]:
            
            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")

                if truthRelation.alignment == None:
                    self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                else:
                    self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)
    def testInducingSwitchError(self):
        relationsDict = HomologyRelation.createAllInclusiveRelationDictFromAlignments(self.alignments, self.contigLengths, "_f")

        orderIndex = 3
        switchStart = 3
        switchEnd = 4
        HomologyRelation.induceSwitchErrorAndUpdateRelationsInNewContig(relationsDict, "ctg1_f", orderIndex, switchStart, switchEnd)

        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 15, '+', "ctg1_f", 1)
        ctg1HomologyBlock3 = HomologyBlock("ctg1", 16, 20, '+', "ctg1_f", 2)
        ctg1HomologyBlock4 = HomologyBlock("ctg1", 21, 22, '+', "ctg1_f", 3)
        ctg1HomologyBlock5 = HomologyBlock("ctg2", 25, 25, '-', "ctg1_f", 4)
        ctg1HomologyBlock6 = HomologyBlock("ctg1", 25, 27, '+', "ctg1_f", 5)
        ctg1HomologyBlock7 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 6)

        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 8, '+', "ctg2_f", 0)
        ctg2HomologyBlock2 = HomologyBlock("ctg2", 9, 18, '+', "ctg2_f", 1)
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 19, 21, '+', "ctg2_f", 2)
        ctg2HomologyBlock4 = HomologyBlock("ctg2", 22, 24, '+', "ctg2_f", 3)
        ctg2HomologyBlock5 = HomologyBlock("ctg1", 23, 24, '-', "ctg2_f", 4)
        ctg2HomologyBlock6 = HomologyBlock("ctg2", 26, 27, '+', "ctg2_f", 5)
        ctg2HomologyBlock7 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f", 6)

        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                  HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                                  HomologyRelation(ctg1HomologyBlock3, None, None, None),
                                  HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock6, getCigarList("2="), '-'),
                                  HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock5, getCigarList("1=1I"), '-'),
                                  HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock4, getCigarList("1D1X1I1="), '-'),
                                  HomologyRelation(ctg1HomologyBlock7, None, None, None)]

        truthRelations["ctg2_f"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                  HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                  HomologyRelation(ctg2HomologyBlock3, None, None, None),
                                  HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock6, None, None),
                                  HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock5, None, None),
                                  HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock4, None, None),
                                  HomologyRelation(ctg2HomologyBlock7, None, None, None)]

        truthRelations["ctg3_f"]= [HomologyRelation(ctg3HomologyBlock1, None, None, None)]

        for ctgName in ["ctg1_f", "ctg2_f", "ctg3_f"]:

            self.assertTrue(ctgName in relationsDict, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(relationsDict[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = relationsDict[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")

                if truthRelation.alignment == None:
                    self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                else:
                    self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testInducingSwitchErrorWithAnnotations(self):
        outputRelations = HomologyRelation.createAllInclusiveRelationDictFromAlignments(self.alignments, self.contigLengths, "_f")

        orderIndex = 3
        switchStart = 3
        switchEnd = 4
        HomologyRelation.induceSwitchErrorAndUpdateRelationsInNewContig(outputRelations, "ctg1_f", orderIndex, switchStart, switchEnd)

        # ctg1_f
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(4, 6)]),
                                                   "annot2": BlockList([(1,3)])}
        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 15, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 4), (8, 9)]),
                                                   "annot2": BlockList([(5, 7)])}
        ctg1HomologyBlock3 = HomologyBlock("ctg1", 16, 20, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1, 3)]),
                                                   "annot2": BlockList([(4, 5)])}
        ctg1HomologyBlock4 = HomologyBlock("ctg1", 21, 22, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 2)])}
        ctg1HomologyBlock5 = HomologyBlock("ctg2", 25, 25, '-', "ctg1_f", 4)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([1,1]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock6 = HomologyBlock("ctg1", 25, 27, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}
        ctg1HomologyBlock7 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 6)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}

        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 8, '+', "ctg2_f", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,8)])}
        ctg2HomologyBlock2 = HomologyBlock("ctg2", 9, 18, '+', "ctg2_f", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 4)]),
                                                   "annot2": BlockList([(5, 10)])}
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 19, 21, '+', "ctg2_f", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(3, 3)]),
                                                   "annot2": BlockList([(1, 2)])}
        ctg2HomologyBlock4 = HomologyBlock("ctg2", 22, 24, '+', "ctg2_f", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5 = HomologyBlock("ctg1", 23, 24, '-', "ctg2_f", 4)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,2)])}
        ctg2HomologyBlock6 = HomologyBlock("ctg2", 26, 27, '+', "ctg2_f", 5)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(1, 2)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock7 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f", 6)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,3)])}
        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([1,10]),
                                                   "annot2": BlockList([])}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock6, getCigarList("2="), '-'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock5, getCigarList("1=1I"), '-'),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock4, getCigarList("1D1X1I1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None)]

        truthRelations["ctg2_f"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                    HomologyRelation(ctg2HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock6, None, None),
                                    HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock5, None, None),
                                    HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock4, None, None),
                                    HomologyRelation(ctg2HomologyBlock7, None, None, None)]

        truthRelations["ctg3_f"]= [HomologyRelation(ctg3HomologyBlock1, None, None, None)]

        for ctgName in ["ctg1_f", "ctg2_f", "ctg3_f"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")

                else:
                    self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)


    def testFillingAnnotationBlockListsFromOriginalContigs(self):
        outputRelations = HomologyRelation.createAllInclusiveRelationDictFromAlignments(self.alignments,
                                                                                      self.contigLengths,
                                                                                      "_f")
        HomologyRelation.fillAnnotationBlockListsFromOriginalContigs(outputRelations,
                                                                     self.annotationBlockListsPerOrigContig,
                                                                     self.contigLengths,
                                                                     "_f")

        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(4, 6)]),
                                                   "annot2": BlockList([(1,3)])}
        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 15, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 4), (8, 9)]),
                                                   "annot2": BlockList([(5, 7)])}
        ctg1HomologyBlock3 = HomologyBlock("ctg1", 16, 20, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1, 3)]),
                                                   "annot2": BlockList([(4, 5)])}
        ctg1HomologyBlock4 = HomologyBlock("ctg1", 21, 27, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 7)])}
        ctg1HomologyBlock5 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}


        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 8, '+', "ctg2_f", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,8)])}
        ctg2HomologyBlock2 = HomologyBlock("ctg2", 9, 18, '+', "ctg2_f", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 4)]),
                                                   "annot2": BlockList([(5, 10)])}
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 19, 21, '+', "ctg2_f", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(3, 3)]),
                                                   "annot2": BlockList([(1, 2)])}
        ctg2HomologyBlock4 = HomologyBlock("ctg2", 22, 27, '+', "ctg2_f", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f", 4)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}

        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,10)]),
                                                   "annot2": BlockList([])}

        #blockAlignment1 = Alignment(f"query_ctg\t10\t0\t10\t+\tref_ctg\t9\t0\t9\t0\t0\t60\tcg:Z:1=1X2I3=1D3=\ttp:A:P")
        #blockAlignment2 = Alignment(f"query_ctg\t11\t0\t11\t+\tref_ctg\t10\t0\t10\t0\t0\t60\tcg:Z:3=2D1X1I1=\ttp:A:P")

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock4, getCigarList("3=2D1X1I1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock5, None, None, None)]
        truthRelations["ctg2_f"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                    HomologyRelation(ctg2HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock4, None, None),
                                    HomologyRelation(ctg2HomologyBlock5, None, None, None)]
        truthRelations["ctg3_f"]= [HomologyRelation(ctg3HomologyBlock1, None, None, None)]

        for ctgName in ["ctg1_f", "ctg2_f", "ctg3_f"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")

def main():
    unittest.main(verbosity=2)

if __name__ == '__main__':
    main()

