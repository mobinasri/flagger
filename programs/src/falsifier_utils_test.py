from block_utils import *
from falsifier_utils import *
import unittest
from project_blocks_multi_thread import makeCigarString

class TestProjection(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """Define alignments for testing falsifier functions"""
        self.alignment1ForMisjoin = Alignment("ctg2\t200\t99\t149\t+\tctg1\t200\t9\t57\t33\t48\t60\tcg:Z:4=1X2I2=1X10D1X2=10I1X1=5D5=1X4=5I3=1X6=\ttp:A:P")
        self.alignment2ForMisjoin = Alignment("ctg4\t200\t0\t37\t+\tctg3\t200\t0\t34\t31\t34\t60\tcg:Z:3=1X2D10=1X5I10=1I5=1D1X\ttp:A:P")
        self.annotationBlockListsPerOrigContigForMisjoin = {"ctg1":{"annot1": BlockList([(1, 28),(50, 200)]),
                                                                    "annot2": BlockList([(29,49)])},
                                                            "ctg2":{"annot1": BlockList([(1, 105), (150, 200)]),
                                                                    "annot2": BlockList([(106,149)])},
                                                            "ctg3":{"annot1": BlockList([(1,4),(101,200)]),
                                                                    "annot2": BlockList([(5,100)])},
                                                            "ctg4":{"annot1": BlockList([(1,50)]),
                                                                    "annot2": BlockList([(51,200)])}
                                                            }
        self.contigLengthsForMisjoin = {"ctg1":200, "ctg2":200, "ctg3":200, "ctg4":200}
        self.alignment1 = Alignment("ctg2\t30\t8\t18\t+\tctg1\t30\t6\t15\t8\t9\t60\tcg:Z:1=1X2I3=1D3=\ttp:A:P")
        self.alignment2 = Alignment("ctg2\t30\t21\t27\t-\tctg1\t30\t20\t27\t5\t7\t60\tcg:Z:3=2D1X1I1=\ttp:A:P")
        self.alignments = [self.alignment1, self.alignment2]
        self.contigLengths = {"ctg1":30, "ctg2":30, "ctg3":10}
        self.annotationBlockListsPerOrigContig = {"ctg1":{"annot1": BlockList([(4, 10),(14, 18)]),
                                                          "annot2": BlockList([(1, 3), (11, 13), (19, 30)])},
                                                  "ctg2":{"annot1": BlockList([(9, 12), (21, 27)]),
                                                          "annot2": BlockList([(1,8), (13, 20), (28, 30)])},
                                                  "ctg3":{"annot1": BlockList([(1,10)]),
                                                          "annot2": BlockList()}
                                                  }
        self.alignment3 = Alignment("ctg2\t68\t5\t65\t+\tctg1\t60\t6\t55\t40\t49\t60\tcg:Z:3=1X2=2I1=2D1X2=3I2=3D2=1X2=4I6=4I2=1X3=1I3=4I1=4D1=2I6=")
        self.annotationBlockListsPerOrigContigForAlignment3 = {"ctg1":{"annot1": BlockList([(1, 12),(30, 34), (51, 60)]),
                                                                       "annot2": BlockList([(13, 29), (35, 50)])},
                                                              "ctg2":{"annot1": BlockList([(1, 68)]),
                                                                      "annot2": BlockList([])}
                                                              }
        self.contigLengthsForAlignment3 = {"ctg1": 60, "ctg2": 68}
        self.contigSequencesForAlignment3 = {"ctg1": "TAAAAAGTGTGTCTCTCTCATATATGGGTACTACTACTACGAGTTTAAAGTAGTAGGGGT",
                                             "ctg2": "ACACAGTGAGTAACACTATACAATCGGATATTACTACGTAATAATACCGAGACCATATTGTAGTATGT"}
        # alignment 4 is same as alignment 3 but in negative orientation
        self.alignment4 = Alignment("ctg2\t68\t3\t63\t-\tctg1\t60\t6\t55\t40\t49\t60\tcg:Z:3=1X2=2I1=2D1X2=3I2=3D2=1X2=4I6=4I2=1X3=1I3=4I1=4D1=2I6=")
        self.annotationBlockListsPerOrigContigForAlignment4 = self.annotationBlockListsPerOrigContigForAlignment3
        self.contigLengthsForAlignment4 =  self.contigLengthsForAlignment3
        self.contigSequencesForAlignment4 = {"ctg1": "TAAAAAGTGTGTCTCTCTCATATATGGGTACTACTACTACGAGTTTAAAGTAGTAGGGGT",
                                             "ctg2": "ACATACTACAATATGGTCTCGGTATTATTACGTAGTAATATCCGATTGTATAGTGTTACTCACTGTGT"}
        print(f"Tests:")

    def testCreatingHomologyRelationsDict(self):
        outputRelationChains = HomologyRelationChains(self.alignments, self.contigLengths, ["ctg1"], "_f")

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

        outputRelations = outputRelationChains.relationChains
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

    def testInducingSwitchErrorWithAnnotations(self):
        outputRelationChains = HomologyRelationChains(self.alignments,
                                                      self.contigLengths,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContig,
                                                                         self.contigLengths,
                                                                         "_f")

        orderIndex = 3
        switchStart = 3
        switchEnd = 4
        switchEffectWindowLength = 1 # though this is not tested here but should be defined
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

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
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock6 = HomologyBlock("ctg1", 25, 27, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}
        ctg1HomologyBlock7 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
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
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1, 3)]),
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
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,10)]),
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

        outputRelations = outputRelationChains.relationChains
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


    def testInducingCollapseErrorWithAnnotations(self):
        outputRelationChains = HomologyRelationChains(self.alignments,
                                                      self.contigLengths,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContig,
                                                                         self.contigLengths,
                                                                         "_f")

        orderIndex = 3
        collapseStart = 3
        collapseEnd = 4
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)

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
        ctg1HomologyBlock5 = HomologyBlock("ctg1", 23, 24, '+', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock6 = HomologyBlock("ctg1", 25, 27, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}
        ctg1HomologyBlock7 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}

        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 8, '+', "ctg2_f.p_1_1", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,8)])}
        ctg2HomologyBlock2 = HomologyBlock("ctg2", 9, 18, '+', "ctg2_f.p_1_1", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 4)]),
                                                   "annot2": BlockList([(5, 10)])}
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 19, 21, '+', "ctg2_f.p_1_1", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(3, 3)]),
                                                   "annot2": BlockList([(1, 2)])}
        ctg2HomologyBlock4 = HomologyBlock("ctg2", 22, 24, '+', "ctg2_f.p_1_1", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1, 3)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5 = HomologyBlock("ctg2", 26, 27, '+', "ctg2_f.p_1_0", 0)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1, 2)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock6 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f.p_1_0", 1)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,3)])}

        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,10)]),
                                                   "annot2": BlockList([])}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock5, getCigarList("2="), '-'),
                                    HomologyRelation(ctg1HomologyBlock5, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock4, getCigarList("1D1X1I1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None)]

        
        truthRelations["ctg2_f.p_1_1"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                       HomologyRelation(ctg2HomologyBlock3, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock6, None, None)]

        truthRelations["ctg2_f.p_1_0"] = [HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock4, None, None),
                                       HomologyRelation(ctg2HomologyBlock6, None, None, None)]

        truthRelations["ctg3_f"]= [HomologyRelation(ctg3HomologyBlock1, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f.p_1_1", "ctg2_f.p_1_0", "ctg3_f"]:

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


    def testInducingDuplicationErrorWithAnnotations(self):
        outputRelationChains = HomologyRelationChains(self.alignments,
                                                      self.contigLengths,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContig,
                                                                         self.contigLengths,
                                                                         "_f")

        orderIndex = 3
        duplicationStart = 3
        duplicationEnd = 4
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

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
        ctg1HomologyBlock5 = HomologyBlock("ctg1", 23, 24, '+', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock6 = HomologyBlock("ctg1", 25, 27, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}
        ctg1HomologyBlock7 = HomologyBlock("ctg1", 28, 30, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1, 3)])}

        # ctg1.DUP_23_24
        ctg1DupHomologyBlock = HomologyBlock("ctg1", 23, 24, '+', "ctg1.DUP_23_24", 0)
        ctg1DupHomologyBlock.annotationBlockLists = {"annot1": BlockList([]),
                                                     "annot2": BlockList([(1,2)])}

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
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1, 3)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5 = HomologyBlock("ctg2", 25, 25, '+', "ctg2_f", 4)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock6 = HomologyBlock("ctg2", 26, 27, '+', "ctg2_f", 5)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(1, 2)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock7 = HomologyBlock("ctg2", 28, 30, '+', "ctg2_f", 6)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,3)])}

        # ctg3_f
        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 10, '+', "ctg3_f", 0)
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,10)]),
                                                   "annot2": BlockList([])}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("1=1X2I3=1D3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock6, getCigarList("2="), '-'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock5, getCigarList("1=1D"), '-'),
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
        truthRelations["ctg1.DUP_23_24"]= [HomologyRelation(ctg1DupHomologyBlock, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f", "ctg3_f", "ctg1.DUP_23_24"]:

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


    def testFillingAnnotationBlockListsFromOriginalContigs(self):
        outputRelationChains = HomologyRelationChains(self.alignments,
                                                      self.contigLengths,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContig,
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

        outputRelations = outputRelationChains.relationChains
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

    def testThreeMisAssembliesWithAnnotationsPositiveAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment3],
                                                      self.contigLengthsForAlignment3,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContigForAlignment3,
                                                                         self.contigLengthsForAlignment3,
                                                                         "_f")
        annotations = ["annot1", "annot2"]
        misAssemblyLength = 5
        minOverlapRatioWithEachAnnotation = 0.5 # min overlap length will be 3
        minMarginLength = 1 # margin should be greater than or equal to switchEffectWindowLength defined later
        outputRelationChains.updateAnnotationBlocksForSampling(annotations,
                                                                       misAssemblyLength,
                                                                       minOverlapRatioWithEachAnnotation,
                                                                       minMarginLength)



        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)

        # ctg1_f
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 19, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([(7,13)])}
        ctg1HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,4)]),
                                                                   "annot2": BlockList([(5,8)])}
        ctg1HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [3],
                                                                     "annot2": [4]}
        ctg1HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 3,
                                                                     "annot2": 4}

        ctg1HomologyBlock3 = HomologyBlock("ctg1", 20, 24, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock3.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock3.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock3.misAssemblyBlockList = BlockList([(1,5,"Dup")])



        ctg1HomologyBlock4 = HomologyBlock("ctg1", 25, 31, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(6,7)]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock4.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock4.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock4.misAssemblyBlockList = BlockList([(7,7,"Sw")])


        ctg1HomologyBlock5 = HomologyBlock("ctg2", 35, 43, '+', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock5.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (9,9,"Sw")])



        ctg1HomologyBlock6 = HomologyBlock("ctg1", 37, 43, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,7)])}
        ctg1HomologyBlock6.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock6.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock6.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        ctg1HomologyBlock7 = HomologyBlock("ctg1", 44, 48, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock7.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock7.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock7.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock7.misAssemblyBlockList = BlockList([(1,5,"Col_Del")])


        ctg1HomologyBlock8 = HomologyBlock("ctg1", 49, 55, '+', "ctg1_f", 7)
        ctg1HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(3,7)]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock8.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,2)]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock8.annotationStartBlockLengthsForSampling = {"annot1": [1],
                                                                     "annot2": []}
        ctg1HomologyBlock8.annotationStartTotalLengthsForSampling = {"annot1": 1,
                                                                     "annot2": 0}



        ctg1HomologyBlock9 = HomologyBlock("ctg1", 56, 60, '+', "ctg1_f", 8)
        ctg1HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock9.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        # ctg1.DUP_20_24
        # falsely duplicated block
        ctg1DupHomologyBlock = HomologyBlock("ctg1", 20, 24, '+', "ctg1.DUP_20_24", 0)
        ctg1DupHomologyBlock.annotationBlockLists = {"annot1": BlockList([]),
                                                     "annot2": BlockList([(1,5)])}
        ctg1DupHomologyBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                     "annot2": BlockList([])}
        ctg1DupHomologyBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                       "annot2": []}
        ctg1DupHomologyBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                       "annot2": 0}
        ctg1DupHomologyBlock.misAssemblyBlockList = BlockList([(1,5,"Dup")])


        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 5, '+', "ctg2_f.p_1_1", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock2 = HomologyBlock("ctg2", 6, 21, '+', "ctg2_f.p_1_1", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,16)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock3 = HomologyBlock("ctg2", 22, 23, '+', "ctg2_f.p_1_1", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,2)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock4 = HomologyBlock("ctg2", 24, 34, '+', "ctg2_f.p_1_1", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1,11)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock4.misAssemblyBlockList = BlockList([(11,11,"Sw")])

        ctg2HomologyBlock5 = HomologyBlock("ctg1", 32, 36, '+', "ctg2_f.p_1_1", 4)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([(4,5)])}
        ctg2HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (5,5,"Sw")])


        ctg2HomologyBlock6 = HomologyBlock("ctg2", 44, 55, '+', "ctg2_f.p_1_1", 5)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(1, 12)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        # the blocks after the collapsed block start here
        ctg2HomologyBlock7 = HomologyBlock("ctg2", 57, 65, '+', "ctg2_f.p_1_0", 0)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock8 = HomologyBlock("ctg2", 66, 68, '+', "ctg2_f.p_1_0", 1)
        ctg2HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([])}

        qBlocks = [ctg2HomologyBlock1, ctg2HomologyBlock2,
                   ctg2HomologyBlock3, ctg2HomologyBlock4,
                   ctg2HomologyBlock5, ctg2HomologyBlock6,
                   ctg2HomologyBlock7, ctg2HomologyBlock8]

        # query blocks will have no sampling blocks based on the current implementation
        for qBlock in qBlocks:
            qBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                           "annot2": BlockList([])}
            qBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                             "annot2": []}
            qBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                             "annot2": 0}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("3=1X2=2I1=2D1X2=3I1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, ctg2HomologyBlock3, getCigarList("1=3D1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock4, getCigarList("1=1X2=4I3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock5, getCigarList("3=4D2="), '+'),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock6, getCigarList("1X3=1I3=4I"), '+'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock8, ctg2HomologyBlock7, getCigarList("1=2I6="), '+'),
                                    HomologyRelation(ctg1HomologyBlock9, None, None, None)]

        # the relations before the collapsed block
        truthRelations["ctg2_f.p_1_1"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                       HomologyRelation(ctg2HomologyBlock3, ctg1HomologyBlock3, None, None),
                                       HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock4, None, None),
                                       HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock5, None, None),
                                       HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock6, None, None)]


        # the relations after the collapsed block
        truthRelations["ctg2_f.p_1_0"] = [HomologyRelation(ctg2HomologyBlock7, ctg1HomologyBlock8, None, None),
                                       HomologyRelation(ctg2HomologyBlock8, None, None, None)]

        truthRelations["ctg1.DUP_20_24"]= [HomologyRelation(ctg1DupHomologyBlock, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f.p_1_1", "ctg2_f.p_1_0", "ctg1.DUP_20_24"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        # check annotation blocks
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            #print(name, truthBlock.orderIndex,blockList.blocks, outputBlock.annotationBlockLists[name].blocks)
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")
                        # check start location blocks for sampling per annotation
                        for name, blockList in truthBlock.annotationStartBlockListsForSampling.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationStartBlockListsForSampling[name]), f"annotation list for sampling is not correct ({name})")
                        # check lengths of the sampling blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationStartBlockLengthsForSampling.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationStartBlockLengthsForSampling[name], f"annotation length list for sampling is not correct ({name})")
                        # check total length of sampling blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationStartTotalLengthsForSampling.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationStartTotalLengthsForSampling[name], f"total annotation length for sampling is not correct ({name})")

                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")
                        self.assertTrue(truthBlock.misAssemblyBlockList.isEqual(outputBlock.misAssemblyBlockList), f"mis-assembly block lists are not correct ({name})")

                    if truthRelation.alignment == None:
                        self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                    else:
                        self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testThreeMisAssembliesWithAnnotationsNegativeAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment4],
                                                      self.contigLengthsForAlignment4,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContigForAlignment4,
                                                                         self.contigLengthsForAlignment4,
                                                                         "_f")
        annotations = ["annot1", "annot2"]
        misAssemblyLength = 5
        minOverlapRatioWithEachAnnotation = 0.5 # min overlap length will be 3
        minMarginLength = 1 # margin should be greater than or equal to switchEffectWindowLength defined later
        outputRelationChains.updateAnnotationBlocksForSampling(annotations,
                                                                       misAssemblyLength,
                                                                       minOverlapRatioWithEachAnnotation,
                                                                       minMarginLength)



        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)

        # ctg1_f
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 19, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([(7,13)])}
        ctg1HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,4)]),
                                                                   "annot2": BlockList([(5,8)])}
        ctg1HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [3],
                                                                     "annot2": [4]}
        ctg1HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 3,
                                                                     "annot2": 4}

        ctg1HomologyBlock3 = HomologyBlock("ctg1", 20, 24, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock3.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock3.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock3.misAssemblyBlockList = BlockList([(1,5,"Dup")])



        ctg1HomologyBlock4 = HomologyBlock("ctg1", 25, 31, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(6,7)]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock4.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock4.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock4.misAssemblyBlockList = BlockList([(7,7,"Sw")])


        ctg1HomologyBlock5 = HomologyBlock("ctg2", 26, 34, '-', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock5.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (9,9,"Sw")])



        ctg1HomologyBlock6 = HomologyBlock("ctg1", 37, 43, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,7)])}
        ctg1HomologyBlock6.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock6.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock6.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        ctg1HomologyBlock7 = HomologyBlock("ctg1", 44, 48, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock7.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock7.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock7.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock7.misAssemblyBlockList = BlockList([(1,5,"Col_Del")])


        ctg1HomologyBlock8 = HomologyBlock("ctg1", 49, 55, '+', "ctg1_f", 7)
        ctg1HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(3,7)]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock8.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,2)]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock8.annotationStartBlockLengthsForSampling = {"annot1": [1],
                                                                     "annot2": []}
        ctg1HomologyBlock8.annotationStartTotalLengthsForSampling = {"annot1": 1,
                                                                     "annot2": 0}



        ctg1HomologyBlock9 = HomologyBlock("ctg1", 56, 60, '+', "ctg1_f", 8)
        ctg1HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock9.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        # ctg1.DUP_20_24
        # falsely duplicated block
        ctg1DupHomologyBlock = HomologyBlock("ctg1", 20, 24, '+', "ctg1.DUP_20_24", 0)
        ctg1DupHomologyBlock.annotationBlockLists = {"annot1": BlockList([]),
                                                     "annot2": BlockList([(1,5)])}
        ctg1DupHomologyBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                     "annot2": BlockList([])}
        ctg1DupHomologyBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                       "annot2": []}
        ctg1DupHomologyBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                       "annot2": 0}
        ctg1DupHomologyBlock.misAssemblyBlockList = BlockList([(1,5,"Dup")])


        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 3, '+', "ctg2_f.p_1_1", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock2 = HomologyBlock("ctg2", 4, 12, '+', "ctg2_f.p_1_1", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}

        # the blocks after the collapse block starts here
        ctg2HomologyBlock3 = HomologyBlock("ctg2", 14, 25, '+', "ctg2_f.p_1_0", 0)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,12)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock3.misAssemblyBlockList = BlockList([(12,12,"Sw")])

        ctg2HomologyBlock4 = HomologyBlock("ctg1", 32, 36, '-', "ctg2_f.p_1_0", 1)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(3,5)]),
                                                   "annot2": BlockList([(1,2)])}
        ctg2HomologyBlock4.misAssemblyBlockList = BlockList([(1,1,"Sw"), (5,5,"Sw")])

        ctg2HomologyBlock5 = HomologyBlock("ctg2", 35, 45, '+', "ctg2_f.p_1_0", 2)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,11)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        ctg2HomologyBlock6 = HomologyBlock("ctg2", 46, 47, '+', "ctg2_f.p_1_0", 3)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(1,2)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock7 = HomologyBlock("ctg2", 48, 63, '+', "ctg2_f.p_1_0", 4)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([(1,16)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock8 = HomologyBlock("ctg2", 64, 68, '+', "ctg2_f.p_1_0", 5)
        ctg2HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}

        qBlocks = [ctg2HomologyBlock1, ctg2HomologyBlock2,
                   ctg2HomologyBlock3, ctg2HomologyBlock4,
                   ctg2HomologyBlock5, ctg2HomologyBlock6,
                   ctg2HomologyBlock7, ctg2HomologyBlock8]

        # query blocks will have no sampling blocks based on the current implementation
        for qBlock in qBlocks:
            qBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                           "annot2": BlockList([])}
            qBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                             "annot2": []}
            qBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                             "annot2": 0}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock7, getCigarList("3=1X2=2I1=2D1X2=3I1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, ctg2HomologyBlock6, getCigarList("1=3D1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock5, getCigarList("1=1X2=4I3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock4, getCigarList("3=4D2="), '+'),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock3, getCigarList("1X3=1I3=4I"), '+'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock8, ctg2HomologyBlock2, getCigarList("1=2I6="), '+'),
                                    HomologyRelation(ctg1HomologyBlock9, None, None, None)]

        # the relations before the collapsed block
        truthRelations["ctg2_f.p_1_1"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock8, None, None)]

        # the relations after the collapsed block
        truthRelations["ctg2_f.p_1_0"] = [HomologyRelation(ctg2HomologyBlock3, ctg1HomologyBlock6, None, None),
                                       HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock5, None, None),
                                       HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock4, None, None),
                                       HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock3, None, None),
                                       HomologyRelation(ctg2HomologyBlock7, ctg1HomologyBlock2, None, None),
                                       HomologyRelation(ctg2HomologyBlock8, None, None, None)]

        truthRelations["ctg1.DUP_20_24"]= [HomologyRelation(ctg1DupHomologyBlock, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f.p_1_1", "ctg2_f.p_1_0", "ctg1.DUP_20_24"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        # check annotation blocks
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")
                        # check start location blocks for sampling per annotation
                        for name, blockList in truthBlock.annotationStartBlockListsForSampling.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationStartBlockListsForSampling[name]), f"annotation list for sampling is not correct ({name})")
                        # check lengths of the sampling blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationStartBlockLengthsForSampling.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationStartBlockLengthsForSampling[name], f"annotation length list for sampling is not correct ({name})")
                        # check total length of sampling blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationStartTotalLengthsForSampling.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationStartTotalLengthsForSampling[name], f"total annotation length for sampling is not correct ({name})")

                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")
                        self.assertTrue(truthBlock.misAssemblyBlockList.isEqual(outputBlock.misAssemblyBlockList), f"mis-assembly block lists are not correct ({name})")

                    if truthRelation.alignment == None:
                        self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                    else:
                        self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testFourMisAssembliesWithAnnotationsPositiveAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment3],
                                                      self.contigLengthsForAlignment3,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContigForAlignment3,
                                                                         self.contigLengthsForAlignment3,
                                                                         "_f")
        annotations = ["annot1", "annot2"]
        misAssemblyLength = 5
        minOverlapRatioWithEachAnnotation = 0.5 # min overlap length will be 3
        minMarginLength = 1 # margin should be greater than or equal to switchEffectWindowLength defined later
        outputRelationChains.updateAnnotationBlocksForSampling(annotations,
                                                                       misAssemblyLength,
                                                                       minOverlapRatioWithEachAnnotation,
                                                                       minMarginLength)



        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)

        # make an erroneous block
        orderIndex = 7
        errorStart = 3
        errorEnd = 6
        outputRelationChains.induceBaseErrorMisAssembly("ctg1_f",
                                                        orderIndex,
                                                        errorStart,
                                                        errorEnd)

        # ctg1_f
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 19, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([(7,13)])}
        ctg1HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,4)]),
                                                                   "annot2": BlockList([(5,8)])}
        ctg1HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [3],
                                                                     "annot2": [4]}
        ctg1HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 3,
                                                                     "annot2": 4}

        ctg1HomologyBlock3 = HomologyBlock("ctg1", 20, 24, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock3.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock3.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock3.misAssemblyBlockList = BlockList([(1,5,"Dup")])



        ctg1HomologyBlock4 = HomologyBlock("ctg1", 25, 31, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(6,7)]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock4.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock4.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock4.misAssemblyBlockList = BlockList([(7,7,"Sw")])


        ctg1HomologyBlock5 = HomologyBlock("ctg2", 35, 43, '+', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock5.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (9,9,"Sw")])



        ctg1HomologyBlock6 = HomologyBlock("ctg1", 37, 43, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,7)])}
        ctg1HomologyBlock6.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock6.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock6.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock6.misAssemblyBlockList = BlockList([(1,1, "Sw")])


        ctg1HomologyBlock7 = HomologyBlock("ctg1", 44, 48, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock7.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock7.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock7.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock7.misAssemblyBlockList = BlockList([(1,5,"Col_Del")])


        ctg1HomologyBlock8 = HomologyBlock("ctg1", 49, 50, '+', "ctg1_f", 7)
        ctg1HomologyBlock8.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock8.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock8.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock8.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}


        ctg1HomologyBlock9 = HomologyBlock("ctg1", 51, 54, '+', "ctg1_f", 8)
        ctg1HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock9.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock9.misAssemblyBlockList = BlockList([(1,4,"Err")])

        ctg1HomologyBlock10 = HomologyBlock("ctg1", 55, 55, '+', "ctg1_f", 9)
        ctg1HomologyBlock10.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                    "annot2": BlockList([])}
        ctg1HomologyBlock10.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                    "annot2": BlockList([])}
        ctg1HomologyBlock10.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                      "annot2": []}
        ctg1HomologyBlock10.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                      "annot2": 0}

        ctg1HomologyBlock11 = HomologyBlock("ctg1", 56, 60, '+', "ctg1_f", 10)
        ctg1HomologyBlock11.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                    "annot2": BlockList([])}
        ctg1HomologyBlock11.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                    "annot2": BlockList([])}
        ctg1HomologyBlock11.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                      "annot2": []}
        ctg1HomologyBlock11.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                      "annot2": 0}

        # ctg1.DUP_20_24
        # falsely duplicated block
        ctg1DupHomologyBlock = HomologyBlock("ctg1", 20, 24, '+', "ctg1.DUP_20_24", 0)
        ctg1DupHomologyBlock.annotationBlockLists = {"annot1": BlockList([]),
                                                     "annot2": BlockList([(1,5)])}
        ctg1DupHomologyBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                     "annot2": BlockList([])}
        ctg1DupHomologyBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                       "annot2": []}
        ctg1DupHomologyBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                       "annot2": 0}
        ctg1DupHomologyBlock.misAssemblyBlockList = BlockList([(1,5,"Dup")])


        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 5, '+', "ctg2_f.p_1_1", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock2 = HomologyBlock("ctg2", 6, 21, '+', "ctg2_f.p_1_1", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,16)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock3 = HomologyBlock("ctg2", 22, 23, '+', "ctg2_f.p_1_1", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,2)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock4 = HomologyBlock("ctg2", 24, 34, '+', "ctg2_f.p_1_1", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1,11)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock4.misAssemblyBlockList = BlockList([(11,11,"Sw")])

        ctg2HomologyBlock5 = HomologyBlock("ctg1", 32, 36, '+', "ctg2_f.p_1_1", 4)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([(4,5)])}
        ctg2HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (5,5,"Sw")])


        ctg2HomologyBlock6 = HomologyBlock("ctg2", 44, 55, '+', "ctg2_f.p_1_1", 5)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(1, 12)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        # the blocks after the collapsed block start here
        ctg2HomologyBlock7 = HomologyBlock("ctg2", 57, 60, '+', "ctg2_f.p_1_0", 0)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock8 = HomologyBlock("ctg2", 61, 64, '+', "ctg2_f.p_1_0", 1)
        ctg2HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock8.misAssemblyBlockList = BlockList([(1,4,"Col_Err")])

        ctg2HomologyBlock9 = HomologyBlock("ctg2", 65, 65, '+', "ctg2_f.p_1_0", 2)
        ctg2HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock10 = HomologyBlock("ctg2", 66, 68, '+', "ctg2_f.p_1_0", 3)
        ctg2HomologyBlock10.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([])}

        qBlocks = [ctg2HomologyBlock1, ctg2HomologyBlock2,
                   ctg2HomologyBlock3, ctg2HomologyBlock4,
                   ctg2HomologyBlock5, ctg2HomologyBlock6,
                   ctg2HomologyBlock7, ctg2HomologyBlock8]

        # query blocks will have no sampling blocks based on the current implementation
        for qBlock in qBlocks:
            qBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                           "annot2": BlockList([])}
            qBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                             "annot2": []}
            qBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                             "annot2": 0}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("3=1X2=2I1=2D1X2=3I1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock3, ctg2HomologyBlock3, getCigarList("1=3D1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock4, getCigarList("1=1X2=4I3="), '+'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock5, getCigarList("3=4D2="), '+'),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock6, getCigarList("1X3=1I3=4I"), '+'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock8, ctg2HomologyBlock7, getCigarList("1=2I1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock9, ctg2HomologyBlock8, None, None),
                                    HomologyRelation(ctg1HomologyBlock10, ctg2HomologyBlock9, getCigarList("1="), '+'),
                                    HomologyRelation(ctg1HomologyBlock11, None, None, None)]

        # the relations before the collapsed block
        truthRelations["ctg2_f.p_1_1"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                       HomologyRelation(ctg2HomologyBlock3, ctg1HomologyBlock3, None, None),
                                       HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock4, None, None),
                                       HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock5, None, None),
                                       HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock6, None, None)]


        # the relations after the collapsed block
        truthRelations["ctg2_f.p_1_0"] = [HomologyRelation(ctg2HomologyBlock7, ctg1HomologyBlock8, None, None),
                                       HomologyRelation(ctg2HomologyBlock8, ctg1HomologyBlock9, None, None),
                                       HomologyRelation(ctg2HomologyBlock9, ctg1HomologyBlock10, None, None),
                                       HomologyRelation(ctg2HomologyBlock10, None, None, None)]

        truthRelations["ctg1.DUP_20_24"]= [HomologyRelation(ctg1DupHomologyBlock, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f.p_1_1", "ctg2_f.p_1_0", "ctg1.DUP_20_24"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        # check annotation blocks
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            #print(name, truthBlock.orderIndex,blockList.blocks, outputBlock.annotationBlockLists[name].blocks)
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")
                        # check start location blocks for sampling per annotation
                        for name, blockList in truthBlock.annotationStartBlockListsForSampling.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationStartBlockListsForSampling[name]), f"annotation list for sampling is not correct ({name})")
                        # check lengths of the sampling blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationStartBlockLengthsForSampling.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationStartBlockLengthsForSampling[name], f"annotation length list for sampling is not correct ({name})")
                        # check total length of sampling blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationStartTotalLengthsForSampling.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationStartTotalLengthsForSampling[name], f"total annotation length for sampling is not correct ({name})")

                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")
                        self.assertTrue(truthBlock.misAssemblyBlockList.isEqual(outputBlock.misAssemblyBlockList), f"mis-assembly block lists are not correct ({name})")

                    if truthRelation.alignment == None:
                        self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                    else:
                        self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testFourMisAssembliesWithAnnotationsNegativeAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment4],
                                                      self.contigLengthsForAlignment4,
                                                      ["ctg1"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContigForAlignment4,
                                                                         self.contigLengthsForAlignment4,
                                                                         "_f")
        annotations = ["annot1", "annot2"]
        misAssemblyLength = 5
        minOverlapRatioWithEachAnnotation = 0.5 # min overlap length will be 3
        minMarginLength = 1 # margin should be greater than or equal to switchEffectWindowLength defined later
        outputRelationChains.updateAnnotationBlocksForSampling(annotations,
                                                                       misAssemblyLength,
                                                                       minOverlapRatioWithEachAnnotation,
                                                                       minMarginLength)



        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)

        # make an erroneous block
        orderIndex = 7
        errorStart = 3
        errorEnd = 6
        outputRelationChains.induceBaseErrorMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       errorStart,
                                                       errorEnd)

        # ctg1_f
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 6, '+', "ctg1_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        ctg1HomologyBlock2 = HomologyBlock("ctg1", 7, 19, '+', "ctg1_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 6)]),
                                                   "annot2": BlockList([(7,13)])}
        ctg1HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,4)]),
                                                                   "annot2": BlockList([(5,8)])}
        ctg1HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [3],
                                                                     "annot2": [4]}
        ctg1HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 3,
                                                                     "annot2": 4}

        ctg1HomologyBlock3 = HomologyBlock("ctg1", 20, 24, '+', "ctg1_f", 2)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock3.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock3.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock3.misAssemblyBlockList = BlockList([(1,5,"Dup")])



        ctg1HomologyBlock4 = HomologyBlock("ctg1", 25, 31, '+', "ctg1_f", 3)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(6,7)]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock4.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock4.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock4.misAssemblyBlockList = BlockList([(7,7,"Sw")])


        ctg1HomologyBlock5 = HomologyBlock("ctg2", 26, 34, '-', "ctg1_f", 4)
        ctg1HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock5.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock5.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock5.misAssemblyBlockList = BlockList([(1,1,"Sw"), (9,9,"Sw")])



        ctg1HomologyBlock6 = HomologyBlock("ctg1", 37, 43, '+', "ctg1_f", 5)
        ctg1HomologyBlock6.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,7)])}
        ctg1HomologyBlock6.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,2)])}
        ctg1HomologyBlock6.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [1]}
        ctg1HomologyBlock6.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 1}
        ctg1HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        ctg1HomologyBlock7 = HomologyBlock("ctg1", 44, 48, '+', "ctg1_f", 6)
        ctg1HomologyBlock7.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,5)])}
        ctg1HomologyBlock7.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock7.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock7.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock7.misAssemblyBlockList = BlockList([(1,5,"Col_Del")])


        ctg1HomologyBlock8 = HomologyBlock("ctg1", 49, 50, '+', "ctg1_f", 7)
        ctg1HomologyBlock8.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,2)])}
        ctg1HomologyBlock8.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock8.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock8.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}


        ctg1HomologyBlock9 = HomologyBlock("ctg1", 51, 54, '+', "ctg1_f", 8)
        ctg1HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock9.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock9.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock9.misAssemblyBlockList = BlockList([(1,4,"Err")])

        ctg1HomologyBlock10 = HomologyBlock("ctg1", 55, 55, '+', "ctg1_f", 9)
        ctg1HomologyBlock10.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock10.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock10.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg1HomologyBlock10.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}

        ctg1HomologyBlock11 = HomologyBlock("ctg1", 56, 60, '+', "ctg1_f", 10)
        ctg1HomologyBlock11.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock11.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                    "annot2": BlockList([])}
        ctg1HomologyBlock11.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                      "annot2": []}
        ctg1HomologyBlock11.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                      "annot2": 0}

        # ctg1.DUP_20_24
        # falsely duplicated block
        ctg1DupHomologyBlock = HomologyBlock("ctg1", 20, 24, '+', "ctg1.DUP_20_24", 0)
        ctg1DupHomologyBlock.annotationBlockLists = {"annot1": BlockList([]),
                                                     "annot2": BlockList([(1,5)])}
        ctg1DupHomologyBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                     "annot2": BlockList([])}
        ctg1DupHomologyBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                       "annot2": []}
        ctg1DupHomologyBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                       "annot2": 0}
        ctg1DupHomologyBlock.misAssemblyBlockList = BlockList([(1,5,"Dup")])


        # ctg2_f
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 3, '+', "ctg2_f.p_1_1", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,3)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock2 = HomologyBlock("ctg2", 4, 4, '+', "ctg2_f.p_1_1", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,1)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock3 = HomologyBlock("ctg2", 5, 8, '+', "ctg2_f.p_1_1", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock3.misAssemblyBlockList = BlockList([(1,4,"Col_Err")])

        ctg2HomologyBlock4 = HomologyBlock("ctg2", 9, 12, '+', "ctg2_f.p_1_1", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}

        # the blocks after the collapse block starts here
        ctg2HomologyBlock5 = HomologyBlock("ctg2", 14, 25, '+', "ctg2_f.p_1_0", 0)
        ctg2HomologyBlock5.annotationBlockLists = {"annot1": BlockList([(1,12)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock5.misAssemblyBlockList = BlockList([(12,12,"Sw")])

        ctg2HomologyBlock6 = HomologyBlock("ctg1", 32, 36, '-', "ctg2_f.p_1_0", 1)
        ctg2HomologyBlock6.annotationBlockLists = {"annot1": BlockList([(3,5)]),
                                                   "annot2": BlockList([(1,2)])}
        ctg2HomologyBlock6.misAssemblyBlockList = BlockList([(1,1,"Sw"), (5,5,"Sw")])

        ctg2HomologyBlock7 = HomologyBlock("ctg2", 35, 45, '+', "ctg2_f.p_1_0", 2)
        ctg2HomologyBlock7.annotationBlockLists = {"annot1": BlockList([(1,11)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock7.misAssemblyBlockList = BlockList([(1,1,"Sw")])


        ctg2HomologyBlock8 = HomologyBlock("ctg2", 46, 47, '+', "ctg2_f.p_1_0", 3)
        ctg2HomologyBlock8.annotationBlockLists = {"annot1": BlockList([(1,2)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock9 = HomologyBlock("ctg2", 48, 63, '+', "ctg2_f.p_1_0", 4)
        ctg2HomologyBlock9.annotationBlockLists = {"annot1": BlockList([(1,16)]),
                                                   "annot2": BlockList([])}

        ctg2HomologyBlock10 = HomologyBlock("ctg2", 64, 68, '+', "ctg2_f.p_1_0", 5)
        ctg2HomologyBlock10.annotationBlockLists = {"annot1": BlockList([(1,5)]),
                                                    "annot2": BlockList([])}

        qBlocks = [ctg2HomologyBlock1, ctg2HomologyBlock2,
                   ctg2HomologyBlock3, ctg2HomologyBlock4,
                   ctg2HomologyBlock5, ctg2HomologyBlock6,
                   ctg2HomologyBlock7, ctg2HomologyBlock8,
                   ctg2HomologyBlock9, ctg2HomologyBlock10]

        # query blocks will have no sampling blocks based on the current implementation
        for qBlock in qBlocks:
            qBlock.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                           "annot2": BlockList([])}
            qBlock.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                             "annot2": []}
            qBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                             "annot2": 0}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock9, getCigarList("3=1X2=2I1=2D1X2=3I1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock3, ctg2HomologyBlock8, getCigarList("1=3D1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock4, ctg2HomologyBlock7, getCigarList("1=1X2=4I3="), '-'),
                                    HomologyRelation(ctg1HomologyBlock5, ctg2HomologyBlock6, getCigarList("3=4D2="), '-'),
                                    HomologyRelation(ctg1HomologyBlock6, ctg2HomologyBlock5, getCigarList("1X3=1I3=4I"), '-'),
                                    HomologyRelation(ctg1HomologyBlock7, None, None, None),
                                    HomologyRelation(ctg1HomologyBlock8, ctg2HomologyBlock4, getCigarList("1=2I1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock9, ctg2HomologyBlock3, None, None),
                                    HomologyRelation(ctg1HomologyBlock10, ctg2HomologyBlock2, getCigarList("1="), '-'),
                                    HomologyRelation(ctg1HomologyBlock11, None, None, None)]

        # the relations before the collapsed block
        truthRelations["ctg2_f.p_1_1"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                       HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock10, None, None),
                                       HomologyRelation(ctg2HomologyBlock3, ctg1HomologyBlock9, None, None),
                                       HomologyRelation(ctg2HomologyBlock4, ctg1HomologyBlock8, None, None)]

        # the relations after the collapsed block
        truthRelations["ctg2_f.p_1_0"] = [HomologyRelation(ctg2HomologyBlock5, ctg1HomologyBlock6, None, None),
                                       HomologyRelation(ctg2HomologyBlock6, ctg1HomologyBlock5, None, None),
                                       HomologyRelation(ctg2HomologyBlock7, ctg1HomologyBlock4, None, None),
                                       HomologyRelation(ctg2HomologyBlock8, ctg1HomologyBlock3, None, None),
                                       HomologyRelation(ctg2HomologyBlock9, ctg1HomologyBlock2, None, None),
                                       HomologyRelation(ctg2HomologyBlock10, None, None, None)]

        truthRelations["ctg1.DUP_20_24"]= [HomologyRelation(ctg1DupHomologyBlock, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f", "ctg2_f.p_1_1", "ctg2_f.p_1_0", "ctg1.DUP_20_24"]:

            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), "Number of relations do not match")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        # check annotation blocks
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")
                        # check start location blocks for sampling per annotation
                        for name, blockList in truthBlock.annotationStartBlockListsForSampling.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationStartBlockListsForSampling[name]), f"annotation list for sampling is not correct ({name})")
                        # check lengths of the sampling blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationStartBlockLengthsForSampling.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationStartBlockLengthsForSampling[name], f"annotation length list for sampling is not correct ({name})")
                        # check total length of sampling blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationStartTotalLengthsForSampling.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationStartTotalLengthsForSampling[name], f"total annotation length for sampling is not correct ({name})")

                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")
                        self.assertTrue(truthBlock.misAssemblyBlockList.isEqual(outputBlock.misAssemblyBlockList), f"mis-assembly block lists are not correct ({name})")

                    if truthRelation.alignment == None:
                        self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                    else:
                        self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testMisjoinWithAnnotationsPositiveAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment1ForMisjoin, self.alignment2ForMisjoin],
                                                      self.contigLengthsForMisjoin,
                                                      ["ctg1", "ctg3"],
                                                      "_f")
        outputRelationChains.fillAnnotationBlockListsFromOriginalContigs(self.annotationBlockListsPerOrigContigForMisjoin,
                                                                         self.contigLengthsForMisjoin,
                                                                         "_f")
        annotations = ["annot1", "annot2"]
        misAssemblyLength = 5
        minOverlapRatioWithEachAnnotation = 0.5 # min overlap length will be 3
        minMarginLength = 1 # margin should be greater than or equal to switchEffectWindowLength defined later
        outputRelationChains.updateAnnotationBlocksForSampling(annotations,
                                                               misAssemblyLength,
                                                               minOverlapRatioWithEachAnnotation,
                                                               minMarginLength)



        # make a misjoin
        orderIndex_1 = 1
        loc_1 = 9
        orderIndex_2 = 0
        loc_2 = 5
        misjoinEffectWindowLength = 1
        outputRelationChains.induceMisjoinMisAssembly("ctg1_f",
                                                     orderIndex_1,
                                                     loc_1,
                                                     "ctg3_f",
                                                     orderIndex_2,
                                                     loc_2,
                                                     misjoinEffectWindowLength)

        # ctg1
        ctg1HomologyBlock1 = HomologyBlock("ctg1", 1, 9, '+', "ctg1_f.Msj_ctg3_f", 0)
        ctg1HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1, 9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock1.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 0}

        ctg1HomologyBlock2 = HomologyBlock("ctg1", 10, 18, '+', "ctg1_f.Msj_ctg3_f", 1)
        ctg1HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1, 9)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,4)]),
                                                                   "annot2": BlockList([])}
        ctg1HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [3],
                                                                     "annot2": []}
        ctg1HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 3,
                                                                     "annot2": 0}
        ctg1HomologyBlock2.annotationBlockListsForSamplingMisjoin = {"annot1": BlockList([(2,8)]),
                                                                     "annot2": BlockList([])}
        ctg1HomologyBlock2.annotationBlockLengthsForSamplingMisjoin = {"annot1": [7],
                                                                       "annot2": []}
        ctg1HomologyBlock2.annotationTotalLengthsForSamplingMisjoin = {"annot1": 7,
                                                                       "annot2": 0}
        ctg1HomologyBlock2.misAssemblyBlockList = BlockList([(9,9,"Msj")])

        ctg1HomologyBlock3 = HomologyBlock("ctg1", 19, 57, '+', "ctg3_f.Msj_ctg1_f", 1)
        ctg1HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,10),(32,39)]),
                                                   "annot2": BlockList([(11,31)])}
        ctg1HomologyBlock3.annotationStartBlockListsForSampling = {"annot1": BlockList([(2,8), (30,34)]),
                                                                   "annot2": BlockList([(9,29)])}
        ctg1HomologyBlock3.annotationStartBlockLengthsForSampling = {"annot1": [7,5],
                                                                     "annot2": [21]}
        ctg1HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 12,
                                                                     "annot2": 21}
        ctg1HomologyBlock3.annotationBlockListsForSamplingMisjoin = {"annot1": BlockList([(2,10), (32,38)]),
                                                                   "annot2": BlockList([(11,31)])}
        ctg1HomologyBlock3.annotationBlockLengthsForSamplingMisjoin = {"annot1": [9,7],
                                                                     "annot2": [21]}
        ctg1HomologyBlock3.annotationTotalLengthsForSamplingMisjoin = {"annot1": 16,
                                                                     "annot2": 21}
        ctg1HomologyBlock3.misAssemblyBlockList = BlockList([(1,1,"Msj")])



        ctg1HomologyBlock4 = HomologyBlock("ctg1", 58, 200, '+', "ctg3_f.Msj_ctg1_f", 2)
        ctg1HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1,143)]),
                                                   "annot2": BlockList([])}
        ctg1HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg1HomologyBlock4.annotationTotalLengthsForSamplingMisnoin = {"annot1": 0,
                                                                       "annot2": 0}

        # ctg2
        ctg2HomologyBlock1 = HomologyBlock("ctg2", 1, 99, '+', "ctg2_f", 0)
        ctg2HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,99)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg2HomologyBlock1.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 0}


        ctg2HomologyBlock2 = HomologyBlock("ctg2", 100, 109, '+', "ctg2_f", 1)
        ctg2HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,6)]),
                                                   "annot2": BlockList([(7,10)])}
        ctg2HomologyBlock2.misAssemblyBlockList = BlockList([(10,10,"Msj")])


        ctg2HomologyBlock3 = HomologyBlock("ctg2", 110, 149, '+', "ctg2_f", 2)
        ctg2HomologyBlock3.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,40)])}
        ctg2HomologyBlock3.misAssemblyBlockList = BlockList([(1,1,"Msj")])


        ctg2HomologyBlock4 = HomologyBlock("ctg2", 150, 200, '+', "ctg2_f", 3)
        ctg2HomologyBlock4.annotationBlockLists = {"annot1": BlockList([(1,51)]),
                                                   "annot2": BlockList([])}
        ctg2HomologyBlock4.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg2HomologyBlock4.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 0}

        # ctg3
        ctg3HomologyBlock1 = HomologyBlock("ctg3", 1, 5, '+', "ctg3_f.Msj_ctg1_f", 0)
        ctg3HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([(5,5)])}
        ctg3HomologyBlock1.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([])}
        ctg3HomologyBlock1.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": []}
        ctg3HomologyBlock1.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg3HomologyBlock1.annotationBlockListsForSamplingMisjoin = {"annot1": BlockList([(2,4)]),
                                                                     "annot2": BlockList([])}
        ctg3HomologyBlock1.annotationBlockLengthsForSamplingMisjoin = {"annot1": [3],
                                                                       "annot2": []}
        ctg3HomologyBlock1.annotationTotalLengthsForSamplingMisjoin = {"annot1": 3,
                                                                       "annot2": 0}
        ctg3HomologyBlock1.misAssemblyBlockList = BlockList([(5,5,"Msj")])


        ctg3HomologyBlock2 = HomologyBlock("ctg3", 6, 34, '+', "ctg1_f.Msj_ctg3_f", 2)
        ctg3HomologyBlock2.annotationBlockLists = {"annot1": BlockList([]),
                                                   "annot2": BlockList([(1,29)])}
        ctg3HomologyBlock2.annotationStartBlockListsForSampling = {"annot1": BlockList([]),
                                                                   "annot2": BlockList([(2,24)])}
        ctg3HomologyBlock2.annotationStartBlockLengthsForSampling = {"annot1": [],
                                                                     "annot2": [23]}
        ctg3HomologyBlock2.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 23}
        ctg3HomologyBlock2.annotationBlockListsForSamplingMisjoin = {"annot1": BlockList([]),
                                                                     "annot2": BlockList([(2,28)])}
        ctg3HomologyBlock2.annotationBlockLengthsForSamplingMisjoin = {"annot1": [],
                                                                       "annot2": [27]}
        ctg3HomologyBlock2.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 27}
        ctg3HomologyBlock2.misAssemblyBlockList = BlockList([(1,1,"Msj")])


        ctg3HomologyBlock3 = HomologyBlock("ctg3", 35, 200, '+', "ctg1_f.Msj_ctg3_f", 3)
        ctg3HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(67,166)]),
                                                   "annot2": BlockList([(1,66)])}
        ctg3HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg3HomologyBlock3.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 0}

        # ctg4
        ctg4HomologyBlock1 = HomologyBlock("ctg4", 1, 4, '+', "ctg4_f", 0)
        ctg4HomologyBlock1.annotationBlockLists = {"annot1": BlockList([(1,4)]),
                                                   "annot2": BlockList([])}
        ctg4HomologyBlock1.misAssemblyBlockList = BlockList([(4,4,"Msj")])


        ctg4HomologyBlock2 = HomologyBlock("ctg4", 5, 37, '+', "ctg4_f", 1)
        ctg4HomologyBlock2.annotationBlockLists = {"annot1": BlockList([(1,33)]),
                                                   "annot2": BlockList([])}
        ctg4HomologyBlock2.misAssemblyBlockList = BlockList([(1,1,"Msj")])


        ctg4HomologyBlock3 = HomologyBlock("ctg4", 38, 200, '+', "ctg4_f", 2)
        ctg4HomologyBlock3.annotationBlockLists = {"annot1": BlockList([(1,13)]),
                                                   "annot2": BlockList([(14,163)])}
        ctg4HomologyBlock3.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                                     "annot2": 0}
        ctg4HomologyBlock3.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                                       "annot2": 0}



        qBlocks = [ctg2HomologyBlock1, ctg2HomologyBlock2,
                   ctg2HomologyBlock3, ctg2HomologyBlock4,
                   ctg4HomologyBlock1, ctg4HomologyBlock2,
                   ctg4HomologyBlock3]

        # query blocks will have no sampling blocks in the current implementation
        for qBlock in qBlocks:
            qBlock.annotationStartTotalLengthsForSampling = {"annot1": 0,
                                                             "annot2": 0}
            qBlock.annotationTotalLengthsForSamplingMisjoin = {"annot1": 0,
                                                               "annot2": 0}

        truthRelations = defaultdict(list)
        truthRelations["ctg1_f.Msj_ctg3_f"] = [HomologyRelation(ctg1HomologyBlock1, None, None, None),
                                               HomologyRelation(ctg1HomologyBlock2, ctg2HomologyBlock2, getCigarList("4=1X2I2=1X1D"), '+'),
                                               HomologyRelation(ctg3HomologyBlock2, ctg4HomologyBlock2, getCigarList("1D10=1X5I10=1I5=1D1X"), '+'),
                                               HomologyRelation(ctg3HomologyBlock3, None, None, None)] 

        # the relations before the collapsed block
        truthRelations["ctg3_f.Msj_ctg1_f"] = [HomologyRelation(ctg3HomologyBlock1, ctg4HomologyBlock1, getCigarList("3=1X1D"), '+'),
                                               HomologyRelation(ctg1HomologyBlock3, ctg2HomologyBlock3, getCigarList("9D1X2=10I1X1=5D5=1X4=5I3=1X6="), '+'),
                                               HomologyRelation(ctg1HomologyBlock4, None, None, None)]

        # the relations after the collapsed block
        truthRelations["ctg2_f"] = [HomologyRelation(ctg2HomologyBlock1, None, None, None),
                                    HomologyRelation(ctg2HomologyBlock2, ctg1HomologyBlock2, None, None),
                                    HomologyRelation(ctg2HomologyBlock3, ctg1HomologyBlock3, None, None),
                                    HomologyRelation(ctg2HomologyBlock4, None, None, None)]

        truthRelations["ctg4_f"] = [HomologyRelation(ctg4HomologyBlock1, ctg3HomologyBlock1, None, None),
                                    HomologyRelation(ctg4HomologyBlock2, ctg3HomologyBlock2, None, None),
                                    HomologyRelation(ctg4HomologyBlock3, None, None, None)]

        outputRelations = outputRelationChains.relationChains
        for ctgName in ["ctg1_f.Msj_ctg3_f", "ctg3_f.Msj_ctg1_f", "ctg2_f", "ctg4_f"]:
            self.assertTrue(ctgName in outputRelations, "Contig does not exist")
            self.assertEqual(len(truthRelations[ctgName]), len(outputRelations[ctgName]), f"Number of relations do not match {ctgName}")

            for i in range(len(truthRelations[ctgName])):
                truthRelation = truthRelations[ctgName][i]
                outputRelation = outputRelations[ctgName][i]

                for truthBlock, outputBlock, name in zip([truthRelation.block, truthRelation.homologousBlock], [outputRelation.block, outputRelation.homologousBlock], ["block", "homologous block"] ):
                    if truthBlock == None:
                        self.assertTrue(outputBlock == None, f"the block is not None ({name})")
                    else:
                        # check annotation blocks
                        for name, blockList in truthBlock.annotationBlockLists.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockLists[name]), f"annotation list is not correct ({name})")
                        # check start location blocks for sampling per annotation
                        for name, blockList in truthBlock.annotationStartBlockListsForSampling.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationStartBlockListsForSampling[name]), f"annotation list for sampling is not correct ({name})")
                        # check lengths of the sampling blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationStartBlockLengthsForSampling.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationStartBlockLengthsForSampling[name], f"annotation length list for sampling is not correct ({name})")
                        # check total length of sampling blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationStartTotalLengthsForSampling.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationStartTotalLengthsForSampling[name], f"total annotation length for sampling is not correct ({name})")

                        # check start location blocks for sampling misjoin per annotation
                        for name, blockList in truthBlock.annotationBlockListsForSamplingMisjoin.items():
                            self.assertTrue(blockList.isEqual(outputBlock.annotationBlockListsForSamplingMisjoin[name]), f"annotation list for sampling misjoin is not correct ({name})")
                        # check lengths of the sampling misjoin blocks per annotation
                        for name, blocksLengthList in truthBlock.annotationBlockLengthsForSamplingMisjoin.items():
                            self.assertListEqual(blocksLengthList, outputBlock.annotationBlockLengthsForSamplingMisjoin[name], f"annotation length list for sampling misjoin is not correct ({name})")
                        # check total length of sampling misjoin blocks per annotation
                        for name, blocksTotalLength in truthBlock.annotationTotalLengthsForSamplingMisjoin.items():
                            self.assertEqual(blocksTotalLength, outputBlock.annotationTotalLengthsForSamplingMisjoin[name], f"total annotation length for sampling misjoin is not correct ({name})")


                        self.assertEqual(truthBlock.origCtg, outputBlock.origCtg, f"origCtg is not correct ({name})")
                        self.assertEqual(truthBlock.origStart, outputBlock.origStart, f"origStart is not correct ({name})")
                        self.assertEqual(truthBlock.origEnd, outputBlock.origEnd, f"origEnd is not correct ({name})")
                        self.assertEqual(truthBlock.origStrand, outputBlock.origStrand, f"origStrand is not correct ({name})")
                        self.assertEqual(truthBlock.newCtg, outputBlock.newCtg, f"newCtg is not correct ({name})")
                        self.assertEqual(truthBlock.orderIndex, outputBlock.orderIndex, f"orderIndex is not correct ({name})")
                        self.assertTrue(truthBlock.misAssemblyBlockList.isEqual(outputBlock.misAssemblyBlockList), f"mis-assembly block lists are not correct ({name})")

                    if truthRelation.alignment == None:
                        self.assertTrue(outputRelation.alignment == None, f"the alignment is not None")
                    else:
                        self.assertListEqual(truthRelation.alignment.cigarList, outputRelation.alignment.cigarList)

    def testSequenceThreeMisAssembliesPositiveAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment3],
                                                      self.contigLengthsForAlignment3,
                                                      ["ctg1"],
                                                      "_f")
        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)


        truthSequences = {"ctg1_f": "TAAAAAGTGTGTCTCTCTCATATATGGGTACTACGTAATACTACGAGTTTAAAGTAGTAGGGGT",
                          "ctg2_f.p_1_1": "ACACAGTGAGTAACACTATACAATCGGATATTACTACTAATACCGAGACCA",
                          "ctg2_f.p_1_0": "ATTGTAGTATGT",
                          "ctg1.DUP_20_24": "ATATA"}

        for newCtg, newSeq in outputRelationChains.yieldNewCtgSequences(self.contigSequencesForAlignment3, singleBaseErrorRate=0.0):
            self.assertEqual(truthSequences[newCtg], newSeq, f"sequence of {newCtg} is not correct")


    def testSequenceThreeMisAssembliesNegativeAlignment(self):
        outputRelationChains = HomologyRelationChains([self.alignment4],
                                                      self.contigLengthsForAlignment4,
                                                      ["ctg1"],
                                                      "_f")
        # make a switch error
        orderIndex = 1
        switchStart = 26
        switchEnd = 30
        switchEffectWindowLength = 1
        outputRelationChains.induceSwitchMisAssembly("ctg1_f",
                                                     orderIndex,
                                                     switchStart,
                                                     switchEnd,
                                                     switchEffectWindowLength)

        # make a false duplication
        orderIndex = 1
        duplicationStart = 14
        duplicationEnd = 18
        outputRelationChains.induceDuplicationMisAssembly("ctg1_f",
                                                          orderIndex,
                                                          duplicationStart,
                                                          duplicationEnd)

        # make a collapsed block
        orderIndex = 5
        collapseStart = 8
        collapseEnd = 12
        outputRelationChains.induceCollapseMisAssembly("ctg1_f",
                                                       orderIndex,
                                                       collapseStart,
                                                       collapseEnd)


        truthSequences = {"ctg1_f": "TAAAAAGTGTGTCTCTCTCATATATGGGTACTACGTAATACTACGAGTTTAAAGTAGTAGGGGT",
                          "ctg2_f.p_1_1": "ACATACTACAAT",
                          "ctg2_f.p_1_0": "TGGTCTCGGTATTAGTAGTAATATCCGATTGTATAGTGTTACTCACTGTGT",
                          "ctg1.DUP_20_24": "ATATA"}

        for newCtg, newSeq in outputRelationChains.yieldNewCtgSequences(self.contigSequencesForAlignment4, singleBaseErrorRate=0.0):
            self.assertEqual(truthSequences[newCtg], newSeq, f"sequence of {newCtg} is not correct")




def main():
    unittest.main(verbosity=2)

if __name__ == '__main__':
    main()

