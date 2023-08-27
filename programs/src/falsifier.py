import sys
import argparse
from block_utils import Alignment
from collections import defaultdict
from copy import deepcopy
from falsifier_utils import *
import numpy as np
import pandas as pd
import json
import re
import os


def induceMultipleMisAssembliesOfTheSameType(relationChains, annotation, misAssemblyType, misAssemblySize, misAssemblyCount, switchEffectWindowSize):
    assert(misAssemblyType in ["Sw", "Err", "Dup", "Col"])
    for i in range(misAssemblyCount):
        newCtg, orderIndex, start, end = relationChains.getRandomMisAssemblyInterval(annotation, misAssemblySize)
        if misAssemblyType == "Sw":
            relationChains.induceSwitchMisAssembly(newCtg,
                                                   orderIndex,
                                                   start,
                                                   end,
                                                   switchEffectWindowSize)
        elif misAssemblyType == "Err":
            relationChains.induceBaseErrorMisAssembly(newCtg,
                                                      orderIndex,
                                                      start,
                                                      end)
        elif misAssemblyType == "Dup":
            relationChains.induceDuplicationMisAssembly(newCtg,
                                                        orderIndex,
                                                        start,
                                                        end)
        elif misAssemblyType == "Col":
            relationChains.induceCollapseMisAssembly(newCtg,
                                                     orderIndex,
                                                     start,
                                                     end)






def getContigLengths(sequences):

    contigLengths = {}
    for name, seq in sequences.items():
        contigLengths[name] = len(seq)

    return contigLengths

def parseAnnotations(annotationsJson):
    # Opening JSON file
    f = open(annotationsJson)

    # returns JSON object as
    # a dictionary
    bedPaths = json.load(f)

    annotations = {}
    for name, path in bedPaths.items():
        annotations[name] = BlockList.parseBed(path)

    return annotations

def parseMisAssemblySizeTable(tsvPath):
    misAssemblySizeTable = pd.read_csv(tsvPath, "\t")
    misAssemblySizeTable.index = misAssemblySizeTable["length_kb"]
    for row in misAssemblySizeTable.index:
        for annotation in misAssemblySizeTable.columns[1:]:
            # convert string into an array of numbers
            # each number is the number of misassmblies that has to be
            # generated from each misassembly type
            # ["Switch error", "Erroneous", "Duplicated", "Collapsed"]
            misAssemblySizeTable.at[row, annotation] = np.array(re.split(' *, *', misAssemblySizeTable.at[row, annotation]),
                                                                 dtype=int)
    return misAssemblySizeTable


def getTotalLengthOfMisAssembliesForAllAnnotations(misAssemblySizeTable):
    """
    Returns a dictionary whose keys are annotation names. Each value is a list of total lengths of the
    requested misassemblies of different sizes.
    :param misAssemblySizeTable: A DataFrame that contains the numbers and the lengths of the requested misassemblies
    """

    totalLengths = defaultdict(list)
    print("Total lengths of requested mis-assemblies")
    for iRow in range(misAssemblySizeTable.shape[0]):
        misAssemblyLengthKb = int(misAssemblySizeTable.at[iRow, "length_kb"])
        print(f"\t Mis-assembly Size: {misAssemblyLengthKb}")
        for annotation in misAssemblySizeTable.columns[1:]:
            counts = misAssemblySizeTable.at[iRow, annotation]
            totLen = sum(counts) * misAssemblyLengthKb * 1e3
            print(f"\t\tAnnotation: {annotation} -> {totLen/1e3} Kb")
            totalLengths[annotation].append(totLen)
    return totalLengths


def checkFeasiblity(relationChains, annotations, misAssemblySizeTable, safetyFactor = 2):
    """
    Checks if there exist enough number of contiguous blocks per annotation to make misassemblies with the
    requested sizes
    Note that it does not guarantee that we won't be out of enough contiguous blocks for creating misassemblies
    since the locations of blocks are selected randomly it may happen that previously created misassemblies
    do not leave enough contiguous blocks for the next misassemblies. It depends on the randomly chosen
    locations of the misassemblies however, running out of blocks would be a rare event if the heuristic
    of this function is passed and it returns True. The higher the value of safetyFactor the less probable it would
    be to run out of blocks.
    :param relationChains:  Chains of relations created out of the alignments
    :param annotations: A list of annotation names
    :param misAssemblySizeTable: A DataFrame that contains the numbers and the sizes of the requested misassemblies
    :param safetyFactor: The total number of bases available in the blocks longer than each requested misassembly size
                         should be more than "safetyFactor" times the total lengths of the misassemblies longer than
                         the related requested size.
                         This check is performed separately per annotation.
    :return: True if feasible otherwise False
    """
    blockSizes = np.array(misAssemblySizeTable["length_kb"], dtype=int)
    totalAvailableLengthsPerAnnotation = relationChains.getTotalLengthOfLongerBlocksForAllAnnotations(annotations, blockSizes)
    totalRequestedLengthsPerAnnotation = getTotalLengthOfMisAssembliesForAllAnnotations(misAssemblySizeTable)

    for annotation in annotations:
        for requestedLen, availableLen in zip(totalRequestedLengthsPerAnnotation[annotation],
                                              totalAvailableLengthsPerAnnotation[annotation]):
            if availableLen < safetyFactor * requestedLen:
                return False
    return True



def main():
    parser = argparse.ArgumentParser(description='Given a high quality diploid assembly and the alignments of hap2 to hap1 induces misassemblies, returns the falsified assembly and the coordinates of the misassembled regions')
    parser.add_argument('--paf', type=str,
                        help='(PAF format) The alignments of hap2 contigs to hap1. It should include the cigar strings.')
    parser.add_argument('--hap1', type=str,
                        help='hap1 fasta file')
    parser.add_argument('--hap2', type=str,
                        help='hap2 fasta file')
    parser.add_argument('--misAssemblyTsv', type=str,
                        help='TSV file that contains the length and number of desired misassemblies per annotation,\
                              the first column is "length_kb" which contains the misassembly sizes in Kb. The rest of the columns \
                              should be named based on the annotation names in the json file. Each entry of the annotaion columns \
                              should be a list of 4 numbers separated by comma like "2,2,2,2"; numbers are for 4 types of misassemblies, \
                              "Switch error", "Erroneous", "Duplication", "Collapsed"')
    parser.add_argument('--annotations', type=str,
                        help='JSON file that contains the annotation names and the paths to the corresponding annotation bed files. The tracks in all the given bed files should not have any overlap.')
    parser.add_argument('--outputDir', type=str, default= "falsifier_outputs",
                        help='Output directory for saving the falsified assembly, the new annotation and mis-assembly coordinates')
    parser.add_argument('--minAlignmentLength', type=int, default = 500000,
                        help='Minimum length of the alignments to be used by the program [Default = 500000 (500Kb)]')
    parser.add_argument('--marginLength', type=int, default = 10000,
                        help='Length of the margins at the ends of each unbroken alignment (a contiguous homology relation) where no mis-assembly is permitted to be created [Default = 10000 (10 Kb)]')
    parser.add_argument('--switchWindowLength', type=int, default = 5000,
                        help='Each switch error is shown as one/multiple bed tracks around the switching point or misjoin. This parameter defines the length of the window labeled as "Msj" on each side of the misjoins created because of a haplotype switch (should be smaller than --marginLength) [Default = 5000 (5 Kb)]')
    parser.add_argument('--minOverlapRatio', type=float, default = 0.75,
                        help='Minimum overlap ratio that each misassembly block should have with the related annotation (should be greater than or equal to 0.5) [Default = 0.75]')
    parser.add_argument('--contigSuffix', type=str, default = "_f",
                        help='The suffix to be added to the names of the new contigs (which could be either intact or with mis-assemblies) [Default = "_f"]')
    parser.add_argument('--safetyFactor', type=float, default = 2.0, help='The factor for checking the feasibility of inducing the misassmblies of requested numbers and sizes [Default = 2.0]')



    # parse arguments
    args = parser.parse_args()
    pafPath = args.paf 
    hap1FastaPath = args.hap1
    hap2FastaPath = args.hap2
    outputDir = args.outputDir
    misAssemblyTsvPath = args.misAssemblyTsv
    annotationsJsonPath = args.annotations
    minAlignmentLength = args.minAlignmentLength
    switchWindowLength = args.switchWindowLength
    minOverlapRatio = args.minOverlapRatio
    marginLength = args.marginLength
    contigSuffix = args.contigSuffix
    safetyFactor = args.safetyFactor


    # parse the alignments from paf file
    alignments = []
    with open(pafPath,"r") as fPaf:
        for line in fPaf:
            alignment = Alignment(line)
            # Filter short alignments based the given minimum length
            alignmentLengthQuery = alignment.contigEnd - alignment.contigStart + 1
            alignmentLengthRef = alignment.chromEnd - alignment.chromStart + 1
            if alignmentLengthQuery >= minAlignmentLength and alignmentLengthRef >= minAlignmentLength:
                alignments.append(alignment)

    # parse hap1 sequences
    hap1Sequences = parseFasta(hap1FastaPath)
    hap1ContigNames = list(hap1Sequences.keys())

    # extend hap1 sequences to contain hap2 sequences too
    diploidSequences = hap1Sequences
    diploidSequences.update(parseFasta(hap2FastaPath))

    diploidContigLengths = getContigLengths(diploidSequences)

    annotationBlockLists = parseAnnotations(annotationsJsonPath)
    annotationNames = list(annotationBlockLists.keys())

    # parse misassembly length/number tsv
    misAssemblySizeTable = parseMisAssemblySizeTable(misAssemblyTsvPath)
    for annotation in misAssemblySizeTable.columns[1:]:
        if annotation not in annotationNames:
            print(f"Error: {annotation} exists in the mis-assembly tsv but its path is not given in the json file!")
            exit()

    # find the hap1 coordinates with exactly one alignment from hap2
    hap1UniqueBlocksPerContig = getBlockListsWithSingleAlignmentPerRefContig(alignments)
    hap1UniqueAlignments = subsetAlignmentsToRefBlocks(alignments, hap1UniqueBlocksPerContig)

    # find the hap2 coordinates with exactly one alignment to hap1 among the
    # alignments which are uniquely mapped to hap1
    hap2UniqueBlocksPerContig = getBlockListsWithSingleAlignmentPerQueryContig(hap1UniqueAlignments)

    # get the alignments that shows 1-to-1 correspondence between hap1 and hap2
    uniqueAlignments = subsetAlignmentsToQueryBlocks(hap1UniqueAlignments, hap2UniqueBlocksPerContig)


    # make relation chains for the whole diploid assembly
    relationChains = HomologyRelationChains(uniqueAlignments,
                                            diploidContigLengths,
                                            hap1ContigNames,
                                            contigSuffix)


    # fill the annotation attributes of the blocks based on the annotations of the assembly
    relationChains.fillAnnotationBlockListsFromOriginalContigs(annotationBlockLists,
                                                               diploidContigLengths,
                                                               contigSuffix)

    if checkFeasiblity(relationChains, annotationNames, misAssemblySizeTable, safetyFactor = safetyFactor):
        print(f"Feasibilty is PASSED with the safety factor of {safetyFactor}")
    else:
        print(f"Feasibilty is NOT PASSED with the safety factor of {safetyFactor}")
        exit()



    misAssemblyTypes = ["Sw", "Err", "Dup", "Col"]
    misAssemblySizesSortedKb = sorted(np.array(misAssemblySizeTable["length_kb"]), reverse=True)
    for misAssemblySizeKb in misAssemblySizesSortedKb:
        # for each mis-assembly size, the start locations for sampling should be updated
        relationChains.updateAnnotationStartLocationsForSampling(annotationNames,
                                                                 misAssemblySizeKb * 1000,
                                                                 minOverlapRatio,
                                                                 marginLength)
        # induce all mis-assemblies of the same size across all annotations
        for annotation in annotationNames:
            misAssemblyCounts = misAssemblySizeTable.at[misAssemblySizeKb, annotation]
            for misAssemblyType, misAssemblyCount in  zip(misAssemblyTypes, misAssemblyCounts) :
                induceMultipleMisAssembliesOfTheSameType(relationChains,
                                                         annotationNames,
                                                         misAssemblyType,
                                                         misAssemblySizeKb * 1000,
                                                         misAssemblyCount,
                                                         switchWindowLength)
                print(f"Mis-assemblies are created: annotation = {annotation} type = {misAssemblyType}, length = {misAssemblySizeKb} Kb, count = {misAssemblyCount}")

    print("Creating mis-assemblies is Done!")

    print("Writing Fasta file for the falsified assembly")
    fastaPath = os.path.join(outputDir,"falsified_asm.fasta")
    relationChains.writeNewContigsToFasta(diploidSequences, fastaPath)
    print(f"The falsified assembly is written to {fastaPath}")

    for annotation in annotationNames:
        bedPath = os.path.join(outputDir, "falsified_asm.{annotation}.bed")
        relationChains.writeAnnotationCoordinatesToBed(annotation, bedPath)
        print(f"The coordinates of the annotation, {annotation} are written to {bedPath}")

    bedPath = os.path.join(outputDir, "falsified_asm.mis_assemblies.bed")
    relationChains.writeMisAssemblyCoordinatesToBed(bedPath)
    print(f"The coordinates of misassemblies are written to {bedPath}")















if __name__ == "__main__": main()