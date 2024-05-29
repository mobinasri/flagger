import sys
import argparse
from block_utils import Alignment
from alignment_utils import *
from collections import defaultdict
from copy import deepcopy
from falsifier_utils import *
import numpy as np
import pandas as pd
import json
import re
import os
import datetime

def induceMultipleMisjoins(relationChains, annotation, misAssemblyCount, misjoinEffectWindowSize):
    misAssemblyType = "Misjoin"
    for i in range(misAssemblyCount):
        iter = 0
        res = False
        while iter < 10:
            if 0 < iter:
                print(f"[{datetime.datetime.now()}] Retrying inducing the failed mis-assembly ({iter+1}/10 retries)")


            newCtg_1, orderIndex_1, loc_1 = relationChains.getRandomMisjoinLocation(annotation, newCtgListToExclude=[])
            if newCtg_1 == None:
                print(f"[{datetime.datetime.now()}] No more blocks remaining for sampling misjoin (Skipping {misAssemblyCount-i}/{misAssemblyCount}) : annotation = {annotation}, misAssemblyType = {misAssemblyType}")
                return i
            newCtg_2, orderIndex_2, loc_2 = relationChains.getRandomMisjoinLocation(annotation, newCtgListToExclude=[newCtg_1])
            if newCtg_2 == None:
                print(f"[{datetime.datetime.now()}] No more blocks remaining for sampling misjoin (Skipping {misAssemblyCount-i}/{misAssemblyCount}) : annotation = {annotation}, misAssemblyType = {misAssemblyType}")
                return i
            print(f"[{datetime.datetime.now()}] Mis-assembly location (misjoin-pair-1):\t {newCtg_1}[{orderIndex_1}]:\t{loc_1}\tlen={relationChains.relationChains[newCtg_1][orderIndex_1].block.origEnd - relationChains.relationChains[newCtg_1][orderIndex_1].block.origStart + 1}")
            print(f"[{datetime.datetime.now()}] Mis-assembly location (misjoin-pair-2):\t {newCtg_2}[{orderIndex_2}]:\t{loc_2}\tlen={relationChains.relationChains[newCtg_2][orderIndex_2].block.origEnd - relationChains.relationChains[newCtg_2][orderIndex_2].block.origStart + 1}")

            res = relationChains.induceMisjoinMisAssembly(newCtg_1,
                                                          orderIndex_1,
                                                          loc_1,
                                                          newCtg_2,
                                                          orderIndex_2,
                                                          loc_2,
                                                          misjoinEffectWindowSize)

            if res == True:
                break
            else:
                print(f"[{datetime.datetime.now()}] Warning: Mis-assembly could not be induced (because of projection).")
                iter += 1
        # if the mis-assembly could not be created after 10 times trying 
        if res == False:
            return i

    return misAssemblyCount

def induceMultipleMisAssembliesOfTheSameType(relationChains, annotation, misAssemblyType, misAssemblySize, misAssemblyCount, switchEffectWindowSize):
    assert(misAssemblyType in ["Sw", "Err", "Dup", "Col"])
    for i in range(misAssemblyCount):
        iter = 0
        res = False
        while iter < 10:
            if 0 < iter:
                print(f"[{datetime.datetime.now()}] Retrying inducing the failed mis-assembly ({iter+1}/10 retries)")
            newCtg, orderIndex, start, end = relationChains.getRandomMisAssemblyInterval(annotation, misAssemblySize)
            if newCtg == None:
                print(f"[{datetime.datetime.now()}] No more blocks remaining for sampling (Skipping {misAssemblyCount-i}/{misAssemblyCount}) : annotation = {annotation}, misAssemblyType = {misAssemblyType}, misAssemblySize = {misAssemblySize}")
                return i
            print(f"[{datetime.datetime.now()}] Mis-assembly location:\t {newCtg}[{orderIndex}]:\t{start}\t{end}\tlen={relationChains.relationChains[newCtg][orderIndex].block.origEnd - relationChains.relationChains[newCtg][orderIndex].block.origStart + 1}")

            if misAssemblyType == "Sw":
                res = relationChains.induceSwitchMisAssembly(newCtg,
                                                             orderIndex,
                                                             start,
                                                             end,
                                                             switchEffectWindowSize)
            elif misAssemblyType == "Err":
                res = relationChains.induceBaseErrorMisAssembly(newCtg,
                                                                orderIndex,
                                                                start,
                                                                end)
            elif misAssemblyType == "Dup":
                res = relationChains.induceDuplicationMisAssembly(newCtg,
                                                                  orderIndex,
                                                                  start,
                                                                  end)
            elif misAssemblyType == "Col":
                res = relationChains.induceCollapseMisAssembly(newCtg,
                                                               orderIndex,
                                                               start,
                                                               end)
            if res == True:
                break
            else:
                print(f"[{datetime.datetime.now()}] Warning: Mis-assembly could not be induced (because of projection).")
                iter += 1
        # if the mis-assembly could not be created after 10 times trying 
        if res == False:
            return i

    return misAssemblyCount





def getContigLengths(sequences):

    contigLengths = {}
    for name, seq in sequences.items():
        contigLengths[name] = len(seq)

    return contigLengths

def parseNumberOfMisjoins(misjoinJson):
    # Opening JSON file
    f = open(misjoinJson)

    # returns JSON object as
    # a dictionary
    misjoinCounts = json.load(f)


    return misjoinCounts

def parseAnnotationsPerContig(annotationsJson):
    # Opening JSON file
    f = open(annotationsJson)

    # returns JSON object as
    # a dictionary
    bedPaths = json.load(f)

    annotationsPerContig = defaultdict(dict)
    annotationNames = []
    for name, path in bedPaths.items():
        annotationNames.append(name)
        blocksForOneAnnotation = BlockList.parseBed(path, saveFourthColumnAsNumeric=False)
        for ctg, blockList in blocksForOneAnnotation.items():
            annotationsPerContig[ctg][name] = blockList

    return annotationsPerContig, annotationNames

def parseMisAssemblySizeTable(tsvPath):
    misAssemblySizeTable = pd.read_csv(tsvPath, sep="\t")
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
    print(f"[{datetime.datetime.now()}] Total lengths of requested mis-assemblies:")
    for row in misAssemblySizeTable.index:
        misAssemblyLengthKb = int(misAssemblySizeTable.at[row, "length_kb"])
        print(f"\t Mis-assembly Size: {misAssemblyLengthKb}Kb")
        for annotation in misAssemblySizeTable.columns[1:]:
            counts = misAssemblySizeTable.at[row, annotation]
            totLen = sum(counts) * misAssemblyLengthKb * 1e3
            print(f"\t\tAnnotation: {annotation} -> {totLen/1e3} Kb")
            totalLengths[annotation].append(totLen)
    return totalLengths


def getTotalCountOfMisAssembliesForAllAnnotations(misAssemblySizeTable):
    """
    Returns a dictionary whose keys are annotation names. Each value is a list of total counts of the
    requested misassemblies of different sizes.
    :param misAssemblySizeTable: A DataFrame that contains the numbers and the lengths of the requested misassemblies
    """

    totalCounts = defaultdict(list)
    print(f"[{datetime.datetime.now()}] Total counts of requested mis-assemblies:")
    for row in misAssemblySizeTable.index:
        misAssemblyLengthKb = int(misAssemblySizeTable.at[row, "length_kb"])
        print(f"\t Mis-assembly Size: {misAssemblyLengthKb}Kb")
        for annotation in misAssemblySizeTable.columns[1:]:
            counts = misAssemblySizeTable.at[row, annotation]
            totCount = sum(counts)
            print(f"\t\tAnnotation: {annotation} -> {totCount}")
            totalCounts[annotation].append(totCount)
    return totalCounts

def checkFeasiblity(relationChains, misAssemblySizeTable):
    """
    Checks if there exist enough number of contiguous blocks per annotation to make misassemblies with the
    requested sizes
    Note that it if this function considers only the worst case scenario if it returns True it is guaranteed that
    all requested misassemblies can be created without running out of annotation blocks. However, since the locations
    of misassemblies are random if it returns False it might happen that all misassemblies can be created. It depends on
    how close the requested number and length of misassembly blocks are to the available blocks.
    :param relationChains:  Chains of relations created out of the alignments
    :param misAssemblySizeTable: A DataFrame that contains the numbers and the sizes of the requested misassemblies
    :return: True if feasible otherwise False
    """

    print(f"[{datetime.datetime.now()}] Feasibility logs:")
    print(f"#annotation\tmisassembly_block_size_kb,\tlower_bound_count\trequested_total_count\tstatus")
    annotations = misAssemblySizeTable.columns[1:]
    flag = True
    for annotation in annotations:
        for row in misAssemblySizeTable.index:
            misAssemblyLengthKb = int(misAssemblySizeTable.at[row, "length_kb"])
            requestedCount = sum(misAssemblySizeTable.at[row, annotation])
            lowerBoundOnCount = relationChains.getLowerBoundOnNumberOfMisassemblies(annotation, misAssemblyLengthKb * 1e3)
            if requestedCount <= lowerBoundOnCount:
                status = "PASSED"
            else:
                status = "NOT_PASSED"
                flag = False
            print(f"{annotation}\t{misAssemblyLengthKb}\t{lowerBoundOnCount}\t{requestedCount}\t{status}")
    return flag


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
                              should be a list of 4 numbers separated by comma like "1,2,2,2"; numbers are for 4 types of misassemblies, \
                              "Switch error", "Erroneous", "Duplication", "Collapsed"')
    parser.add_argument('--misjoinJson', type=str,
                        help='JSON file that contains the number of misjoins (actually pairs of misjoins) per annotation,\
                              Each key is an annotation name and each value is the number of pairs of misjoins that \
                              should be created in the corresponding annotation')
    parser.add_argument('--annotationsJson', type=str,
                        help='JSON file that contains the annotation names and the paths to the corresponding annotation bed files. The tracks in all the given bed files should not have any overlap.')
    parser.add_argument('--outputDir', type=str, default= "falsifier_outputs",
                        help='Output directory for saving the falsified assembly, the new annotation and mis-assembly coordinates')
    parser.add_argument('--minAlignmentLength', type=int, default = 50000,
                        help='Minimum length of the alignments to be used by the program [Default = 50000 (50Kb)]')
    parser.add_argument('--marginLength', type=int, default = 2000,
                        help='Length of the margins at the ends of each unbroken alignment (a contiguous homology relation) where no mis-assembly is permitted to be created [Default = 1000 (1 Kb)]')
    parser.add_argument('--switchFlagWindowLength', type=int, default = 1000,
                        help='Each switch error is shown as one/multiple bed tracks around the switching point. This parameter defines the length of the window labeled as "Sw" on each side of switching points created because of haplotype switches (should be smaller than --marginLength) [Default = 5000 (5 Kb)]')
    parser.add_argument('--misjoinFlagWindowLength', type=int, default = 5000,
                        help='Each misjoin is shown as one/multiple bed tracks around the misjoin. This parameter defines the length of the window labeled as "Msj" on each side of a misjoin between non-homologous contigs (should be smaller than --marginLength) [Default = 5000 (5 Kb)]')
    parser.add_argument('--minOverlapRatio', type=float, default = 0.75,
                        help='Minimum overlap ratio that each misassembly block should have with the related annotation (should be greater than or equal to 0.5) [Default = 0.75]')
    parser.add_argument('--contigSuffix', type=str, default = "_f",
                        help='The suffix to be added to the names of the new contigs (which could be either intact or with mis-assemblies) [Default = "_f"]')
    parser.add_argument('--singleBaseErrorRate', type=float, default = 0.05, help='The rate of single-base errors that will be induced in "Err" blocks [Default = 0.05]')
    parser.add_argument('--maxGapLength', type=int, default = 500, help='Split alignments into smaller alignments with no gaps longer than this parameter[Default = 500]')



    # parse arguments
    args = parser.parse_args()
    pafPath = args.paf 
    hap1FastaPath = args.hap1
    hap2FastaPath = args.hap2
    outputDir = args.outputDir
    misAssemblyTsvPath = args.misAssemblyTsv
    misjoinJsonPath = args.misjoinJson
    annotationsJsonPath = args.annotationsJson
    minAlignmentLength = args.minAlignmentLength
    switchFlagWindowLength = args.switchFlagWindowLength
    misjoinFlagWindowLength = args.misjoinFlagWindowLength
    minOverlapRatio = args.minOverlapRatio
    marginLength = args.marginLength
    contigSuffix = args.contigSuffix
    maxGapLength = args.maxGapLength
    singleBaseErrorRate = args.singleBaseErrorRate


    # parse the alignments from paf file
    alignments = []
    with open(pafPath,"r") as fPaf:
        for line in fPaf:
            alignment = Alignment(line)
            # Filter short alignments based the given minimum length
            alignmentLengthQuery = alignment.contigEnd - alignment.contigStart + 1
            alignmentLengthRef = alignment.chromEnd - alignment.chromStart + 1
            if alignmentLengthQuery >= minAlignmentLength and \
                    alignmentLengthRef >= minAlignmentLength and \
                    alignment.isPrimary:
                # split alignments into alignments with short gaps
                for alignmentWithShortGap in splitIntoAlignmentsWithShortGaps(alignment, maxGapLength):
                    alignments.append(alignmentWithShortGap)

    # parse hap1 sequences
    hap1Sequences = parseFasta(hap1FastaPath)
    hap1ContigNames = list(hap1Sequences.keys())

    # extend hap1 sequences to contain hap2 sequences too
    diploidSequences = hap1Sequences.copy()
    diploidSequences.update(parseFasta(hap2FastaPath))

    diploidContigLengths = getContigLengths(diploidSequences)

    print(f"[{datetime.datetime.now()}] Prased contigs and their sizes")
    for contigName, contigSize in diploidContigLengths.items():
        print(contigName, contigSize)
    totalGenomeSizeKb = sum(diploidContigLengths.values()) / 1e3
    print(f"[{datetime.datetime.now()}] Diploid genome size is {totalGenomeSizeKb} Kb")

    annotationBlockLists, annotationNames = parseAnnotationsPerContig(annotationsJsonPath)

    # parse misassembly length/number tsv
    misAssemblySizeTable = parseMisAssemblySizeTable(misAssemblyTsvPath)
    for annotation in misAssemblySizeTable.columns[1:]:
        if annotation not in annotationNames:
            print(f"[{datetime.datetime.now()}] Error: {annotation} exists in the mis-assembly tsv but its path is not given in the json file!")
            exit()
    annotationsForCreatingMisAssembly = list(misAssemblySizeTable.columns[1:])

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

    # get the total length of all annotation blocks  (both ref and query) without considering the 1-to-1 mapping
    # this will be used for calculating the actual misassembly rate per annotation in the final falsified assembly
    # note that blockSize is set to [0] which means do not filter annotation blocks based on their size and get
    # all of them
    annotationLengths = relationChains.getTotalLengthOfLongerBlocksForAllAnnotations(annotationNames,
                                                                                    None,
                                                                                    onlyRefInHomology=False)
    annotationLengths["whole_genome"] = totalGenomeSizeKb * 1e3

    if checkFeasiblity(relationChains, misAssemblySizeTable):
        print(f"[{datetime.datetime.now()}] Feasibilty is PASSED.")
    else:
        print(f"[{datetime.datetime.now()}] Feasibilty is NOT PASSED. Please use more contiguous alignments/annotations tracks or request for shorter (or fewer) misassemblies.")
        exit()


    # update sampling locations for misjoins
    # misassembly length is set to zero since it does not matter
    relationChains.updateAnnotationBlocksForSampling(annotationNames,
                                                     0,
                                                     minOverlapRatio,
                                                     marginLength)
    print(f"[{datetime.datetime.now()}] Creating misjoins")
    numberOfMisjoinsPerAnnotation = parseNumberOfMisjoins(misjoinJsonPath)

    for annotation in annotationNames:
        if annotation in numberOfMisjoinsPerAnnotation:
            misAssemblyCount = numberOfMisjoinsPerAnnotation[annotation]
            res = induceMultipleMisjoins(relationChains,
                                         annotation,
                                         misAssemblyCount,
                                         misjoinFlagWindowLength)
            print(f"[{datetime.datetime.now()}] ({res}/{misAssemblyCount}) Mis-assemblies could be created successfully:\tannotation = {annotation},\ttype = misjoin")



    misAssemblyTypes = ["Sw", "Err", "Dup", "Col"]
    # first create large misassemblies since finding intervals for smaller misassemblies is easier
    misAssemblySizesSortedKb = sorted(np.array(misAssemblySizeTable["length_kb"]), reverse=True)
    total_successful = 0
    total_requested = 0
    totalMisAssembledBasesKbByAnnotationAndType = {annotation: {"Sw":0, "Err":0, "Dup":0, "Col":0}
                                                   for annotation in annotationsForCreatingMisAssembly}
    totalMisAssembledBasesKbByAnnotationAndType["whole_genome"] = {"Sw":0, "Err":0, "Dup":0, "Col":0}
    for misAssemblySizeKb in misAssemblySizesSortedKb:
        # for each mis-assembly size, the start locations for sampling should be updated
        relationChains.updateAnnotationBlocksForSampling(annotationNames,
                                                                 misAssemblySizeKb * 1000,
                                                                 minOverlapRatio,
                                                                 marginLength)
        # induce all mis-assemblies of the same size across all annotations
        for annotation in annotationsForCreatingMisAssembly:
            misAssemblyCounts = misAssemblySizeTable.at[misAssemblySizeKb, annotation]
            for misAssemblyType, misAssemblyCount in  zip(misAssemblyTypes, misAssemblyCounts) :
                #print(relationChains.newCtgAnnotationWeightsForSampling[annotation])
                res = induceMultipleMisAssembliesOfTheSameType(relationChains,
                                                               annotation,
                                                               misAssemblyType,
                                                               misAssemblySizeKb * 1000,
                                                               misAssemblyCount,
                                                               switchFlagWindowLength)
                print(f"[{datetime.datetime.now()}] ({res}/{misAssemblyCount}) Mis-assemblies could be created successfully:\tannotation = {annotation},\ttype = {misAssemblyType},\tlength = {misAssemblySizeKb} Kb")
                total_successful += res
                total_requested += misAssemblyCount
                if misAssemblyType == "Sw":
                    totalMisAssembledBasesKbByAnnotationAndType["whole_genome"][misAssemblyType] += 2 * res * misAssemblySizeKb
                    totalMisAssembledBasesKbByAnnotationAndType[annotation][misAssemblyType] += 2 * res * misAssemblySizeKb
                else:
                    totalMisAssembledBasesKbByAnnotationAndType["whole_genome"][misAssemblyType] += 1 * res * misAssemblySizeKb
                    totalMisAssembledBasesKbByAnnotationAndType[annotation][misAssemblyType] += 1 * res * misAssemblySizeKb

    print(f"[{datetime.datetime.now()}] Creating mis-assemblies is Done! ({total_successful}/{total_requested}) mi-assemblies could be created successfully.")

    print(f"[{datetime.datetime.now()}] Writing actual misassembly rates whole genome and per annotation in the final falsified assembly.")
    misAssemblyRateTextPath = os.path.join(outputDir, f"falsified_asm.misassembly_rates.txt")
    with open(misAssemblyRateTextPath, "w") as f:
        f.write(f"#annotation\tmisassembly_type\ttotal_misassembly_size_kb\ttotal_annotation_size_kb\tmisassembly_rate_percent\n")
        # Print mis-assembly rate per annotation (both per misassembly type and altogether)
        for annotation in annotationsForCreatingMisAssembly + ["whole_genome"]:
            totalMisAssembledBasesKb = 0
            totalAnnotationLengthKb = annotationLengths[annotation] / 1e3
            for misAssemblyType, misAssemblySizeKb in totalMisAssembledBasesKbByAnnotationAndType[annotation].items():
                misAssemblyRate = misAssemblySizeKb / totalAnnotationLengthKb
                f.write(f"{annotation}\t{misAssemblyType}\t{misAssemblySizeKb}\t{totalAnnotationLengthKb}\t{misAssemblyRate * 100:0.3f}\n")
                totalMisAssembledBasesKb += misAssemblySizeKb
            misAssemblyRateTotal = totalMisAssembledBasesKb / totalAnnotationLengthKb
            f.write(f"{annotation}\tALL\t{totalMisAssembledBasesKb}\t{totalAnnotationLengthKb}\t{misAssemblyRateTotal * 100:0.3f}\n")


    os.makedirs(outputDir, exist_ok = True)
    print(f"[{datetime.datetime.now()}] Writing Fasta file for the falsified assembly")
    diploidFastaPath = os.path.join(outputDir,"falsified_asm.dip.fa")
    hap1FastaPath = os.path.join(outputDir,"falsified_asm.hap1.fa")
    hap2FastaPath = os.path.join(outputDir,"falsified_asm.hap2.fa")
    relationChains.writeNewContigsToFasta(diploidSequences, diploidFastaPath, hap1FastaPath, hap2FastaPath, singleBaseErrorRate)
    print(f"[{datetime.datetime.now()}] The falsified assembly is written to {diploidFastaPath} (hap1={hap1FastaPath}, hap2={hap2FastaPath})")

    for annotation in annotationNames:
        bedPath = os.path.join(outputDir, f"falsified_asm.{annotation}.bed")
        relationChains.writeAnnotationCoordinatesToBed(annotation, bedPath)
        print(f"[{datetime.datetime.now()}] The coordinates of the annotation, {annotation} are written to {bedPath}")

    bedPath = os.path.join(outputDir, "falsified_asm.mis_assemblies.bed")
    relationChains.writeMisAssemblyCoordinatesToBed(bedPath)
    print(f"[{datetime.datetime.now()}] The coordinates of misassemblies are written to {bedPath}")



if __name__ == "__main__": main()
