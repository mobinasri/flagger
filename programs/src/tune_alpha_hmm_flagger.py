import sys
import argparse
from collections import defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
import json
import re
import os
import subprocess
import shlex
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
from smt.applications import EGO
from smt.surrogate_models import KRG
from smt.utils.design_space import DesignSpace

DIMENSION = 10

def getRandomX(lowerBound, upperBound):
    dim = DIMENSION
    x = (upperBound - lowerBound) * np.random.random_sample(dim) + lowerBound
    return x

def getStartPoints(lowerBound, upperBound, numberOfPoints):
    x = np.array([getRandomX(lowerBound, upperBound) for i in range(numberOfPoints)])
    return x

def getXLimits(lowerBound, upperBound):
    return np.array([[lowerBound, upperBound] for i in range(DIMENSION)])

def convertAlphaMatrixToX(alphaMatrix):
    x = []
    x.append(alphaMatrix[0,0])
    x.append(alphaMatrix[0,2])
    x.append(alphaMatrix[1,1])
    x.append(alphaMatrix[1,2])
    x.append(alphaMatrix[2,0])
    x.append(alphaMatrix[2,1])
    x.append(alphaMatrix[2,2])
    x.append(alphaMatrix[2,3])
    x.append(alphaMatrix[3,2])
    x.append(alphaMatrix[3,3])
    return np.array(x)

def convertXToAlphaMatrix(x):
    alphaMatrix = np.zeros((4, 4))
    alphaMatrix[0,0] = x[0]
    alphaMatrix[0,2] = x[1]
    alphaMatrix[1,1] = x[2]
    alphaMatrix[1,2] = x[3]
    alphaMatrix[2,0] = x[4]
    alphaMatrix[2,1] = x[5]
    alphaMatrix[2,2] = x[6]
    alphaMatrix[2,3] = x[7]
    alphaMatrix[3,2] = x[8]
    alphaMatrix[3,3] = x[9]
    return alphaMatrix

def getInputPrefix(inputPath):
    elems = inputPath.split('.')
    if elems[-1] == 'gz' and 2 < len(elems):
        return '.'.join(elems[:-2])
    else:
        return '.'.join(elems[:-1])

def saveAlphaMatrixInTSV(alphaMatrix, tsvPath):
    np.savetxt(tsvPath, alphaMatrix, delimiter='\t', fmt="%.3f")

def getParametersString(alphaTsvPath, inputPath, outputDir, modelType, otherParamsString):
    return f"--alpha {alphaTsvPath} --input {inputPath} --outputDir {outputDir} --modelType {modelType} {otherParamsString}"

def getOverlapBasedScore(outputDir, annotationLabel, sizeLabel):
    tsvPath = [os.path.join(outputDir, a) for a in os.listdir(outputDir) if "final.benchmarking.tsv" in a][0]
    table = pd.read_csv(tsvPath,sep="\t")
    table = table.rename(columns={"#Metric_Type": "Metric_Type"})
    table = table.loc[lambda table: table.Metric_Type == "overlap_based"]
    table = table.loc[lambda table: table.Category_Name == annotationLabel]
    table = table.loc[lambda table: table.Size_Bin_Name == sizeLabel]
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN_NO_HAP"]
    return float(table["F1-Score"].item())

def getBaseLevelScore(outputDir, annotationLabel, sizeLabel):
    tsvPath = [os.path.join(outputDir, a) for a in os.listdir(outputDir) if "final.benchmarking.tsv" in a][0]
    table = pd.read_csv(tsvPath,sep="\t")
    table = table.rename(columns={"#Metric_Type": "Metric_Type"})
    table = table.loc[lambda table: table.Metric_Type == "base_level"]
    table = table.loc[lambda table: table.Category_Name == annotationLabel]
    table = table.loc[lambda table: table.Size_Bin_Name == sizeLabel]
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN_NO_HAP"]
    return float(table["F1-Score"].item())

def getContiguityScore(outputDir, annotationLabel, sizeLabel):
    tsvPath = [os.path.join(outputDir, a) for a in os.listdir(outputDir) if "final.benchmarking.auN_ratio.tsv" in a][0]
    table = pd.read_csv(tsvPath,sep="\t")
    table = table.loc[lambda table: table.Category_Name == annotationLabel]
    table = table.loc[lambda table: table.Size_Bin_Name == sizeLabel]
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN"]
    return 100 * float(table["auN_Ratio"].item())


def functionToMinimize(x):
    y = []
    for i in range(len(x)):
        y.append([functionToMinimizeInternal(x[i])])
    return np.array(y)

# run hmm_flagger with the given alpha matrix on all input files
# and return the average of all combined scores
def functionToMinimizeInternal(x):
    functionToMinimizeInternal.pointIndex += 1
    outputDir = functionToMinimizeInternal.paramsDict['outputDir']
    inputPathTrainList = functionToMinimizeInternal.paramsDict['inputPathTrainList']
    inputPathTestList = functionToMinimizeInternal.paramsDict['inputPathTestList']
    modelType = functionToMinimizeInternal.paramsDict['modelType']
    otherParamsString = functionToMinimizeInternal.paramsDict['otherParamsString']
    maxJobs = functionToMinimizeInternal.paramsDict['maxJobs']
    annotationLabel = functionToMinimizeInternal.paramsDict['annotationLabel']
    sizeLabel = functionToMinimizeInternal.paramsDict['sizeLabel']
    pointIndex = functionToMinimizeInternal.pointIndex

    alphaMatrix = convertXToAlphaMatrix(x)
    os.makedirs(f"{outputDir}/optimization_point_{pointIndex}", exist_ok=True)
    alphaTsvPath = f"{outputDir}/optimization_point_{pointIndex}/alpha_mat_point_{pointIndex}.tsv"
    saveAlphaMatrixInTSV(alphaMatrix, alphaTsvPath)

    internalOutputDirTrainList = []
    paramsStringList = []
    for inputPath in inputPathTrainList:
        prefix = getInputPrefix(inputPath)
        internalOutputDir = f"{outputDir}/optimization_point_{pointIndex}/{prefix}"
        os.makedirs(internalOutputDir, exist_ok=True)
        internalOutputDirTrainList.append(internalOutputDir)
        paramsString = getParametersString(alphaTsvPath, inputPath, internalOutputDir, modelType, otherParamsString)
        paramsStringList.append(paramsString)

    internalOutputDirTestList = []
    for inputPath in inputPathTestList:
        prefix = getInputPrefix(inputPath)
        internalOutputDir = f"{outputDir}/optimization_point_{pointIndex}/{prefix}"
        os.makedirs(internalOutputDir, exist_ok=True)
        internalOutputDirTestList.append(internalOutputDir)
        paramsString = getParametersString(alphaTsvPath, inputPath, internalOutputDir, modelType, otherParamsString)
        paramsStringList.append(paramsString)

    i = functionToMinimizeInternal.pointIndex
    n = functionToMinimizeInternal.numberOfStartPoints
    nIter = functionToMinimizeInternal.numberOfIterations
    isStartPoint = False
    if i <= n:
        print(f"[{datetime.datetime.now()}] Point Index (for EGO) = {i} (Start point {i} out of {n}) ", file=sys.stderr)
        isStartPoint = True
    else:
        print(f"[{datetime.datetime.now()}] Point Index (for EGO) = {i} (Iteration {i - n} out of {nIter})", file=sys.stderr)

    print(f"[{datetime.datetime.now()}] Initiating {len(paramsStringList)} hmm_flagger jobs for training/testing", file=sys.stderr)
    runHMMFlaggerForList(paramsStringList, maxJobs)

    combinedScoreTrainList = []
    score1TrainList = []
    score2TrainList = []
    score3TrainList = []
    for inputIndex, internalOutputDir in enumerate(internalOutputDirTrainList):
        score1 = getOverlapBasedScore(internalOutputDir, annotationLabel, sizeLabel)
        score2 = getBaseLevelScore(internalOutputDir, annotationLabel, sizeLabel)
        score3 = getContiguityScore(internalOutputDir, annotationLabel, sizeLabel)
        score1TrainList.append(score1)
        score2TrainList.append(score2)
        score3TrainList.append(score3)
        combined = (score1 + score2 + score3) / 3.0
        print(f"[{datetime.datetime.now()}] [train input index = {inputIndex}] OverlapBased = {score1:.3f}, BaseLevel = {score2:.3f}, Contiguity = {score3:.3f}, Combined = {combined:.3f}", file=sys.stderr)
        combinedScoreTrainList.append(combined)

    combinedScoreTestList = []
    score1TestList = []
    score2TestList = []
    score3TestList = []
    for inputIndex, internalOutputDir in enumerate(internalOutputDirTestList):
        score1 = getOverlapBasedScore(internalOutputDir, annotationLabel, sizeLabel)
        score2 = getBaseLevelScore(internalOutputDir, annotationLabel, sizeLabel)
        score3 = getContiguityScore(internalOutputDir, annotationLabel, sizeLabel)
        score1TestList.append(score1)
        score2TestList.append(score2)
        score3TestList.append(score3)
        combined = (score1 + score2 + score3) / 3.0
        print(f"[{datetime.datetime.now()}] [test input index = {inputIndex}] OverlapBased = {score1:.3f}, BaseLevel = {score2:.3f}, Contiguity = {score3:.3f}, Combined = {combined:.3f}", file=sys.stderr)
        combinedScoreTestList.append(combined)

    finalScoreTrain = sum(combinedScoreTrainList) / len(combinedScoreTrainList)
    if 0 < len(combinedScoreTestList):
        finalScoreTest = sum(combinedScoreTestList) / len(combinedScoreTestList)
    else:
        finalScoreTest = 0.0

    print(f"[{datetime.datetime.now()}] Alpha = ", file=sys.stderr)
    np.savetxt(sys.stderr, alphaMatrix, delimiter='\t', fmt="%.3f")
    print(f"Final Score Train = {finalScoreTrain:.3f}", file=sys.stderr)
    if 0 < len(combinedScoreTestList):
        print(f"Final Score Test = {finalScoreTest:.3f}", file=sys.stderr)
    sys.stderr.flush()

    functionToMinimizeInternal.allScoresTrain.append(["start" if isStartPoint else "iteration", finalScoreTrain, score1TrainList, score2TrainList, score3TrainList, combinedScoreTrainList])

    if 0 < len(combinedScoreTestList):
        functionToMinimizeInternal.allScoresTest.append(["start" if isStartPoint else "iteration", finalScoreTest, score1TestList, score2TestList, score3TestList, combinedScoreTestList])

    # times -1 since EGO algorithm minimizes the objective function
    return -1 * finalScoreTrain


def runHMMFlaggerForList(paramsStringList, maxJobs):
    success = 0
    with ThreadPoolExecutor(max_workers=maxJobs) as executor:
        futures = [executor.submit(runHMMFlagger, paramsString) for paramsString in paramsStringList]
        for future in as_completed(futures):
            # retrieve the result
            result = future.result()
            ##print(result[0])
            if result == 0:
                success += 1
    #print(f"[{datetime.datetime.now()}] Successfully ran {len(paramsStringList)} jobs for hmm_flagger", file=sys.stderr)
    sys.stderr.flush()

def runHMMFlagger(paramsString):
    """
    :param paramsString: hmm_flagger parameters
    """

    cmdString = f"hmm_flagger_new {paramsString}"
    x = subprocess.run(shlex.split(cmdString), capture_output=True)

    return x.returncode



def main():
    parser = argparse.ArgumentParser(
        description='This program takes a list of coverage/bin files and tune alpha matrix for hmm_flagger using Efficient Global Optimization (EGO). More info about EGO https://smt.readthedocs.io/en/latest/_src_docs/applications/ego.html')
    parser.add_argument('--inputFilesTrain', type=str,
                        help='Comma separated list of coverage/bin files that will be used for tuning/training alpha matrix')
    parser.add_argument('--inputFilesTest', type=str,
                        help='(Optional) Comma separated list of coverage/bin files that will be used for testing each computed alpha matrix. These files should be independent from the training set (--inputFilesTrain) for measuring underfitting/overfitting. (They can be other chromosomes)')
    parser.add_argument('--outputDir', type=str, default = "tune_alpha",
                        help='Output directory')
    parser.add_argument('--otherParamsText', type=str, default="",
                        help='Text file with one line that contains the optional parameters other than input file, model type and output dir')
    parser.add_argument('--numberOfStartPoints', type=int, default=10,
                        help='EGO algorithm needs a set of start points to compute an initial estimation of the objective function. These points will be selected randomly in the feasibility space of the alpha matrix. (Default=10)')
    parser.add_argument('--lowerBound', type=float, default=0.0,
                        help='Lower bound on each entry of alpha matrix (Default=0.0)')
    parser.add_argument('--upperBound', type=float, default=0.8,
                        help='Upper bound on each entry of alpha matrix (Default=0.8)')
    parser.add_argument('--iterations', type=int, default=50,
                        help='Number of iterations for EGO algorithm (After computing start points) (Default=50)')
    parser.add_argument('--modelType', type=str, default='gaussian',
                        help='modelType can be either gaussian or trunc_exp_gaussian')
    parser.add_argument('--annotationLabel', type=str, default="whole_genome",
                        help='Annotation label whose score will be used for optimizing alpha matrix (Default=whole_genome)')
    parser.add_argument('--sizeLabel', type=str, default="ALL_SIZES",
                        help='Size label whose score will be used for optimizing alpha matrix (Default=ALL_SIZES)')
    parser.add_argument('--maxJobs', type=int, default=8,
                        help='Maximum number of hmm_flagger jobs that can be run at the same time. It is useful to adjust this number based on total number of cores available and the number of threads that each hmm_flagger job will take.')
    parser.add_argument('--egoCriterion', type=str, default="EI",
                        help='Criterion for EGO. It can be one of EI, SBO or LCB. https://smt.readthedocs.io/en/latest/_src_docs/applications/ego.html#implementation-notes (Default = EI)')

    args = parser.parse_args()
    lowerBound = args.lowerBound
    upperBound = args.upperBound
    numberOfStartPoints = args.numberOfStartPoints
    numberOfIterations = args.iterations
    criterion = args.egoCriterion

    functionToMinimizeInternal.pointIndex = 0
    functionToMinimizeInternal.allScoresTrain = []
    functionToMinimizeInternal.allScoresTest = []
    functionToMinimizeInternal.numberOfStartPoints = numberOfStartPoints
    functionToMinimizeInternal.numberOfIterations = numberOfIterations
    functionToMinimizeInternal.paramsDict = {}
    functionToMinimizeInternal.paramsDict["inputPathTrainList"] = args.inputFilesTrain.strip().split(',')
    functionToMinimizeInternal.paramsDict["inputPathTestList"] = args.inputFilesTest.strip().split(',')
    trainSize = len(functionToMinimizeInternal.paramsDict["inputPathTrainList"])
    testSize = len(functionToMinimizeInternal.paramsDict["inputPathTestList"])
    functionToMinimizeInternal.paramsDict["outputDir"] = args.outputDir
    functionToMinimizeInternal.paramsDict["modelType"] = args.modelType
    functionToMinimizeInternal.paramsDict["maxJobs"] = args.maxJobs
    functionToMinimizeInternal.paramsDict["annotationLabel"] = args.annotationLabel
    functionToMinimizeInternal.paramsDict["sizeLabel"] = args.sizeLabel

    if args.otherParamsText.strip() != "":
        with open(args.otherParamsText.strip(), "r") as f:
            lines = f.readlines()
            if 0 < len(lines):
                functionToMinimizeInternal.paramsDict["otherParamsString"] = lines[0]


    startX = getStartPoints(lowerBound, upperBound, numberOfStartPoints)

    xLimits = getXLimits(lowerBound, upperBound)

    randomState = 42  # for reproducibility
    designSpace = DesignSpace(xLimits, random_state=randomState)
    
    criterion = "EI"  #'EI' or 'SBO' or 'LCB'

    ego = EGO(
        n_iter=numberOfIterations,
        criterion=criterion,
        xdoe=startX,
        surrogate=KRG(design_space=designSpace, print_global=False),
        random_state=randomState,
        n_parallel = 1
    )

    x_opt, y_opt, _, x_data, y_data = ego.optimize(fun=functionToMinimize)
    y_opt = -1 * y_opt

    optimumAlphaMatrix = convertXToAlphaMatrix(x_opt) 
    print(f"[{datetime.datetime.now()}] Optimum in Alpha =",file=sys.stderr)
    np.savetxt(sys.stderr, optimumAlphaMatrix, delimiter='\t', fmt="%.3f")
    print(f"with f(x)={y_opt}", file=sys.stderr)

    os.makedirs(args.outputDir, exist_ok=True)
    alphaTsvPath = os.path.join(args.outputDir, "alpha_optimum.tsv")

    print(f"[{datetime.datetime.now()}] Saving optimum alpha matrix in {alphaTsvPath}", file=sys.stderr)
    saveAlphaMatrixInTSV(optimumAlphaMatrix, alphaTsvPath)

    optimumIndices = []
    scoresTrainTsvPath = os.path.join(args.outputDir, "scores_train.tsv")
    print(f"[{datetime.datetime.now()}] Saving scores table for train files in {scoresTrainTsvPath}", file=sys.stderr)
    with open(scoresTrainTsvPath, "w") as f:
        # write header line
        f.write("Point_Index\tPoint_Type\tFinal_Score\tIs_Optimum_By_Train")
        for inputIndex in range(1, trainSize + 1):
            f.write("\tOverlap_Based_{inputIndex}\tBase_Level_{inputIndex}\tAuN_Based_{inputIndex}\tCombined_{inputIndex}")
        f.write("\n")

        for i, scores in enumerate(functionToMinimizeInternal.allScoresTrain):
            div = abs(y_opt[0] - scores[1]) / (y_opt[0])
            if div < 1e-4:
                isOptimum = "YES"
                optimumIndices.append(i)
            else:
                isOptimum = "NO"

            f.write(f"{i + 1}\t{scores[0]}\t{scores[1]:.3f}\t{isOptimum}")
            for overlapBasedScore, baseLevelScore, auNScore, combinedScore in zip(scores[2], scores[3], scores[4], scores[5]):
                f.write(f"\t{overlapBasedScore:.3f}\t{baseLevelScore:.3f}\t{auNScore:.3f}\t{combinedScore:.3f}")
            f.write("\n")


    if 0 < testSize:
        scoresTestTsvPath = os.path.join(args.outputDir, "scores_test.tsv")
        print(f"[{datetime.datetime.now()}] Saving scores table for test files in {scoresTestTsvPath}", file=sys.stderr)
        with open(scoresTestTsvPath, "w") as f:
            # write header line
            f.write("Point_Index\tPoint_Type\tFinal_Score\tIs_Optimum_By_Train")
            for inputIndex in range(1, testSize + 1):
                f.write("\tOverlap_Based_{inputIndex}\tBase_Level_{inputIndex}\tAuN_Based_{inputIndex}\tCombined_{inputIndex}")
            f.write("\n")

            for i, scores in enumerate(functionToMinimizeInternal.allScoresTest):
                isOptimum = "YES" if i in optimumIndices else "NO"

                f.write(f"{i + 1}\t{scores[0]}\t{scores[1]:.3f}\t{isOptimum}")
                for overlapBasedScore, baseLevelScore, auNScore, combinedScore in zip(scores[2], scores[3], scores[4], scores[5]):
                    f.write(f"\t{overlapBasedScore:.3f}\t{baseLevelScore:.3f}\t{auNScore:.3f}\t{combinedScore:.3f}")
                f.write("\n")

    print(f"[{datetime.datetime.now()}] Done!", file=sys.stderr)

main()
