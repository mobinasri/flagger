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
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN"]
    return float(table["F1-Score"].item())

def getBaseLevelScore(outputDir, annotationLabel, sizeLabel):
    tsvPath = [os.path.join(outputDir, a) for a in os.listdir(outputDir) if "final.benchmarking.tsv" in a][0]
    table = pd.read_csv(tsvPath,sep="\t")
    table = table.rename(columns={"#Metric_Type": "Metric_Type"})
    table = table.loc[lambda table: table.Metric_Type == "base_level"]
    table = table.loc[lambda table: table.Category_Name == annotationLabel]
    table = table.loc[lambda table: table.Size_Bin_Name == sizeLabel]
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN"]
    return float(table["F1-Score"].item())

def getContiguityScore(outputDir, annotationLabel, sizeLabel):
    tsvPath = [os.path.join(outputDir, a) for a in os.listdir(outputDir) if "final.benchmarking.auN_ratio.tsv" in a][0]
    table = pd.read_csv(tsvPath,sep="\t")
    table = table.loc[lambda table: table.Category_Name == annotationLabel]
    table = table.loc[lambda table: table.Size_Bin_Name == sizeLabel]
    table = table.loc[lambda table: table.Label == "HARMONIC_MEAN"]
    return float(table["auN_Ratio"].item())

def getCombinedScore(outputDir, annotationLabel, sizeLabel):
    score1 = getOverlapBasedScore(outputDir, annotationLabel, sizeLabel)
    score2 = getBaseLevelScore(outputDir, annotationLabel, sizeLabel)
    score3 = getContiguityScore(outputDir, annotationLabel, sizeLabel)
    return (score1 + score2 + score3) / 3


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
    inputPathList = functionToMinimizeInternal.paramsDict['inputPathList']
    modelType = functionToMinimizeInternal.paramsDict['modelType']
    otherParamsString = functionToMinimizeInternal.paramsDict['otherParamsString']
    maxJobs = functionToMinimizeInternal.paramsDict['maxJobs']
    annotationLabel = functionToMinimizeInternal.paramsDict['annotationLabel']
    sizeLabel = functionToMinimizeInternal.paramsDict['sizeLabel']
    pointIndex = functionToMinimizeInternal.pointIndex

    alphaMatrix = convertXToAlphaMatrix(x)
    os.makedirs(f"{outputDir}/optimization_point_{pointIndex}", exist_ok=True)
    alphaTsvPath = f"{outputDir}/optimization_point_{pointIndex}/alpha_mat.tsv"
    saveAlphaMatrixInTSV(alphaMatrix, alphaTsvPath)

    internalOutputDirList = []
    paramsStringList = []
    for inputPath in inputPathList:
        prefix = getInputPrefix(inputPath)
        internalOutputDir = f"{outputDir}/optimization_point_{pointIndex}/{prefix}"
        os.makedirs(internalOutputDir, exist_ok=True)
        internalOutputDirList.append(internalOutputDir)
        paramsString = getParametersString(alphaTsvPath, inputPath, internalOutputDir, modelType, otherParamsString)
        paramsStringList.append(paramsString)

    runHMMFlaggerForList(paramsStringList, maxJobs)

    combinedScoreList = []
    for internalOutputDir in internalOutputDirList:
        combinedScoreList.append(getCombinedScore(internalOutputDir, annotationLabel, sizeLabel))

    # times -1 since EGO algorithm minimizes the objective function
    finalScore = -1 * sum(combinedScoreList) / len(combinedScoreList)

    print(f"[{datetime.datetime.now()}] Alpha = {alphaMatrix} , Final Score = {finalScore}", file=sys.stderr)
    sys.stderr.flush()

    return finalScore


def runHMMFlaggerForList(paramsStringList, maxJobs):
    success = 0
    with ThreadPoolExecutor(max_workers=maxJobs) as executor:
        futures = [executor.submit(runHMMFlagger, paramsString) for paramsString in paramsStringList]
        for future in as_completed(futures):
            # retrieve the result
            result = future.result()
            print(result[0])
            if result[1] == 0:
                success += 1
    print(f"[{datetime.datetime.now()}] Successfully ran {len(paramsStringList)} jobs for hmm_flagger", file=sys.stderr)
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
    parser.add_argument('--inputFiles', type=str,
                        help='Comma separated list of coverage/bin files that will be used for tuning alpha matrix')
    parser.add_argument('--outputDir', type=str, default = "tune_alpha",
                        help='Output directory')
    parser.add_argument('--otherParamsText', type=str, default="",
                        help='Text file with one line that contains the optional parameters other than input file, model type and output dir')
    parser.add_argument('--numberOfStartPoints', type=int, default=10,
                        help='EGO algorithm needs a set of start points to compute an initial estimation of the objective function. These points will be selected randomly in the feasibility space of the alpha matrix.')
    parser.add_argument('--lowerBound', type=float, default=0.0,
                        help='Lower bound on each entry of alpha matrix')
    parser.add_argument('--upperBound', type=float, default=0.8,
                        help='Upper bound on each entry of alpha matrix')
    parser.add_argument('--iterations', type=int, default=50,
                        help='Number of iterations for EGO algorithm (After computing start points)')
    parser.add_argument('--modelType', type=str, default='gaussian',
                        help='modelType can be either gaussian or trunc_exp_gaussian')
    parser.add_argument('--annotationLabel', type=str, default="whole_genome",
                        help='Annotation label whose score will be used for optimizing alpha matrix')
    parser.add_argument('--sizeLabel', type=str, default="ALL_SIZES",
                        help='Size label whose score will be used for optimizing alpha matrix')
    parser.add_argument('--maxJobs', type=str, default=8,
                        help='Maximum number of hmm_flagger jobs that can be run at the same time. It is useful to adjust this number based on total number of cores available and the number of threads that each hmm_flagger job will take.')

    args = parser.parse_args()
    lowerBound = args.lowerBound
    upperBound = args.upperBound
    numberOfStartPoints = args.numberOfStartPoints
    numberOfIterations = args.iterations

    functionToMinimizeInternal.pointIndex = 0
    functionToMinimizeInternal.paramsDict = {}
    functionToMinimizeInternal.paramsDict["inputPathList"] = args.inputFiles.strip().split(',')
    functionToMinimizeInternal.paramsDict["outputDir"] = args.outputDir
    functionToMinimizeInternal.paramsDict["modelType"] = args.modelType
    functionToMinimizeInternal.paramsDict["maxJobs"] = args.maxJobs
    functionToMinimizeInternal.paramsDict["annotationLabel"] = args.annotationLabel
    functionToMinimizeInternal.paramsDict["sizeLabel"] = args.sizeLabel

    if args.otherParamsText.strip() is not "":
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
    print("Minimum in x={} with f(x)={}".format(x_opt, y_opt))