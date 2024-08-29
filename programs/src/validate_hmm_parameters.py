import numpy as np
import pandas as pd
import argparse
from simulate_coverage_data import parseTransitionMatrixPerRegion, parseEmissionParametersPerRegion  

def yieldRelativeDiffForOneDist(emissionParametersTruth, emissionParametersQuery):
    distName = emissionParametersTruth["Dist"]
    if distName == "Negative Binomial" or distName == "Gaussian":
        params = ["Mean", "Var", "Weight"]
    elif distName == "Truncated Exponential":
        params = ["Mean", "Trunc_Point"]
    else:
        print("[Error] Distribution column should be from these three options; 'Negative Binomial', 'Gaussian', 'Truncated Exponential'")
        exit()

    for param in params:
        for queryValue, truthValue in zip(emissionParametersQuery[param], emissionParametersTruth[param]):
            if truthValue == 0.0:
                relDiff = 0.0
            else:
                relDiff = abs(queryValue - truthValue) / abs(truthValue)
            yield distName, param, relDiff




def printRelativeDiffForEmission(emissionParametersPerRegionTruth, emissionParametersPerRegionQuery, numberOfRegions):
    print(f"\n#Emission Parameters:\n")
    print(f"#Region\tState\tDistribution\tParameter\tRel_Diff\tLower_0.1")
    for region in range(numberOfRegions):
        for stateIndex, state in enumerate(["Err", "Dup", "Hap", "Col"]):
            for distName, param, relDiff in yieldRelativeDiffForOneDist(emissionParametersPerRegionTruth[region][stateIndex], 
                                                                        emissionParametersPerRegionQuery[region][stateIndex]):
                if relDiff < 0.1:
                    converged = "YES"
                else:
                    converged = "NO"
                print(f"{region}\t{state}\t{distName}\t{param}\t{relDiff:.3f}\t{converged}")

def printRelativeDiffForTransition(transitionMatrixPerRegionTruth, transitionMatrixPerRegionQuery, numberOfRegions):
    print(f"\n#Transition Probabilites:\n")
    print(f"#Region\tStateFrom\tStateTo\tRel_Diff\tLower_0.1")
    for region in range(numberOfRegions):
        diffMatrix = np.abs(transitionMatrixPerRegionTruth[region] - transitionMatrixPerRegionQuery[region])
        truthMatrixNonZero = np.abs(transitionMatrixPerRegionTruth[region])
        diffMatrix[truthMatrixNonZero == 0.0] = 0.0
        truthMatrixNonZero[truthMatrixNonZero == 0.0] = 1.0
        relDiffMatrix = diffMatrix / truthMatrixNonZero
        for stateIndex1, state1 in enumerate(["Err", "Dup", "Hap", "Col"]):
            for stateIndex2, state2 in enumerate(["Err", "Dup", "Hap", "Col"]):
                relDiff = relDiffMatrix[stateIndex1][stateIndex2]
                if relDiff < 0.1:
                    converged = "YES"
                else:
                    converged = "NO"
                print(f"{region}\t{state1}\t{state2}\t{relDiff:.3f}\t{converged}")


def main():
    parser = argparse.ArgumentParser(description='Calculate relative differences between truth HMM parameter values used for generating simulated coverage data and the (query) parameter values fit by EM algorithm in hmm_test program.')
    parser.add_argument('--pathToEmissionTruth', type=str,
                    help='Path to the tsv file that contains TRUTH emission parameters for different states and regions.')
    parser.add_argument('--pathToTransitionTruth', type=str,
                    help='Path to the tsv file that contains TRUTH transition matrices for different regions.')
    parser.add_argument('--pathToEmissionQuery', type=str,
                    help='Path to the tsv file that contains QUERY emission parameters for different states and regions.')
    parser.add_argument('--pathToTransitionQuery', type=str,
                    help='Path to the tsv file that contains QUERY transition matrices for different regions.')

    # Fetch the arguments
    args = parser.parse_args()
    pathToEmissionTruth = args.pathToEmissionTruth
    pathToTransitionTruth = args.pathToTransitionTruth
    pathToEmissionQuery = args.pathToEmissionQuery
    pathToTransitionQuery = args.pathToTransitionQuery

    # parse emission parameters and transition matrices
    emissionParametersPerRegionTruth = parseEmissionParametersPerRegion(pathToEmissionTruth)
    transitionMatrixPerRegionTruth = parseTransitionMatrixPerRegion(pathToTransitionTruth)

    emissionParametersPerRegionQuery = parseEmissionParametersPerRegion(pathToEmissionQuery)
    transitionMatrixPerRegionQuery = parseTransitionMatrixPerRegion(pathToTransitionQuery)

    numberOfRegions = len(emissionParametersPerRegionTruth)

    # print relative change
    printRelativeDiffForEmission(emissionParametersPerRegionTruth, emissionParametersPerRegionQuery, numberOfRegions)
    printRelativeDiffForTransition(transitionMatrixPerRegionTruth, transitionMatrixPerRegionQuery, numberOfRegions)

if __name__ == "__main__": main()
