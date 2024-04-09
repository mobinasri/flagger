import numpy as np
import pandas as pd
from numpy.random import choice
import argparse
from scipy.stats import truncexpon

def parseTransitionMatrixPerRegion(pathToTsv):
    transitionTable = pd.read_csv(pathToTsv, sep="\t")
    #print(transitionTable)
    numberOfRegions = len(set(transitionTable["#Region"]))

    transitionMatrixPerRegion = []
    for region in range(numberOfRegions):
        offset = region * 5
        transitionMatrix = np.array(transitionTable.iloc[offset:offset + 4, 2:6])
        row_sums = np.abs(transitionMatrix).sum(axis=1)
        transitionMatrix /= row_sums[:, np.newaxis] # normalize by row sum
        #print(transitionMatrix)
        transitionMatrixPerRegion.append(transitionMatrix)

    return transitionMatrixPerRegion


def parseEmissionParametersPerRegion(pathToTsv):
    emissionTable = pd.read_csv(pathToTsv, sep="\t")
    indices = [emissionTable["#State"][i] + "_" + emissionTable["Parameter"][i] for i in range(emissionTable.shape[0])]
    emissionTable.index = indices
    numberOfRegions = emissionTable.shape[1] - 4

    emissionParametersPerRegion = []
    for region in range(numberOfRegions):
        emissionParametersPerState = []
        for state in ["Err", "Dup", "Hap", "Col"]:
            parameters = {"Dist": None, "Mean":[], "Var":[], "Weight":[], "Comp":[]}
            parameters["Dist"] = emissionTable[f"Distribution"][f"{state}_Mean"]
            parameters["Comp"] = int(emissionTable[f"Components"][f"{state}_Mean"])
            parameters["Mean"] = np.array(emissionTable[f"Values_Region_{region}"][f"{state}_Mean"].split(","),dtype=np.float64)
            if parameters["Dist"] != "Truncated Exponential":
                parameters["Var"] = np.array(emissionTable[f"Values_Region_{region}"][f"{state}_Var"].split(","),dtype=np.float64)
                parameters["Weight"] = np.array(emissionTable[f"Values_Region_{region}"][f"{state}_Weight"].split(","),dtype=np.float64)
            else:
                parameters["Trunc_Point"] = np.array(emissionTable[f"Values_Region_{region}"][f"{state}_Trunc_Point"].split(","),dtype=np.float64)
            emissionParametersPerState.append(parameters)
        emissionParametersPerRegion.append(emissionParametersPerState)
    return emissionParametersPerRegion


def getThetaAndR(mean, var):
  theta = mean / var
  r = mean ** 2 / (var - mean)
  return theta, r

def generateNegativeBinomialObservation(means , variances, weights):
    comp = choice(np.arange(len(weights)), 1, p=weights)[0]
    theta, r = getThetaAndR(means[comp], variances[comp])
    obs = np.random.negative_binomial(r, theta, 1)[0]
    return obs

def generateGaussianObservation(means , variances, weights):
    comp = choice(np.arange(len(weights)), 1, p=weights)[0]
    mean = means[comp]
    sigma = np.sqrt(variances[comp])
    obs = np.round(np.random.normal(mean, sigma, 1)[0])
    return 0 if obs <= 0 else obs

def generateTruncatedExponentialObservation(mean, truncPoint):
    obs = np.round(truncexpon.rvs(scale=mean, loc=0, b=truncPoint, size=1)[0])
    return 0 if obs <= 0 else obs


def generateOneObservation(emissionParameters):
    distName = emissionParameters["Dist"]
    if distName == "Negative Binomial":
        return generateNegativeBinomialObservation(emissionParameters["Mean"], emissionParameters["Var"], emissionParameters["Weight"])
    elif distName == "Gaussian":
        return generateGaussianObservation(emissionParameters["Mean"], emissionParameters["Var"], emissionParameters["Weight"])
    elif distName == "Truncated Exponential":
        return generateTruncatedExponentialObservation(emissionParameters["Mean"][0], emissionParameters["Trunc_Point"][0])
    else:
        print("[Error] Distribution column should be from these three options; 'Negative Binomial', 'Gaussian', 'Truncated Exponential'")
        exit()

def getRegionAndStateIndex(transitionMatrixPerRegion, preRegion, preState, regionChangeRate):
    numberOfRegions = len(transitionMatrixPerRegion)
    numberOfStates = transitionMatrixPerRegion[0].shape[0]
    if preRegion == None:
        # get one region and state randomly
        region = choice(np.arange(numberOfRegions), 1)[0]
    else:
        if 1 < numberOfRegions: 
            regionWeights = np.ones(numberOfRegions) * 1 / (numberOfRegions - 1) * regionChangeRate
            regionWeights[preRegion] = 1.0 - regionChangeRate
        else:
            regionWeights = np.ones(1)
        # get the new region
        region = choice(np.arange(numberOfRegions), 1, p=regionWeights)[0]
    
    # if the region has not changed
    if preRegion == region:
        state = choice(np.arange(numberOfStates), 1, p=transitionMatrixPerRegion[region][preState])[0]
        return region, state
    else:
        state = np.random.randint(numberOfStates, size=1)[0] # uniformly random

    return region, state

def generateObservations(transitionMatrixPerRegion, emissionParametersPerRegion, numberOfObservations, regionChangeRate):
    region = None
    state = None
    observations = []
    regions = []
    states = []
    for i in range(numberOfObservations):
        region, state = getRegionAndStateIndex(transitionMatrixPerRegion,
                                               preRegion=region, 
                                               preState=state,
                                               regionChangeRate=regionChangeRate)
        obs = generateOneObservation(emissionParametersPerRegion[region][state])
        observations.append(obs)
        regions.append(region)
        states.append(state)
    return regions, states, observations

def writeObservations(regions, states, observations, pathToWrite):
    with open(pathToWrite,"w") as f:
        f.write(f"#Total\tHigh_Mapq\tHighly_Clipped\tRegion\tState\n")
        for pos in range(len(regions)):
            obs = observations[pos]
            state = states[pos]
            region = regions[pos]
            f.write(f"{obs}\t0.0\t0.0\t{region}\t{state}\n")

def main():
    parser = argparse.ArgumentParser(description='Simulate coverage data for running hmm_test.')
    parser.add_argument('--pathToEmission', type=str,
                    help='Path to the tsv file that contains emission parameters for different states and regions.')
    parser.add_argument('--pathToTransition', type=str,
                    help='Path to the tsv file that contains transition matrices for different regions.')
    parser.add_argument('--pathOutput', type=str, default="observations.tsv",
                    help='Path for the output tsv file that contains the simulated observations, states and regions. (Default= "observations.tsv")')
    parser.add_argument('--numberOfObservations', type=int, default=10000,
                    help='Total number of observations for simulation.(Default = 10000)')
    parser.add_argument('--regionChangeRate', type=float, default=0.001,
                    help='Rate of changing regions (will be ignored if there is only one region). (Default= 0.001)')

    # Fetch the arguments
    args = parser.parse_args()
    pathOutput = args.pathOutput
    pathToEmission = args.pathToEmission
    pathToTransition = args.pathToTransition
    numberOfObservations = args.numberOfObservations
    regionChangeRate = args.regionChangeRate

    # parse emission parameters and transition matrices
    emissionParametersPerRegion = parseEmissionParametersPerRegion(pathToEmission)
    transitionMatrixPerRegion = parseTransitionMatrixPerRegion(pathToTransition)

    # simulate observations
    regions, states, observations = generateObservations(transitionMatrixPerRegion, 
                                                         emissionParametersPerRegion, 
                                                         numberOfObservations=numberOfObservations, 
                                                         regionChangeRate=regionChangeRate)
    # write observations, states and regions in a tsv file
    writeObservations(regions, states, observations, pathToWrite=pathOutput)

if __name__ == "__main__": main()
