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


def generateOneObservation(emissionParameters, alpha, preObservation):
    distName = emissionParameters["Dist"]
    if distName == "Negative Binomial":
        return generateNegativeBinomialObservation(emissionParameters["Mean"], emissionParameters["Var"], emissionParameters["Weight"])
    elif distName == "Gaussian":
        if preObservation is not None and alpha is not None:
            mean = emissionParameters["Mean"] * (1 - alpha) + preObservation * alpha
        else:
            mean = emissionParameters["Mean"]
        return generateGaussianObservation(mean, emissionParameters["Var"], emissionParameters["Weight"])
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

def generateObservations(transitionMatrixPerRegion, emissionParametersPerRegion, numberOfObservations, regionChangeRate, alphaMatrix):
    region = None
    state = None
    observations = []
    regions = []
    states = []
    preObs = None
    preState = None
    for i in range(numberOfObservations):
        region, state = getRegionAndStateIndex(transitionMatrixPerRegion,
                                               preRegion=region, 
                                               preState=state,
                                               regionChangeRate=regionChangeRate)
        if i == 0 or preState == None or preObs == None or alphaMatrix == None:
            obs = generateOneObservation(emissionParametersPerRegion[region][state], 0, None)
        else:
            obs = generateOneObservation(emissionParametersPerRegion[region][state], alphaMatrix[preState][state], preObs)
        observations.append(obs)
        regions.append(region)
        states.append(state)
        preObs = obs
        preState = state
    return regions, states, observations


def createAlphaMatrix(alphaString):
    if alphaString == None or len(alphaString) == 0:
        return None
    alphaValues = [float(num) for num in alphaString.strip().split(',')]
    numberOfStates = len(alphaValues) - 1
    alphaMatrix = np.zeros((numberOfStates, numberOfStates))

    alphaTrans = alphaValues[-1]
    alphaMatrix.fill(alphaTrans)

    for state in range(numberOfStates):
        alphaMatrix[state][state] = alphaValues[state]
    return alphaMatrix

def writeObservationsIntoCov(regions, states, observations, regionCoverages, contigLengths, pathToWrite):
    numberOfRegions = len(regionCoverages)
    numberOfObservations  = len(regions)
    with open(pathToWrite,"w") as f:
        f.write(f"#annotation:len:{2+len(contigLengths)}\n")
        f.write(f"#annotation:name:0:no_annotation\n")
        f.write(f"#annotation:name:1:whole_genome\n")
        for i in range(len(contigLengths)):
            f.write(f"#annotation:name:{2 + i}:TEST_CONTIG_{i}\n")
        f.write(f"#region:len:{numberOfRegions}\n")
        for i in range(numberOfRegions):
            f.write(f"#region:coverage:{i}:{regionCoverages[i]}\n")
        f.write(f"#label:len:4\n")
        f.write(f"#truth:true\n")
        f.write(f"#prediction:false\n")

        start = 0
        for i, contigLength in enumerate(contigLengths):
            f.write(f">TEST_CONTIG_{i} {contigLength}\n")
            for pos in range(start, start + contigLength):
                obs = observations[pos]
                truthState = states[pos]
                region = regions[pos]
                startInContig = pos - start + 1
                endInContig = startInContig
                #coverage, high_mapq_coverage, high_clip_coverage, annotation indices, region index, truth
                f.write(f"{startInContig}\t{endInContig}\t{obs}\t0.0\t0.0\t1,{2+i}\t{region}\t{truthState}\n")
            start += contigLength

def main():
    parser = argparse.ArgumentParser(description='Simulate coverage data for running hmm_flagger.')
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
    parser.add_argument('--contigLengths', type=str, default="",
                        help='A comma separated list of numbers. The sum of number of should be equal to the number of observations. For example for --numberOfObservations 100 use can pass --contigLengths 30,40,30  (Default= one contig covering all observations)')
    parser.add_argument('--alpha', type=str, default="",
                        help='The dependency factors of the current emission density to the previous emission (Only works for Gaussian). It should be a comma-separated string of 5 numbers for these states respectively err,dup,hap,col,trans. (trans is for transitioning from one state to a different one) [Default = all alpha factors set to 0]')


# Fetch the arguments
    args = parser.parse_args()
    pathOutput = args.pathOutput
    pathToEmission = args.pathToEmission
    pathToTransition = args.pathToTransition
    numberOfObservations = args.numberOfObservations
    regionChangeRate = args.regionChangeRate
    contigLengthsStr = args.contigLengths
    alphaString = args.alpha

    alphaMatrix = createAlphaMatrix(alphaString)

    contigLengths = [int(i) for i in contigLengthsStr.strip().split(',')]
    if len(contigLengths) == 0:
        contigLengths = [numberOfObservations]

    if sum(contigLengths) != numberOfObservations:
        print("Error: total length of contigs does not match the number of observations.")
        exit()



    # parse emission parameters and transition matrices
    emissionParametersPerRegion = parseEmissionParametersPerRegion(pathToEmission)
    transitionMatrixPerRegion = parseTransitionMatrixPerRegion(pathToTransition)

    # simulate observations
    regions, states, observations = generateObservations(transitionMatrixPerRegion, 
                                                         emissionParametersPerRegion, 
                                                         numberOfObservations=numberOfObservations, 
                                                         regionChangeRate=regionChangeRate,
                                                         alphaMatrix=alphaMatrix)
    # hap index is 2
    # we want Mean parameter
    # 0 means the first component (hap has just one component)
    regionCoverages = [params[2]["Mean"][0] for params in emissionParametersPerRegion]

    # write observations, states and regions in a tsv file
    writeObservationsIntoCov(regions, states, observations, regionCoverages, contigLengths, pathToWrite=pathOutput)

if __name__ == "__main__": main()