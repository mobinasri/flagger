import sys
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description='Get the coverage file and calculate mean and sd')
    parser.add_argument('--countsInput', type=str,
                    help='Input counts file')
    parser.add_argument('--meanOutput', type=str,
                    help='File path to save coverage mean')
    parser.add_argument('--sdOutput', type=str,
                    help='File path to save coverage standard deviation')
    args = parser.parse_args()
    countsPath = args.countsInput
    meanOutputPath = args.meanOutput
    sdOutputPath = args.sdOutput

    # Calculate the mean
    totalSize = 0
    totalCoverage = 0
    with open(countsPath,'r') as f:
        for line in f:
            attrbs = line.strip().split()
            coverage = int(attrbs[0])
            count = int(attrbs[1])
            totalCoverage += coverage * count
            totalSize += count
    coverageMean = totalCoverage / totalSize

    # Calculate the standard deviation
    totalSquaredCoverage = 0
    with open(countsPath,'r') as f:
        for line in f:
            attrbs = line.strip().split()
            coverage = int(attrbs[0])
            count = int(attrbs[1])
            # Exclude very large coverages in calculating standard deviation
            if (coverage < (2.5 * coverageMean)):
                totalSquaredCoverage += ((coverage - coverageMean) ** 2) * count
    coverageSD = math.sqrt(totalSquaredCoverage/totalSize)

    with open(meanOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageMean))
    with open(sdOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageSD))




main()
