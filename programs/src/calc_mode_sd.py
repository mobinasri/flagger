import sys
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description='Get the coverage file and calculate mode and sd')
    parser.add_argument('--countsInput', type=str,
                    help='Input counts file')
    parser.add_argument('--minCoverage', type=int, default=5,
                    help='min coverage value for calculating mode coverage')
    parser.add_argument('--modeOutput', type=str,
                    help='File path to save coverage mode')
    parser.add_argument('--sdOutput', type=str,
                    help='File path to save coverage standard deviation')
    args = parser.parse_args()
    countsPath = args.countsInput
    minCoverage = args.minCoverage
    modeOutputPath = args.modeOutput
    sdOutputPath = args.sdOutput

    # Calculate the mod
    totalSize = 0
    totalCoverage = 0
    coverageMode = 0
    maxCount  = 0
    with open(countsPath,'r') as f:
        for line in f:
            attrbs = line.strip().split()
            coverage = int(attrbs[0])
            count = int(attrbs[1])
            totalSize += count
            totalCoverage += count * coverage 
            if coverage < minCoverage:
                continue
            if maxCount < count:
                coverageMode = coverage 
                maxCount = count

    coverageMean = totalCoverage / totalSize
    # Calculate the standard deviation
    totalSquaredCoverage = 0
    totalSize = 0
    with open(countsPath,'r') as f:
        for line in f:
            attrbs = line.strip().split()
            coverage = int(attrbs[0])
            count = int(attrbs[1])
            # Exclude very large coverages in calculating standard deviation
            if (coverage < (2.5 * coverageMode)):
                totalSquaredCoverage += ((coverage - coverageMean) ** 2) * count
                totalSize += count
    coverageSD = math.sqrt(totalSquaredCoverage/totalSize)

    with open(modeOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageMode))
    with open(sdOutputPath,'w') as f:
        f.write("{:.2f}\n".format(coverageSD))




main()
