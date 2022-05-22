from coverage_frequency_distribution import CoverageDistribution
from collections import defaultdict
import sys
import argparse
import operator

def main():
    parser = argparse.ArgumentParser(description='Fit coverage model')
    parser.add_argument('--counts', type=str, help='coverage counts file')
    parser.add_argument('--output', type=str, help='probability table output')
    args = parser.parse_args()
    countsPath = args.counts
    outputPath = args.output

    counts = defaultdict(int)
    with open(countsPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            cov = int(attrbs[0])
            freq = int(attrbs[1])
            counts[cov] = freq
    
    # fit the model
    cov_model = CoverageDistribution.fit(counts, tol = 1e-4, max_iters=1000)
    # find the fitted distribution
    X = list(range(max(counts)+1))
    Y = [cov_model.pdf(x) for x in X]
    hap_cov = max(counts.items(), key=operator.itemgetter(1))[0]
    scale_factor = counts[hap_cov] / Y[hap_cov]
    Y = [scale_factor * y for y in Y]
    # print the probability table
    with open(outputPath,"w+") as f:
        f.write("#coverage\tfreq\tfitted\terror\thaploid\tcollapsed\n")
        for cov in range(len(counts)):
            prob_err = cov_model.probability_erroneous(cov)
            prob_correct = cov_model.probability_haploid(cov)
            prob_collapsed = cov_model.probability_collapsed(cov)
            f.write("{:d}\t{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(cov, counts[cov], Y[cov], prob_err, prob_correct, prob_collapsed))

if __name__ == "__main__":
    main()

