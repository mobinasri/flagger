from coverage_frequency_distribution_extra_constrained_poisson import CoverageDistribution
from collections import defaultdict
import sys
import argparse
import operator

def main():
    parser = argparse.ArgumentParser(description='Fit coverage model')
    parser.add_argument('--counts', type=str, help='coverage counts file')
    parser.add_argument('--cov', type=float, help='expected haploid coverage (Optional, will be inferred from the distribution if not used)', default=None)
    parser.add_argument('--output', type=str, help='probability table output')
    args = parser.parse_args()
    countsPath = args.counts
    haploidCov = args.cov
    outputPath = args.output

    counts = defaultdict(int)
    with open(countsPath,"r") as f:
        for line in f:
            attrbs = line.strip().split()
            cov = int(attrbs[0])
            freq = int(attrbs[1])
            counts[cov] = freq
   
    if sum(counts.values()) == 0:
        with open(outputPath,"w+") as f:
            f.write("#coverage\tfreq\tfit\terror\textra\thaploid\tcollapsed\n")
            for cov in range(len(counts)):
                f.write("{:d}".format(cov) + "\t0" * 6 + "\n")
        return
    # fit the model
    cov_model = CoverageDistribution.fit(counts, tol = 1e-4, max_iters=500, start_point_cov_per_ploidy= haploidCov)
    # find the fitted distribution
    X = list(range(max(counts)+1))
    Y = [cov_model.pdf(x) for x in X]
    hap_cov = max(counts.items(), key=operator.itemgetter(1))[0]
    scale_factor = counts[hap_cov] / Y[hap_cov]
    Y = [scale_factor * y for y in Y]
    # print the probability table
    with open(outputPath,"w+") as f:
        f.write("#coverage\tfreq\tfit\terror\tduplicated\thaploid\tcollapsed\n")
        for cov in range(len(counts)):
            prob_err = cov_model.probability_erroneous(cov)
            prob_extra = cov_model.probability_extra(cov)
            prob_correct = cov_model.probability_haploid(cov)
            prob_collapsed = cov_model.probability_collapsed(cov)
            f.write("{:d}\t{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(cov, counts[cov], Y[cov], prob_err, prob_extra, prob_correct, prob_collapsed))

if __name__ == "__main__":
    main()

