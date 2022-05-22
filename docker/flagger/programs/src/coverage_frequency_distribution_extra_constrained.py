#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 14:29:21 2021

@author: Jordan
"""

# Some statistics helper functions
from copy import deepcopy
from math import log, exp, sqrt, pi, ceil
def exp_pdf(x, lam):
    if x < 0.0:
        return 0.0
    else:
        return lam * exp(-lam * x)   
def trunc_exp_pdf(x, lam, b):
    if x < 0.0 or x > b:
        return 0.0
    else:
        return exp_pdf(x, lam) / (1.0 - exp(-lam * b)) 
def norm_pdf(x, mu, sigma):
    z = (x - mu) / sigma
    return 1.0 / (sigma * sqrt(2.0 * pi)) * exp(-z * z / 2.0)

# adapted from https://en.wikipedia.org/wiki/Golden-section_search
def golden_section_search(f, a, b, tol = 1e-6):
    
    invphi = (sqrt(5.0) - 1.0) / 2.0  # 1 / phi
    invphi2 = (3.0 - sqrt(5.0)) / 2.0 # 1 / phi^2
    
    h = b - a
    if h <= tol:
        return (b + a) / 2.0
    
    n = int(ceil(log(tol / h) / log(invphi)))
    
    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)
    
    for k in range(n-1):
        if yc > yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)
    
    if yc > yd:
        return (a + d) / 2.0
    else:
        return (c + b) / 2.0

##
# A class that represents a model of coverage frequencies in a phased diploid
# assembly, which is fit using data from a coverage frequency histogram in the 
# CoverageDistribution.fit() method
##
class CoverageDistribution:
    
    ############################
    # The downstream interface #
    ############################
    
    # Return the probability that an element in the assembly with coverage c is
    # actually absent in the underlying genome
    def probability_erroneous(self, c):
        return self.component_probabilities(c)[1]
    
    # Return the probability that an element in the assembly with coverage c is
    # present in the underlying genome and haploid (indicating that is it
    # correctly assembled)
    def probability_haploid(self, c):
        comps = self.component_probabilities(c)
        # comps[2] is the probability of the main haploid component
        return comps[2]
    def probability_extra(self, c):
        comps = self.component_probabilities(c)
        # comps[0] is the probability of the extra component
        return comps[0]
    
    # Return the probability that an element in the assembly with coverage c is
    # actually consists of multiple elements in the underlying genome collapsed
    # (note that it may be that another copy of this element is where the 
    # collapsing actually ocurred)
    def probability_collapsed(self, c):
        comps = self.component_probabilities(c)
        return 1 - comps[0] - comps[1] - comps[2]
    
    # Fit and return a CoverageDistribution using a coverage frequency histogram,
    # which should be provided as either a collections.Counter or dict whose keys
    # are coverages and whose values are the number of bases that have that
    # coverage
    @staticmethod
    def fit(spectrum, tol = 1e-6, max_iters = 10000):
        
        # copy the spectrum so we can modify it without breaking expectations
        spectrum = deepcopy(spectrum)
        
        # get a coarse estimate that we can use as a starting point for EM
        cov_dist = CoverageDistribution.heuristic_fit(spectrum)
        
        for iteration in range(max_iters):
            
            print(f"iteration {iteration}")
            print("\terr scale {}".format(cov_dist.err_scale))
            print("\tcov per ploidy {}".format(cov_dist.cov_per_ploidy))
            print("\tvar per ploidy {}".format(cov_dist.var_per_ploidy))
            #print("\tcov extra {}".format(cov_dist.cov_extra))
            #print("\tvar extra {}".format(cov_dist.var_extra))
            print("\tweight extra {}".format(cov_dist.mixture_weight_extra))
            print("\tlow ploidy weights {}".format(", ".join("{:.3e}".format(v) for v in cov_dist.mixture_weights[:5])))
            print("\thigh ploidy weights {}".format(", ".join("{:.3e}".format(v) for v in cov_dist.mixture_weights[-8:])))
            #c = sorted(spectrum)[:200]
            #plt.plot(c, [cov_dist.pdf(x) for x in c])
            #plt.plot(c, [spectrum[x] * cov_dist.pdf(round(cov_dist.cov_per_ploidy)) / spectrum[round(cov_dist.cov_per_ploidy)]  for x in c])
            #plt.show()
            assignment = {}
            assignment_extra = {}
            for cov in spectrum: 
                probs = cov_dist.component_probabilities(cov_dist.cov_adj(cov))
                assignment[cov] = probs[1:]
                assignment_extra[cov] = probs[0]
            new_mixture_weights = []
            
            weight_denom = sum(spectrum.values())
                        
            # gather sufficient statistics for error scale
            eff_N = 0.0
            eff_sum = 0.0
            weight_numer = 0.0
            for cov in spectrum:
                freq = spectrum[cov]
                prob = assignment[cov][0]
                eff_sum += freq * prob * cov_dist.cov_adj(cov)
                eff_N += freq * prob
                weight_numer += prob * freq
                          
            # maximize likelihood for a truncated exponential
            b = cov_dist.err_trunc_point()
            llike = lambda lam: eff_N * log(lam) - eff_N * log(1.0 - exp(-lam * b)) - eff_sum * lam
            rate = golden_section_search(llike, 0, b, tol / 100.0)
            new_error_scale = 1.0 / rate

            # component weight is weighted proportion assigned to the error component
            new_mixture_weights.append(weight_numer / weight_denom)


            # the haploid coverage is the weighted mean
            numer = 0.0
            denom = 0.0
            weight_numer = 0.0
            for cov in spectrum:
                # for extra mode
                freq = spectrum[cov]
                prob = assignment_extra[cov]
                numer += freq * prob * cov_dist.cov_adj(cov)
                denom += freq * prob * 0.5
                weight_numer += prob * freq
            new_mixture_weight_extra = weight_numer / weight_denom

            for ploidy in range(1, cov_dist.max_ploidy() + 1):
                weight_numer = 0.0
                for cov in spectrum:
                    freq = spectrum[cov]
                    prob = assignment[cov][ploidy]
                    numer += freq * prob * cov_dist.cov_adj(cov)
                    denom += prob * freq * ploidy
                    weight_numer += prob * freq
                new_mixture_weights.append(weight_numer / weight_denom)
             
            new_cov_per_ploidy = numer / denom
            
            # the haploid std deviation is the sqrt of the weighted average
            # square deviation, adjusted for the fact that the variance of
            # ploidy P should be P times the variance of the haploid
            numer = 0.0
            denom = 0.0
            for cov in spectrum:
                freq = spectrum[cov]
                # for extra mode
                prob = assignment_extra[cov]
                dev = cov_dist.cov_adj(cov) - 0.5 * new_cov_per_ploidy
                numer += freq * prob * dev * dev / 0.5
                denom += freq * prob
                for ploidy in range(1, cov_dist.max_ploidy() + 1):
                    prob = assignment[cov][ploidy]
                    dev = cov_dist.cov_adj(cov) - ploidy * new_cov_per_ploidy
                    numer += freq * prob * dev * dev / ploidy
                    denom += freq * prob
                    
            new_var_per_ploidy = numer / denom

            # for extra mode
            new_cov_extra = new_cov_per_ploidy * 0.5
           
            numer = 0.0
            denom = 0.0
            for cov in spectrum:
                freq = spectrum[cov]
                prob = assignment_extra[cov]
                dev = cov_dist.cov_adj(cov) - new_cov_extra
                numer += freq * prob * dev * dev
                denom += freq * prob
            new_var_extra = numer /denom

            # the next iteration's distrubition
            new_cov_dist = CoverageDistribution(new_cov_per_ploidy, new_var_per_ploidy,
                                                new_error_scale, new_mixture_weights,
                                                cov_dist.err_trunc_point(), new_cov_extra, new_var_extra, new_mixture_weight_extra)
            
            converged = cov_dist.converged(new_cov_dist, tol)
            cov_dist = new_cov_dist
            if converged:
                break
            
        return cov_dist
    
    ##############################
    # Internal/debugging methods #
    ##############################
    
    # TODO: generalize this for truncated distributions at the max value?
    
    @staticmethod
    def heuristic_fit(spectrum):
        
        max_cov = max(spectrum)
        
        cov_at_max_frequency = None
        max_frequency = -1
        has_increased = False
        freq_l = list(spectrum.values())
        max_spectrum = max(freq_l[1:])
        for cov_value in sorted(spectrum)[1:]:
            if cov_value - 1 in spectrum:
                if spectrum[cov_value] > spectrum[cov_value - 1]:
                    has_increased = True
            if spectrum[cov_value] > (0.75 * max_spectrum) and has_increased:
                max_frequency = spectrum[cov_value]
                cov_at_max_frequency = cov_value
        print(max_spectrum, cov_at_max_frequency) 
        cov_per_ploidy = cov_at_max_frequency
        var_per_ploidy = cov_per_ploidy
        
        cov_extra = cov_at_max_frequency / 2
        var_extra = cov_extra

        max_ploidy = int(max_cov / cov_per_ploidy)
        while max_ploidy * cov_per_ploidy - 3.0 * sqrt(max_ploidy * cov_per_ploidy) < max_cov:
            max_ploidy += 1
        
        # start with an arbitrary, small scale
        err_scale = 0.01
        
        # the denominator along with pseudocounts
        total_count = sum(spectrum[cov] for cov in spectrum) + max_ploidy + 1
        
        # approximate the mixture weights simplistically by hard-assigning all data
        # points to their nearest component
        mixture_weights = []
        # weight of the exponential component
        count = 1
        for cov in range(0, int(cov_extra / 2)):
            if cov in spectrum:
                count += spectrum[cov]
        mixture_weights.append(count/ total_count)
        # weight of the extra component
        count = 1
        for cov in range(int(cov_extra / 2), int((cov_extra + cov_per_ploidy)/2)):
            if cov in spectrum:
                count += spectrum[cov]
        mixture_weight_extra = count / total_count
        # weight of the main haploid component
        count = 1
        for cov in range(int((cov_extra + cov_per_ploidy)/2), int(3 * cov_per_ploidy / 2)):
            if cov in spectrum:
                count += spectrum[cov]
        mixture_weights.append(count / total_count)
        # other components
        for i in range(2, max_ploidy + 1):
            count = 1
            for cov in range(int((i - 0.5) * cov_per_ploidy), int((i + 0.5) * cov_per_ploidy)):
                if cov in spectrum:
                    count += spectrum[cov]
            mixture_weights.append(count / total_count)
        
        return CoverageDistribution(cov_per_ploidy, var_per_ploidy, err_scale, mixture_weights, cov_per_ploidy, cov_extra, var_extra, mixture_weight_extra)
    
    def __init__(self, cov_per_ploidy, var_per_ploidy, err_scale, mixture_weights, err_trunc_point, cov_extra, var_extra, mixture_weight_extra,
                 pseudo_cov = .25):
        self.cov_per_ploidy = cov_per_ploidy
        self.var_per_ploidy = var_per_ploidy
        self.err_scale = err_scale
        self.mixture_weights = mixture_weights
        self.err_trunc = err_trunc_point
        self.pseudo_cov = pseudo_cov
        # "_extra" shows the characterestics of the ultra-long component
        self.cov_extra = cov_extra
        self.var_extra = var_extra
        self.mixture_weight_extra = mixture_weight_extra 
    
    def pdf(self, x):
        return sum(self.likelihood(ploidy, x) for ploidy in range(-1, self.max_ploidy() + 1))
    
    def max_ploidy(self):
        return len(self.mixture_weights) - 1
    
    def cov_adj(self, cov):
        return cov + self.pseudo_cov
    
    def likelihood(self, ploidy, cov):
        l = 0.0
        if ploidy == -1: # for the extra component
            l = norm_pdf(cov, self.cov_extra, sqrt(self.var_extra))
            return self.mixture_weight_extra * l
        if ploidy == 0:
            # error distribution, switch to rate parameterization
            l = trunc_exp_pdf(cov, 1.0 / self.err_scale, self.err_trunc_point())
        else:
            l = norm_pdf(cov, self.cov_per_ploidy * ploidy, sqrt(self.var_per_ploidy * ploidy))
        return self.mixture_weights[ploidy] * l
    
    def err_trunc_point(self):
        return self.err_trunc
    
    def component_probabilities(self, cov):
        likelihoods = [self.likelihood(ploidy, cov) for ploidy in range(-1, self.max_ploidy() + 1)]
        Z = sum(likelihoods)
        if Z == 0:
            Z = 1
        for i in range(len(likelihoods)):
            likelihoods[i] /= Z
        return likelihoods
    
    def converged(self, other, tol):
        if other.cov_per_ploidy != 0.0 and abs(self.cov_per_ploidy / other.cov_per_ploidy - 1.0) > tol:
            return False
        if other.var_per_ploidy != 0.0 and abs(sqrt(self.var_per_ploidy / other.var_per_ploidy) - 1.0) > tol:
            return False
        if other.err_scale != 0.0 and abs(self.err_scale / other.err_scale - 1.0) > tol:
            return False
        max_ploidy_of_interest = 2
        for wt, other_wt, ploidy in zip(self.mixture_weights, other.mixture_weights, range(len(self.mixture_weights))):
            if ploidy > max_ploidy_of_interest:
                break
            if other_wt and abs(wt / other_wt - 1.0) > tol:
                return False
        if other.cov_extra != 0.0 and abs(self.cov_extra / other.cov_extra - 1.0) > tol:
            return False
        if other.var_extra != 0.0 and abs(self.var_extra / other.var_extra - 1.0) > tol:
            return False
        if other.mixture_weight_extra != 0.0 and abs(self.mixture_weight_extra / other.mixture_weight_extra - 1.0) > tol:
            return False
        return True
    
    
if __name__ == "__main__":
    
    from random import expovariate, normalvariate
    from collections import defaultdict
    import matplotlib.pyplot as plt
    
#    ###################################
#    # Generate some semi-realistic data
#    ###################################
#    
#    # the simulation parameters
#    N = 100000
#    var_inflation = 1.5
#    max_cov = 255
#    hap_cov = 25
#    err_frac = .01
#    hap_frac = .95
#    dip_frac = .01
#    high_copy_total_frac = .03
#    err_scale = 1.0
#    high_copy_fracs = [exp(-cov / hap_cov) for cov in range(3 * hap_cov, 2 * max_cov, hap_cov)]
#    Z = sum(high_copy_fracs)
#    for i in range(len(high_copy_fracs)):
#        high_copy_fracs[i] *= high_copy_total_frac / Z
#    mixture_fracs = [hap_frac, dip_frac] + high_copy_fracs
#    # do the simulation
#    counts = defaultdict(int)
#    for i in range(round(err_frac * N)):
#        x = int(expovariate(err_scale) + 1)
#        if x <= max_cov:
#            counts[x] += 1
#    for j in range(len(mixture_fracs)):
#        for i in range(round(mixture_fracs[j] * N)):
#            x = max(1, int(round(normalvariate((j + 1) * hap_cov, sqrt((j + 1) * var_inflation * hap_cov)))))
#            if x <= max_cov:
#                counts[x] += 1
#                
    
    ###################################
    # Load data
    ###################################
    
    counts_fp = "HG01358.paternal.ont.sorted.counts"
    counts = {}
    with open(counts_fp) as f:
        for line in f:
            cov, count = map(int, line.strip().split())
            counts[cov] = count
            
            
    covs = sorted(counts)[:200]
    plt.plot(covs, [counts[c] for c in covs])
        
    ###################################
    # Fit the model
    ###################################
    
    cov_model = CoverageDistribution.fit(counts, tol = 1e-2, max_iters=10)

    hap_cov = int(CoverageDistribution.heuristic_fit(counts).cov_per_ploidy)

    # look at how good the fit is:
    covs = sorted(counts)[:min(200, len(counts))]
    plt.plot(covs, [counts[f] for f in covs])
    Y = [cov_model.pdf(x) for x in covs]
    scale_factor = counts[hap_cov] / Y[covs.index(hap_cov)]
    Y = [scale_factor * y for y in Y]
    plt.plot(covs, Y)
    plt.show()
    
    ###################################
    # Query some probabilities
    ###################################
    
    for cov in range(0, 50, 5):
        prob_err = cov_model.probability_erroneous(cov)
        prob_correct = cov_model.probability_haploid(cov)
        prob_collapsed = cov_model.probability_collapsed(cov)
        
        print("coverage level {}: prob error {:.4f}, hap {:.4f}, coll {:.4f},".format(cov, prob_err, prob_correct, prob_collapsed))
