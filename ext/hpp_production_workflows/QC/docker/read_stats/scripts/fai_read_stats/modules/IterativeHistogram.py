import math
import numpy
import sys


class IterativeHistogram:
    def __init__(self, start, stop, n_bins, unbounded_upper_bin=False, unbounded_lower_bin=False, include_upper_edge=True):
        self.start = start
        self.stop = stop
        self.n_bins = n_bins
        self.histogram = numpy.zeros(n_bins)
        self.bin_size = (stop - start)/n_bins
        self.edges = [start + self.bin_size*i for i in range(n_bins+1)]

        self.unbounded_upper_bin = unbounded_upper_bin
        self.unbounded_lower_bin = unbounded_lower_bin
        self.include_upper_edge = include_upper_edge

    def get_bin(self, x):
        # find index of bin by normalizing and centering the value w.r.t. bin edges and add value to that bin
        bin_index = int(math.floor((x - self.start)/self.bin_size))

        if x == self.stop and self.include_upper_edge:
            bin_index = self.n_bins - 1

        if self.unbounded_lower_bin and x < self.start:
            bin_index = 0

        if self.unbounded_upper_bin and x > self.stop:
            bin_index = self.n_bins - 1

        return bin_index

    def update(self, x):
        bin_index = self.get_bin(x)

        if 0 <= bin_index <= (self.n_bins-1):
            self.histogram[bin_index] += 1

        # print(x, bin_index)
        # print(self.edges)
        # print(self.histogram)

    def get_histogram(self):
        return self.histogram

    def get_normalized_histogram(self):
        total = sum(self.histogram)
        normalized_histogram = self.histogram/numpy.sum(self.histogram)

        return normalized_histogram


if __name__ == "__main__":
    # test the iterative histogram
    iterative_histogram = IterativeHistogram(start=0, stop=10, n_bins=10)
    iterative_histogram.update(0)        # 1
    iterative_histogram.update(-1)       # None
    iterative_histogram.update(10)       # 10
    iterative_histogram.update(9.99999)  # 10
    iterative_histogram.update(10.0001)  # None
    iterative_histogram.update(0.5)      # 1
    iterative_histogram.update(1.5)      # 2
    iterative_histogram.update(1.0)      # 2
    iterative_histogram.update(1.99999)  # 2

    # ^ expect [2,3,0,0,0,0,0,0,0,2]

    print(iterative_histogram.get_histogram())

    iterative_histogram = IterativeHistogram(start=0, stop=1.0, n_bins=10)
    iterative_histogram.update(0)         # 1
    iterative_histogram.update(-0.1)      # None
    iterative_histogram.update(1.0)       # 10
    iterative_histogram.update(0.999999)  # 10
    iterative_histogram.update(1.00001)   # None
    iterative_histogram.update(0.05)      # 1
    iterative_histogram.update(0.15)      # 2
    iterative_histogram.update(0.10)      # 2
    iterative_histogram.update(0.199999)  # 2

    # ^ expect [2,3,0,0,0,0,0,0,0,2]

    print(iterative_histogram.get_histogram())

    iterative_histogram = IterativeHistogram(start=1, stop=2.0, n_bins=10)
    iterative_histogram.update(1 + 0)         # 1
    iterative_histogram.update(1 + -0.1)      # None
    iterative_histogram.update(1 + 1.0)       # 10
    iterative_histogram.update(1 + 0.999999)  # 10
    iterative_histogram.update(1 + 1.00001)   # None
    iterative_histogram.update(1 + 0.05)      # 1
    iterative_histogram.update(1 + 0.15)      # 2
    iterative_histogram.update(1 + 0.10)      # 2
    iterative_histogram.update(1 + 0.199999)  # 2

    # ^ expect [2,3,0,0,0,0,0,0,0,2]

    print(iterative_histogram.get_histogram())

    iterative_histogram = IterativeHistogram(start=-0.5, stop=0.5, n_bins=10)
    iterative_histogram.update(-0.5 + 0)         # 1
    iterative_histogram.update(-0.5 + -0.1)      # None
    iterative_histogram.update(-0.5 + 1.0)       # 10 right edge
    iterative_histogram.update(-0.5 + 0.999999)  # 10
    iterative_histogram.update(-0.5 + 1.00001)   # None
    iterative_histogram.update(-0.5 + 0.05)      # 1
    iterative_histogram.update(-0.5 + 0.15)      # 2
    iterative_histogram.update(-0.5 + 0.10)      # 2 ... in near-edge cases float division may shift left a bin
    iterative_histogram.update(-0.5 + 0.199999)  # 2

    # DON'T USE THIS CLASS FOR BINS WITH SIZE that NEARS FLOAT PRECISION
    # ^ expect [2,3,0,0,0,0,0,0,0,2]

    print(iterative_histogram.get_histogram())
    print(iterative_histogram.get_normalized_histogram())
    print(sum(iterative_histogram.get_normalized_histogram()))

    iterative_histogram = IterativeHistogram(start=0, stop=1.0, n_bins=10, unbounded_lower_bin=True, unbounded_upper_bin=True)
    iterative_histogram.update(0)         # 1
    iterative_histogram.update(-0.1)      # 1
    iterative_histogram.update(1.0)       # 10
    iterative_histogram.update(0.999999)  # 10
    iterative_histogram.update(1.00001)   # 10
    iterative_histogram.update(0.05)      # 1
    iterative_histogram.update(0.15)      # 2
    iterative_histogram.update(0.10)      # 2
    iterative_histogram.update(0.199999)  # 2

    # ^ expect [3,3,0,0,0,0,0,0,0,3]

    print(iterative_histogram.get_histogram())
    print(iterative_histogram.get_normalized_histogram())
    print(sum(iterative_histogram.get_normalized_histogram()))


