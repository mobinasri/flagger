from matplotlib import pyplot
import numpy
import sys
import os


def plot_nx(length_frequencies, total_length, output_dir):

    figure = pyplot.figure()
    axes = pyplot.axes()

    legend_names = list()

    x1 = 0
    y_prev = None

    x_coords = list()
    y_coords = list()

    for length,frequency in length_frequencies:
        for i in range(frequency):
            y = length
            width = float(length) / float(total_length)
            x2 = x1 + width

            if y_prev is not None:
                x_coords.extend([x1, x1])
                y_coords.extend([y_prev, y])

            x_coords.extend([x1, x2])
            y_coords.extend([y, y])

            x1 = x2
            y_prev = y

    if y_coords[-1] != 0:
        y_coords.append(0)
        x_coords.append(x_coords[-1])

    axes.plot(x_coords, y_coords, linewidth=0.6)

    axes.axvline(0.5, linestyle="--", alpha=0.3, linewidth=0.7, zorder=-1)

    axes.set_xlim([0, 1])

    axes.set_title("Nx")
    axes.set_xlabel("Cumulative coverage (normalized to 1)")
    axes.set_ylabel("Length")

    path = os.path.join(output_dir, "Nx.png")
    sys.stderr.write("SAVING FIGURE: %s\n" % path)
    figure.savefig(path, dpi=200)

    pyplot.close()


def plot_iterative_histogram(iterative_histogram, output_dir):
    figure = pyplot.figure()
    axes = pyplot.axes()

    bounds = numpy.array(iterative_histogram.edges)

    center = (bounds[:-1] + bounds[1:]) / 2

    axes.bar(center, iterative_histogram.histogram, width=iterative_histogram.bin_size, align="center")

    axes.set_xlabel("Read length (bp)")
    axes.set_ylabel("Frequency")

    path = os.path.join(output_dir, "histogram.png")
    sys.stderr.write("SAVING FIGURE: %s\n" % path)
    figure.savefig(path, dpi=200)

    pyplot.close()
