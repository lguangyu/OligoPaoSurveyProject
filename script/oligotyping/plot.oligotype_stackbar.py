#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import re
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type = str,
		help = "input oligotype percent table")
	ap.add_argument("-p", "--plot", type = str, default = "-",
		metavar = "png",
		help = "output plot image, or write to stdout by default '-'")
	# parse and refine arsg
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer
	return args


def load_data(fname):
	raw_text = numpy.loadtxt(fname, dtype = object, delimiter = "\t")
	oligotypes	= raw_text[0, 1:]
	samples		= raw_text[1:, 0]
	percents	= raw_text[1:, 1:].astype(float)
	return oligotypes, samples, percents


def setup_layout(figure, n_samples):
	layout = dict()
	layout["figure"] = figure

	# margins
	left_margin_inch	= 0.3
	right_margin_inch	= 2.3
	top_margin_inch		= 0.1
	bottom_margin_inch	= 1.5

	# bar width
	bar_width_inch		= 0.2
	bar_height_inch		= 3.0
	axes_width_inch		= bar_width_inch * n_samples
	axes_height_inch	= bar_height_inch

	# figure dimensions
	figure_width_inch	= left_margin_inch + axes_width_inch + right_margin_inch
	figure_height_inch	= top_margin_inch + axes_height_inch + bottom_margin_inch
	figure.set_size_inches(figure_width_inch, figure_height_inch)

	# create axes
	axes_left	= left_margin_inch / figure_width_inch
	axes_bottom	= bottom_margin_inch / figure_height_inch
	axes_width	= axes_width_inch / figure_width_inch
	axes_height	= axes_height_inch / figure_height_inch
	axes = figure.add_axes([axes_left, axes_bottom, axes_width, axes_height])
	layout["axes"] = axes
	# axes style
	#for sp in axes.spines.values():
	#	sp.set_visible(False)
	axes.set_facecolor("#E8E8F8")
	axes.tick_params(left = False, labelleft = False, bottom = False)

	return layout


def parse_sample_date(sample_name):
	m = re.search(r"(\d{6})", sample_name)
	return (m.group(1) if m is not None else sample_name)


def plot_oligotype_stack_bar(png, oligotypes, samples, percents):
	n_samples, n_oligotypes = len(samples), len(oligotypes)
	assert percents.shape == (n_samples, n_oligotypes)
	# create layout
	figure = matplotlib.pyplot.figure()
	layout = setup_layout(figure, n_samples)

	# sort sampels
	sample_argsort = numpy.argsort([parse_sample_date(s) for s in samples])

	# plot
	colors = matplotlib.pyplot.get_cmap("Set3").colors\
		+ matplotlib.pyplot.get_cmap("Pastel2").colors
	axes = layout["axes"]
	bottom = numpy.zeros(n_samples, dtype = float)
	handles = list()
	x = numpy.arange(n_samples) + 0.5
	for i in range(n_oligotypes):
		height = percents[sample_argsort, i]
		bar_plot = axes.bar(x, height, width = 0.8, bottom = bottom,
			align = "center", edgecolor = "#404040", linewidth = 0.5,
			facecolor = colors[i], label = oligotypes[i])
		bottom += height
		handles.append(bar_plot)

	# legend
	axes.legend(handles = handles, loc = 2, ncol = 2, fontsize = 8,
		bbox_to_anchor = (1.002, 1.02))

	# axis misc
	axes.set_xlim(0, n_samples)
	axes.set_ylim(0.0, 100.0)
	#axes.set_xlabel("Samples")
	axes.set_ylabel("Oligotypes")
	axes.set_xticks(x)
	axes.set_xticklabels(samples[sample_argsort], rotation = 90, fontsize = 8,
		horizontalalignment = "center", verticalalignment = "top")

	matplotlib.pyplot.savefig(png, dpi = 300)
	matplotlib.pyplot.close()
	return


def main():
	args = get_args()
	# load data
	oligotypes, samples, percents = load_data(args.input)
	# plot
	plot_oligotype_stack_bar(args.plot, oligotypes, samples, percents)
	return


if __name__ == "__main__":
	main()
