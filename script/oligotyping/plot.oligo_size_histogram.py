#!/usr/bin/env python3

import argparse
import io
import matplotlib
import matplotlib.pyplot
import numpy
import os
import sys

import mpllayout


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("oligo_output", type=str,
		help="oligotyping output directory")
	ap.add_argument("--plot", "-p", type=str, default="-",
		metavar="png",
		help="output plot image [stdout]")
	ap.add_argument("--dpi", type=int, default=300,
		metavar="int",
		help="output plot image dpi [300]")

	# parse and refine args
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer

	return args


class OligoSize(object):
	@classmethod
	def matrix_count_file_from_dir(cls, path: str) -> str:
		return os.path.join(path, "MATRIX-COUNT.txt")

	@property
	def count_data(self) -> numpy.ndarray:
		return self._count_data

	@count_data.setter
	def count_data(self, value):
		value = numpy.asarray(value, dtype=int)
		if value.ndim != 2:
			raise ValueError("count_data must be a 2-dimensional array")
		self._count_data = value
		return

	@property
	def total_counts(self) -> numpy.ndarray:
		return self.count_data.sum(axis=0)

	def __init__(self, count_data: numpy.ndarray, *ka, **kw):
		super().__init__(*ka, **kw)
		self.count_data = count_data
		return

	@classmethod
	def from_oligo_count_table(cls, path: str):
		if not os.path.isfile(path):
			raise IOError("path '%s' is not a file, make sure the oligotyping "
				"finished completely" % path)
		raw = numpy.loadtxt(path, dtype=object, delimiter="\t")
		data = raw[1:, 1:].astype(int)
		new = cls(count_data=data)
		return new

	@classmethod
	def from_oligo_output_dir(cls, path: str):
		if not os.path.isdir(path):
			raise IOError("path '%s' is not a directory" % path)
		fname = cls.matrix_count_file_from_dir(path)
		return cls.from_oligo_count_table(fname)

	def _plot_size_distribution_setup_layout(self) -> dict:
		lc = mpllayout.LayoutCreator(
			left_margin=0.7,
			right_margin=0.2,
			bottom_margin=0.7,
			top_margin=0.2,
		)

		ax = lc.add_frame("axes")
		ax.set_anchor("bottomleft")
		ax.set_size(6, 3)

		layout = lc.create_figure_layout()
		ax = layout["axes"]

		# apply axes styles
		for sp in ax.spines.values():
			sp.set_visible(False)
		ax.set_facecolor("#f8f8ff")
		ax.tick_params(
			left=True, labelleft=True,
			right=False, labelright=False,
			bottom=True, labelbottom=True,
			top=False, labeltop=False,
		)

		return layout

	def plot_size_distribution(self, png, *, dpi=300) -> None:
		layout = self._plot_size_distribution_setup_layout()
		figure = layout["figure"]

		total_counts = self.total_counts
		log_counts = numpy.log10(total_counts)
		log_10bins = numpy.log10(numpy.arange(1, 11))

		# xmax, in log-scale, ceiling
		xlog_max = int(numpy.ceil(log_counts.max()))

		# bin edges in log-scale
		log_bins = numpy.hstack(
			[log_10bins + i for i in range(xlog_max)] + [xlog_max]
		)

		# plot
		axes = layout["axes"]
		axes.hist(log_counts, bins=log_bins,
			align="mid",
			edgecolor="none",
			facecolor="#4040ff80"
		)

		# misc
		axes.set_xlim(0, xlog_max)
		axes.xaxis.set_major_formatter(
			matplotlib.ticker.FormatStrFormatter("10$^{%u}$")
		)
		axes.set_xticks(numpy.arange(xlog_max + 1, dtype=float))
		axes.set_xticks(log_bins, minor=True)
		for t in axes.xaxis.get_ticklabels():
			t.set_rotation(90)

		axes.set_xlabel("oligo size")
		axes.set_ylabel("count")


		# save figure and clean up
		figure.savefig(png, dpi=dpi)
		matplotlib.pyplot.close()
		return


def main():
	args = get_args()
	o = OligoSize.from_oligo_output_dir(args.oligo_output)
	o.plot_size_distribution(args.plot, dpi=args.dpi)
	return


if __name__ == "__main__":
	main()

