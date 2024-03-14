#!/usr/bin/env python3

import argparse
import collections
import io
import itertools
import matplotlib
import matplotlib.pyplot
import numpy
import os
import re
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("dirname", type=str,
		help="directory name that contains all blastdbcmd taxonomy tables")
	ag = ap.add_mutually_exclusive_group()
	ag.add_argument("-x", "--scan-ext", type=str, metavar="str",
		help="scan for all files with this extension in <dirname> to process "
			"(exclusive with -l/--file-list) [blastdbcmd]")
	ag.add_argument("-l", "--file-list", type=str, metavar="file",
		help="provide a list of files in <dirname> instead of scanning by "
			"extension (exclusive with -x/--scan-ext)")
	ap.add_argument("-d", "--delimiter", type=str, default="\t",
		metavar="str",
		help="delimiter in input/output tables [<tab>]")
	ap.add_argument("-k", "--key-field", type=int, default=2,
		metavar="int",
		help="the taxnomy value field number (0-based) in input files [0]")
	ap.add_argument("-n", "--num-legend-taxons", type=int, default=20,
		metavar="int",
		help="number of taxons to show in legend, increase this number too much"
			" will cause inevitable color confliction [20]")
	ap.add_argument("-t", "--table", type=str, default="-",
		metavar="file",
		help="output table [stdout]")
	ap.add_argument("-p", "--plot", type=str,
		metavar="png",
		help="draw a stackbar image plot only if an output file set")

	# parse and refine arsg
	args = ap.parse_args()
	if (args.scan_ext is None) and (args.file_list is None):
		args.scan_ext = "blastdbcmd"
	if args.table == "-":
		args.table = sys.stdout

	return args


def get_fp(f, *ka, factory=open, **kw) -> io.IOBase:
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument of get_fp() must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


def read_oligo_tax_count(fname, key_field, *, delimiter="\t") \
		-> collections.Counter:
	with open(fname, "r") as fp:
		tax = [i.split(delimiter)[key_field] for i in fp.read().splitlines()]
	return collections.Counter(tax)


def _iter_file_by_scan_ext(dirname, scan_ext: str) -> iter:
	for i in os.scandir(dirname):
		if i.name.endswith(scan_ext):
			yield i
	return


def _iter_file_by_file_list(dirname, file_list: str) -> iter:
	with open(file_list, "r") as fp:
		file_set = set(fp.read().splitlines())

	for i in os.scandir(dirname):
		if i.name in file_set:
			yield i
	return


def read_oligo_tax_count_in_dir(dirname, *, scan_ext=None, file_list=None,
		key_field=0, **kw) -> dict:
	if (scan_ext is None) and (file_list is None):
		raise ValueError("must provide either scan_ext or file_list")
	if (scan_ext is not None) and (file_list is not None):
		raise ValueError("cannot provide both scan_ext and file_list")
	if scan_ext is not None:
		file_iter = _iter_file_by_scan_ext(dirname, scan_ext)
	else:
		file_iter = _iter_file_by_file_list(dirname, file_list)

	ret = dict()
	for i in file_iter:
		oligo = int(re.search(r"^(\d+)", i.name).group(1))
		ret[oligo] = read_oligo_tax_count(i.path, key_field, **kw)
	return ret


def write_output_table(f, oligo_tax: dict, *, delimiter="\t") -> None:
	# get unique taxons from Counter's keys
	# get sorted as in list
	oligos = sorted(oligo_tax.keys())
	uniq_tax = sorted(set(itertools.chain(*oligo_tax.values())))
	total_counts = {k: v.total() for k, v in oligo_tax.items()}
	with get_fp(f, "w") as fp:
		header_line = delimiter.join([""] + ["OPU_%03u" % i for i in oligos])
		print(header_line, file=fp)
		for t in uniq_tax:
			vals = [oligo_tax[i][t] / (total_counts[i] or 1) for i in oligos]
			line = delimiter.join([t] + [str(i) for i in vals])
			print(line, file=fp)
	return


def setup_layout(figure, n_oligo):
	layout = dict()
	layout["figure"] = figure

	# margins
	left_margin_inch = 0.3
	right_margin_inch = 2.3
	top_margin_inch = 0.1
	bottom_margin_inch = 1.5

	# bar width
	bar_width_inch = 0.2
	bar_height_inch = 3.0
	axes_width_inch = bar_width_inch * n_oligo
	axes_height_inch = bar_height_inch

	# figure dimensions
	figure_width_inch = left_margin_inch + axes_width_inch + right_margin_inch
	figure_height_inch = top_margin_inch + axes_height_inch + bottom_margin_inch
	figure.set_size_inches(figure_width_inch, figure_height_inch)

	# create axes
	axes_left = left_margin_inch / figure_width_inch
	axes_bottom = bottom_margin_inch / figure_height_inch
	axes_width = axes_width_inch / figure_width_inch
	axes_height = axes_height_inch / figure_height_inch
	axes = figure.add_axes([axes_left, axes_bottom, axes_width, axes_height])
	layout["axes"] = axes
	# axes style
	# for sp in axes.spines.values():
	# sp.set_visible(False)
	axes.set_facecolor("#E8E8F8")
	axes.tick_params(left=False, labelleft=False, bottom=False)

	return layout


def get_oligo_tax_bootstrap_matrix(oligo_tax: dict):
	# get unique oligos list
	u_oligo = sorted(oligo_tax.keys())

	# get unique taxons list
	tax_iter = itertools.chain(*[v.keys() for v in oligo_tax.values()])
	u_tax = sorted(set(tax_iter))

	# get oligo tax bootstrap matrix
	bs_mat = numpy.zeros((len(u_oligo), len(u_tax)), dtype=float)
	for i, o in enumerate(u_oligo):
		bs_total = oligo_tax[o].total()
		for j, t in enumerate(u_tax):
			bs_mat[i, j] = oligo_tax[o][t] / (bs_total or 1.0)

	return u_oligo, u_tax, bs_mat


class LegendColors(list):
	@staticmethod
	def get_default_color_pool() -> list:
		colors = matplotlib.pyplot.get_cmap("Set3").colors \
			+ matplotlib.pyplot.get_cmap("Pastel2").colors \
			+ matplotlib.pyplot.get_cmap("Dark2").colors
		colors = [matplotlib.colors.to_hex(i) for i in colors]
		colors = LegendColors._drop_replicates(colors)
		return colors

	@staticmethod
	def _drop_replicates(vec: list):
		seen = set()
		ret = list()
		for i in vec:
			if i not in seen:
				ret.append(i)
				seen.add(i)
		return ret

	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		if not self:
			self.extend(self.get_default_color_pool())
		self.check_elements()
		return

	def check_elements(self) -> None:
		for i in super().__iter__():
			if not matplotlib.colors.is_color_like(i):
				raise ValueError("'%s' cannot be interpreted as color")
		return

	def __getitem__(self, index):
		"""
		allow cyclic index beyond the list length
		"""
		return super().__getitem__(index % len(self))

	def __iter__(self):
		return itertools.cycle(super().__iter__())


def plot_oligo_tax_stackbar(png, oligo_tax: dict, num_legend_taxons: int = 20):
	if png is None:
		return

	u_oligo, u_tax, bs_mat = get_oligo_tax_bootstrap_matrix(oligo_tax)

	n_oligo, n_tax = len(u_oligo), len(u_tax)
	assert bs_mat.shape == (n_oligo, n_tax)

	# create layout
	figure = matplotlib.pyplot.figure()
	layout = setup_layout(figure, n_oligo)

	# plot
	colors = LegendColors()
	axes = layout["axes"]
	bottom = numpy.zeros(n_oligo, dtype=float)
	handles = list()
	x = numpy.arange(n_oligo) + 0.5
	for h, t, c in zip(bs_mat.T, u_tax, colors):
		if len(handles) >= num_legend_taxons:
			break
		bar = axes.bar(x, h, width=0.8, bottom=bottom,
			align="center", edgecolor="#404040", linewidth=0.5,
			facecolor=c, label=t)
		bottom += h
		handles.append(bar)

	# check if 'others' needs to be added
	if len(handles) < n_tax:
		h = bs_mat[:, len(handles):].sum(axis=1)
		bar = axes.bar(x, h, width=0.8, bottom=bottom,
			align="center", edgecolor="#404040", linewidth=0.5,
			facecolor="#ffffff", label="others")
		bottom += h
		handles.append(bar)

	# legend
	axes.legend(handles=handles, loc=2, fontsize=8,
		bbox_to_anchor=(1.002, 1.02))

	# axis misc
	axes.set_xlim(0, n_oligo)
	axes.set_ylim(0.0, 1.0)
	axes.set_ylabel("Taxonomy bootstrap")
	axes.set_xticks(x)
	axes.set_xticklabels(["OLIGO_%03u" % i for i in u_oligo],
		rotation=90, fontsize=8,
		horizontalalignment="center", verticalalignment="top")

	matplotlib.pyplot.savefig(png, dpi=300)
	matplotlib.pyplot.close()
	return


def main():
	args = get_args()
	# load data
	oligo_tax = read_oligo_tax_count_in_dir(args.dirname,
		scan_ext=args.scan_ext, file_list=args.file_list,
		key_field=args.key_field, delimiter=args.delimiter)
	# table output
	write_output_table(args.table, oligo_tax, delimiter=args.delimiter)
	# plot output
	plot_oligo_tax_stackbar(args.plot, oligo_tax,
		num_legend_taxons=args.num_legend_taxons)
	return


if __name__ == "__main__":
	main()
