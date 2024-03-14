#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot
import numpy
import scipy.cluster
import scipy.spatial.distance
import sklearn.metrics
import sys
import typing
# custom lib
import mpllayout
import pylib


def read_oligo_list(fname: str = None) -> typing.Optional[list]:
	if fname is None:
		return None
	with pylib.util.get_fp(fname, "r") as fp:
		ret = fp.read().splitlines()
	return ret


def create_layout(heatmap_r):
	lc = mpllayout.LayoutCreator(
		left_margin=0.5,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=0.7,
	)

	colorbar = lc.add_frame("colorbar")
	colorbar.set_anchor("bottomleft")
	colorbar.set_size(1.5, 0.5)

	heatmap = lc.add_frame("heatmap")
	heatmap.set_anchor("bottomleft", ref_frame=colorbar,
		ref_anchor="topright", offsets=(0.2, 0.1))
	heatmap.set_size(6.0, 0.15 * heatmap_r)

	catebar = lc.add_frame("catebar")
	catebar.set_anchor("topleft", ref_frame=heatmap,
		ref_anchor="bottomleft", offsets=(0.0, -0.1))
	catebar.set_size(6.0, 0.3)

	dendro_top = lc.add_frame("dendro_top")
	dendro_top.set_anchor("bottomleft", ref_frame=heatmap,
		ref_anchor="topleft", offsets=(0, 0.1))
	dendro_top.set_anchor("topright", ref_frame=heatmap,
		ref_anchor="topright", offsets=(0, 1.6))

	dendro_right = lc.add_frame("dendro_right")
	dendro_right.set_anchor("bottomleft", ref_frame=heatmap,
		ref_anchor="bottomright", offsets=(0.1, 0))
	dendro_right.set_anchor("topright", ref_frame=heatmap,
		ref_anchor="topright", offsets=(1.6, 0))

	# create layout
	layout = lc.create_figure_layout()
	layout["dendro_top_i2d"] =\
		dendro_top.get_width() / dendro_top.get_height()
	layout["dendro_right_i2d"] =\
		dendro_right.get_height() / dendro_right.get_width()

	for n in ["catebar", "dendro_top", "dendro_right"]:
		ax = layout[n]
		for sp in ax.spines.values():
			sp.set_visible(False)
		ax.set_facecolor("#f0f0f8")

	layout["heatmap"].tick_params(
		left=False, labelleft=True,
		right=False, labelright=False,
		bottom=False, labelbottom=False,
		top=False, labeltop=False)
	layout["catebar"].tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=False, labelbottom=False,
		top=False, labeltop=False)
	layout["dendro_top"].tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=False, labelbottom=False,
		top=False, labeltop=False)
	layout["dendro_right"].tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=False, labelbottom=False,
		top=False, labeltop=False)

	return layout


def get_adjusted_dmax(dmax, n_leaves, i2d_ratio) -> float:
	return dmax / (1 - i2d_ratio / (2 * n_leaves))


def plot(png, table, *, oligo_list=None, cbar_label=None, dpi=300):
	if oligo_list is not None:
		table = table.filter_cols(by_tags=oligo_list)

	data = table.data.T
	nrow, ncol = data.shape

	# top hca and dendrogram
	top_dmat = sklearn.metrics.pairwise_distances(data.T, metric="euclidean")
	sq_top_dmat = scipy.spatial.distance.squareform(top_dmat,
		checks=False)
	top_linkage = scipy.cluster.hierarchy.linkage(sq_top_dmat,
		method="complete", metric="precomputed", optimal_ordering=True)
	top_dendro = scipy.cluster.hierarchy.dendrogram(top_linkage,
		orientation="top", no_plot=True)

	# right hca and dendrogram
	right_dmat = sklearn.metrics.pairwise_distances(data, metric="euclidean")
	sq_right_dmat = scipy.spatial.distance.squareform(right_dmat,
		checks=False)
	right_linkage = scipy.cluster.hierarchy.linkage(sq_right_dmat,
		method="complete", metric="precomputed", optimal_ordering=True)
	right_dendro = scipy.cluster.hierarchy.dendrogram(right_linkage,
		orientation="right", no_plot=True)

	# create layout
	layout = create_layout(nrow)
	figure = layout["figure"]

	# plot heatmap
	ax = layout["heatmap"]
	plot_data = data
	abs_max = numpy.abs(data).max()
	pcolor = ax.pcolor(plot_data[numpy.ix_(right_dendro["leaves"],
		top_dendro["leaves"])], cmap="bwr",
		vmin=-abs_max, vmax=+abs_max)
	ax.set_xlim(0, ncol)
	ax.set_ylim(0, nrow)
	ax.set_yticks(numpy.arange(nrow) + 0.5)
	ax.set_yticklabels([table.cols[i] for i in right_dendro["leaves"]],
		fontsize=10,
		horizontalalignment="right", verticalalignment="center")
	# ax.set_xlabel(label2, fontsize = 14, weight = "semibold")
	# ax.set_ylabel(label1, fontsize = 14, weight = "semibold")

	# colorbar
	ax = layout["colorbar"]
	cbar = figure.colorbar(pcolor, cax=ax, orientation="horizontal")
	cbar.outline.set_visible(False)
	cbar.set_label(cbar_label, fontsize=12)

	# catebar
	ax = layout["catebar"]
	for i, v in enumerate(top_dendro["leaves"]):
		sample = table.rows[v]
		wwtp = pylib.supp.SAMPLE_WWTP[sample]
		cate = pylib.supp.WWTP_CATE[wwtp]
		color = pylib.supp.CATE_COLOR[cate]
		p = matplotlib.patches.Rectangle((i, 0), 1, 1,
			edgecolor="none", facecolor=color)
		ax.add_patch(p)
	ax.set_xlim(0, ncol)
	ax.set_ylim(0, 1)
	# catebar legend
	handles = list()
	for k in pylib.supp.CATE_COLOR.def_key_order:
		color = pylib.supp.CATE_COLOR[k]
		p = matplotlib.patches.Rectangle((0, 0), 0, 0,
			edgecolor="none", facecolor=color, label=k)
		handles.append(p)
	ax.legend(handles=handles, loc=9, bbox_to_anchor=[0.50, -0.00], ncol=4,
		frameon=False, handlelength=0.75, fontsize=10, title="Community type")

	# top dendrogram
	ax = layout["dendro_top"]
	for xys in zip(top_dendro["icoord"], top_dendro["dcoord"]):
		line = matplotlib.lines.Line2D(*xys, linestyle="-", linewidth=1.0,
			color="#202080", zorder=3)
		ax.add_line(line)
	ax.set_xlim(0, 10 * data.shape[1])
	ax.set_ylim(0, get_adjusted_dmax(numpy.max(top_dendro["dcoord"]),
		ncol, i2d_ratio=layout["dendro_top_i2d"]))

	# right dendrogram
	ax = layout["dendro_right"]
	for xys in zip(right_dendro["dcoord"], right_dendro["icoord"]):
		line = matplotlib.lines.Line2D(*xys, linestyle="-", linewidth=1.0,
			color="#202080", zorder=3)
		ax.add_line(line)
	ax.set_xlim(0, get_adjusted_dmax(numpy.max(right_dendro["dcoord"]),
		nrow, i2d_ratio=layout["dendro_right_i2d"]))
	ax.set_ylim(0, 10 * data.shape[0])

	# savefig and clean-up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str,
		metavar="oligo-matrix",
		help="oligo relative abundance matrices")
	ap.add_argument("--oligo-list", "-l", type=str,
		metavar="txt",
		help="only plot the listed oligos in file, (default: plot all)")
	ap.add_argument("--plot", "-p", type=str, default="-",
		metavar="png",
		help="output plot image (default: stdout)")
	ap.add_argument("--colorbar-label", type=str,
		metavar="str",
		help="show this label on colorbar (default: no)")
	ap.add_argument("--dpi", type=int, default=600,
		metavar="int",
		help="output image dpi (default: 600)")
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer
	return args


def main():
	args = get_args()
	table = pylib.table.RelaAbundTable.from_file(args.input)
	# read_oligo_list show return None if args.read_oligo_list is None
	oligo_list = read_oligo_list(args.oligo_list)
	plot(args.plot, table, oligo_list=oligo_list, dpi=args.dpi,
		cbar_label=args.colorbar_label)
	return


if __name__ == "__main__":
	main()
