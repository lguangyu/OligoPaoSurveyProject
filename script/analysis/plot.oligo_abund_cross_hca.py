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


def create_layout(heatmap_r, heatmap_c):
	lc = mpllayout.LayoutCreator(
		left_margin=2.0,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=2.0,
	)

	heatmap = lc.add_frame("heatmap")
	heatmap.set_anchor("bottomleft")
	heatmap.set_size(0.25 * heatmap_c, 0.25 * heatmap_r)

	colorbar = lc.add_frame("colorbar")
	colorbar.set_anchor("topright", ref_frame=heatmap,
		ref_anchor="bottomleft", offsets=(-0.2, -0.2))
	colorbar.set_size(1.5, 0.5)

	dendro_top = lc.add_frame("dendro_top")
	dendro_top.set_anchor("bottomleft", ref_frame=heatmap,
		ref_anchor="topleft", offsets=(0, 0.1))
	dendro_top.set_anchor("topright", ref_frame=heatmap,
		ref_anchor="topright", offsets=(0, 2.1))

	dendro_right = lc.add_frame("dendro_right")
	dendro_right.set_anchor("bottomleft", ref_frame=heatmap,
		ref_anchor="bottomright", offsets=(0.1, 0))
	dendro_right.set_anchor("topright", ref_frame=heatmap,
		ref_anchor="topright", offsets=(2.1, 0))

	# create layout
	layout = lc.create_figure_layout()
	layout["dendro_top_i2d"] =\
		dendro_top.get_width() / dendro_top.get_height()
	layout["dendro_right_i2d"] =\
		dendro_right.get_height() / dendro_right.get_width()

	for n in ["heatmap", "dendro_top", "dendro_right"]:
		ax = layout[n]
		for sp in ax.spines.values():
			sp.set_visible(False)
		ax.set_facecolor("#f0f0f8")

	layout["heatmap"].tick_params(
		left=False, labelleft=True,
		right=False, labelright=False,
		bottom=False, labelbottom=True,
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


def plot(png, *, metric="cosine", table1, table2, label1="", label2="",
		list1=None, list2=None, dpi=300):
	if list1 is not None:
		table1 = table1.filter_cols(by_tags=list1)
	if list2 is not None:
		table2 = table2.filter_cols(by_tags=list2)

	# makesure the heatmap always width>height
	if table1.n_cols > table2.n_cols:
		table1, table2 = table2, table1
		label1, label2 = label2, label1

	# calculate heatmap matrix
	metric_meth = pylib.corr_dist.CORR_DIST_FUNC.get(metric)
	data = metric_meth(table1.data.T, table2.data.T)

	# table1->table2 hca and dendrogram (top)
	top_dmat = sklearn.metrics.pairwise_distances(data.T, metric="euclidean")
	sq_top_dmat = scipy.spatial.distance.squareform(top_dmat,
		checks=False)
	top_linkage = scipy.cluster.hierarchy.linkage(sq_top_dmat,
		method="complete", metric="precomputed", optimal_ordering=True)
	top_dendro = scipy.cluster.hierarchy.dendrogram(top_linkage,
		orientation="top", no_plot=True)

	# table2->table1 hca and dendrogram (right)
	right_dmat = sklearn.metrics.pairwise_distances(data, metric="euclidean")
	sq_right_dmat = scipy.spatial.distance.squareform(right_dmat,
		checks=False)
	right_linkage = scipy.cluster.hierarchy.linkage(sq_right_dmat,
		method="complete", metric="precomputed", optimal_ordering=True)
	right_dendro = scipy.cluster.hierarchy.dendrogram(right_linkage,
		orientation="right", no_plot=True)

	# create layout
	layout = create_layout(*data.shape)
	figure = layout["figure"]

	# plot heatmap
	ax = layout["heatmap"]
	plot_data = metric_meth.to_plot_data(data)
	pcolor = ax.pcolor(plot_data[numpy.ix_(right_dendro["leaves"],
		top_dendro["leaves"])], cmap=metric_meth.cmap,
		vmin=metric_meth.vmin, vmax=metric_meth.vmax)
	ax.set_xlim(0, data.shape[1])
	ax.set_ylim(0, data.shape[0])
	ax.set_xticks(numpy.arange(data.shape[1]) + 0.5)
	ax.set_xticklabels(table2.cols[top_dendro["leaves"]],
		fontsize=10, rotation=90,
		horizontalalignment="center", verticalalignment="top")
	ax.set_yticks(numpy.arange(data.shape[0]) + 0.5)
	ax.set_yticklabels(table1.cols[right_dendro["leaves"]],
		fontsize=10,
		horizontalalignment="right", verticalalignment="center")
	ax.set_xlabel(label2, fontsize=14, weight="semibold")
	ax.set_ylabel(label1, fontsize=14, weight="semibold")

	# colorbar
	ax = layout["colorbar"]
	cbar = figure.colorbar(pcolor, cax=ax, orientation="horizontal")
	cbar.outline.set_visible(False)
	cbar.set_label(metric_meth.name_str, fontsize=12)

	# top dendrogram
	ax = layout["dendro_top"]
	for xys in zip(top_dendro["icoord"], top_dendro["dcoord"]):
		line = matplotlib.lines.Line2D(*xys, linestyle="-", linewidth=1.0,
			color="#202080", zorder=3)
		ax.add_line(line)
	ax.set_xlim(0, 10 * data.shape[1])
	ax.set_ylim(0, get_adjusted_dmax(numpy.max(top_dendro["dcoord"]),
		data.shape[1], i2d_ratio=layout["dendro_top_i2d"]))

	# right dendrogram
	ax = layout["dendro_right"]
	for xys in zip(right_dendro["dcoord"], right_dendro["icoord"]):
		line = matplotlib.lines.Line2D(*xys, linestyle="-", linewidth=1.0,
			color="#202080", zorder=3)
		ax.add_line(line)
	ax.set_xlim(0, get_adjusted_dmax(numpy.max(right_dendro["dcoord"]),
		data.shape[0], i2d_ratio=layout["dendro_right_i2d"]))
	ax.set_ylim(0, 10 * data.shape[0])

	# savefig and clean-up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("--metric", "-m", type=str,
		default=pylib.corr_dist.CORR_DIST_FUNC.default_key,
		choices=pylib.corr_dist.CORR_DIST_FUNC.get_key_list(),
		help="metric to calculate correlation heatmap (default: %s)"
			% pylib.corr_dist.CORR_DIST_FUNC.default_key)

	ag = ap.add_argument_group("1st input")
	ag.add_argument("--oligo-table-1", "-1", type=str, required=True,
		metavar="oligo-abund",
		help="oligo relative abundance/compositional data table (required)")
	ag.add_argument("--oligo-label-1", type=str, required=True,
		metavar="str",
		help="label of the input table, will be shown aside axis (required)")
	ag.add_argument("--oligo-list-1", type=str,
		metavar="txt",
		help="only calculate correlation with and visualize the listed oligos "
			"(default: all oligos)")

	ag = ap.add_argument_group("2nd input")
	ag.add_argument("--oligo-table-2", "-2", type=str, required=True,
		metavar="oligo-abund",
		help="oligo relative abundance/compositional data table (required)")
	ag.add_argument("--oligo-label-2", type=str, required=True,
		metavar="str",
		help="label of the input table, will be shown aside axis (required)")
	ag.add_argument("--oligo-list-2", type=str,
		metavar="txt",
		help="only calculate correlation with and visualize the listed oligos "
			"(default: all oligos)")

	ag = ap.add_argument_group("output")
	ag.add_argument("--plot", "-p", type=str, default="-",
		metavar="png",
		help="output plot image (default: stdout)")
	ag.add_argument("--dpi", type=int, default=600,
		metavar="int",
		help="output image dpi (default: 600)")

	# parse and refine args
	args = ap.parse_args()
	if args.plot == "-":
		args.plot = sys.stdout.buffer

	return args


def main():
	args = get_args()
	# load data
	table1 = pylib.table.RelaAbundTable.from_file(args.oligo_table_1)
	label1 = args.oligo_label_1
	list1 = read_oligo_list(args.oligo_list_1)

	table2 = pylib.table.RelaAbundTable.from_file(args.oligo_table_2)
	label2 = args.oligo_label_2
	list2 = read_oligo_list(args.oligo_list_2)

	plot(args.plot, metric=args.metric, dpi=args.dpi,
		table1=table1, label1=label1, list1=list1,
		table2=table2, label2=label2, list2=list2,
	)
	return


if __name__ == "__main__":
	main()
