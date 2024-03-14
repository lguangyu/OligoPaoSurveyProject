#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot
import matplotlib.patches
import numpy
import scipy.spatial
import sklearn.decomposition
import sys
# custom lib
import mpllayout
import pylib


def create_layout():
	lc = mpllayout.LayoutCreator(
		left_margin=1.0,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=0.7,
	)

	ax = lc.add_frame("pca")
	ax.set_anchor("bottomleft")
	ax.set_size(5, 5)

	# create layout
	layout = lc.create_figure_layout()

	ax = layout["pca"]
	for sp in ax.spines.values():
		sp.set_visible(True)
	# ax.set_facecolor("#f0f0f8")
	ax.tick_params(
		left=False, labelleft=True,
		right=False, labelright=False,
		bottom=False, labelbottom=True,
		top=False, labeltop=False)
	ax.grid(linestyle="-", linewidth=1.0, color="#d0d0d0", zorder=1)

	return layout


def add_wwtp_hull(axes, table, trans_pos):
	for wwtp, cate in pylib.supp.WWTP_CATE.items():
		color = pylib.supp.CATE_COLOR[cate]

		# filter positions
		pos = trans_pos.compress([pylib.supp.SAMPLE_WWTP[i] == wwtp
			for i in table.rows], axis=0
		)

		if len(pos) == 0:
			continue

		# use convex hull instead of elipse
		hull = scipy.spatial.ConvexHull(pos)
		hull_xy = pos[hull.vertices]
		patch = matplotlib.patches.Polygon(hull_xy, closed=True, alpha=0.5,
			linestyle="--", facecolor=color + "40", edgecolor=color, zorder=2)
		axes.add_patch(patch)

		# add text
		# xy = (hull_xy.max(axis=0) + hull_xy.min(axis=0)) / 2
		xy = numpy.median(pos, axis=0)
		axes.text(*xy, pylib.supp.WWTP_DISPLAY[wwtp],
			fontsize=10, color="#000000",
			horizontalalignment="center", verticalalignment="center")
	return


def plot(png, table, method, *, dpi=300):
	# decomposition
	m = pylib.decomposition.DECOMP_METHODS.get(method)
	pos = m.fit_transform(table.data)

	# create layout
	layout = create_layout()

	# plot pca
	edgecolors = [
		pylib.supp.CATE_COLOR[pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[i]]]
		for i in table.rows
	]
	ax = layout["pca"]
	ax.scatter(pos[:, 0], pos[:, 1], s=40,  # c = "#4040ff",
		marker="o", linewidth=1.5, facecolors="#ffffff80",
		edgecolors=edgecolors, zorder=3)

	# add ellipses
	add_wwtp_hull(ax, table, pos)

	# label
	# for s, xy, c in zip(table.rows, pos, edgecolors):
	# ax.text(*xy, s, fontsize = 6, color = c + "40")

	# misc
	ax.set_xlabel(m.xlabel_str, fontsize=12, fontweight="semibold")
	ax.set_ylabel(m.ylabel_str, fontsize=12, fontweight="semibold")

	# save fig and clean up
	layout["figure"].savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		metavar="oligo-abund",
		help="oligo relative abundance matrix (default: stdin)")
	ap.add_argument("--method", "-m", type=str,
		default=pylib.decomposition.DECOMP_METHODS.default_key,
		choices=pylib.decomposition.DECOMP_METHODS.get_key_list(),
		help="decomposition method (default: %s)"
			% pylib.decomposition.DECOMP_METHODS.default_key)
	ap.add_argument("--plot", "-p", type=str, default="-",
		metavar="png",
		help="output plot image (default: stdout)")
	ap.add_argument("--dpi", type=int, default=600,
		metavar="int",
		help="output image dpi (default: 600)")
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.plot == "-":
		args.plot = sys.stdout.buffer
	return args


def main():
	args = get_args()
	table = pylib.table.RelaAbundTable.from_file(args.input)
	plot(args.plot, table, args.method, dpi=args.dpi)
	return


if __name__ == "__main__":
	main()
