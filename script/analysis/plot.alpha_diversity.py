#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import numpy
import sys
# custom lib
import mpllayout
import pylib


def create_layout(ncate):
	lc = mpllayout.LayoutCreator(
		left_margin=0.8,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=0.7,
	)

	violin = lc.add_frame("violin")
	violin.set_anchor("bottomleft")
	violin.set_size(1.2 * ncate, 2.0)

	# create layout
	layout = lc.create_figure_layout()

	for n in ["violin"]:
		ax = layout[n]
		for sp in ax.spines.values():
			sp.set_visible(True)
		# ax.set_facecolor("#f0f0f8")
		ax.set_facecolor("#ffffff")
		ax.tick_params(
			left=True, labelleft=True,
			right=False, labelright=False,
			bottom=False, labelbottom=True,
			top=False, labeltop=False)
		ax.grid(linestyle="-", linewidth=1.0, color="#d0d0d0", zorder=1)

	return layout


def plot(png, table, *, dpi=300):
	ncate = len(pylib.supp.CATE_COLOR)

	# create layout
	layout = create_layout(ncate)
	figure = layout["figure"]

	# alpha diversity violin plot
	ax = layout["violin"]
	for i, c in enumerate(pylib.supp.CATE_COLOR.def_key_order):
		color = pylib.supp.CATE_COLOR[c]
		# filter alpha diversity data
		mask = [
			pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[v]] == c
			for v in table.rows
		]
		d = table.data[mask].squeeze()
		violin = ax.violinplot([d], positions=[i + 0.5], vert=True, widths=0.8)
		for patch in violin["bodies"]:
			patch.set(facecolor=color + "40", edgecolor=color, alpha=None,
				zorder=2)
		for n in ["cmaxes", "cbars", "cmins"]:
			violin[n].set(linewidth=1.0, color=color, zorder=3)

	ax.set_xlim(0, ncate)
	ax.set_xticks(numpy.arange(ncate) + 0.5)
	ax.set_xticklabels(pylib.supp.CATE_COLOR.def_key_order, fontsize=12)
	ax.set_ylabel("Oligo alpha diversity\n(Shannon)", fontsize=12)

	# savefig and clean-up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str,
		metavar="alpha-diversity-table",
		help="oligo alpha diversity table")
	ap.add_argument("--plot", "-p", type=str, default="-",
		metavar="png",
		help="output plot image (default: stdout)")
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
	plot(args.plot, table, dpi=args.dpi)
	return


if __name__ == "__main__":
	main()
