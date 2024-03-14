#!/usr/bin/env python3

import argparse
import collections
import json
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import matplotlib.ticker
import numpy
import sys
import typing
# custom lib
import mpllayout
import pylib


def load_json(f):
	with pylib.util.get_fp(f, "r") as fp:
		ret = json.load(fp)
	return ret


def create_layout():
	lc = mpllayout.LayoutCreator(
		left_margin=0.8,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=1.0,
	)

	ax = lc.add_frame("axes")
	ax.set_anchor("bottomleft")
	ax.set_size(5.0, 2.5)

	# create layout
	layout = lc.create_figure_layout()

	ax = layout["axes"]
	for sp in ax.spines.values():
		sp.set_visible(True)
	ax.set_facecolor("#ffffff")
	ax.tick_params(
		left=True, labelleft=True,
		right=False, labelright=False,
		bottom=False, labelbottom=True,
		top=False, labeltop=False)

	return layout


def group_rows_by_wwtp(table) -> dict:
	ret = collections.defaultdict(lambda: collections.defaultdict(list))
	for i, r in enumerate(table.rows):
		wwtp = pylib.supp.SAMPLE_WWTP[r]
		cate = pylib.supp.WWTP_CATE[wwtp]
		ret[cate][wwtp].append(i)
	return ret


def plot(cfg, table, otu_list):
	# create layout
	layout = create_layout()
	figure = layout["figure"]

	# group data by wwtp
	cate_group = group_rows_by_wwtp(table)

	# plot
	ax = layout["axes"]
	x = 0
	xpad = 0.1
	xw = (1 - xpad * 2) / len(otu_list)
	wwtp_full_list = list()
	cate_list = sorted(cate_group.keys())
	for cate in cate_list:
		wwtp_group = cate_group[cate]
		wwtp_list = sorted(wwtp_group.keys())
		wwtp_full_list.extend(wwtp_list)
		for wwtp in wwtp_list:
			for i, otu in enumerate(otu_list):
				color = otu["color"]
				d = table.filter_rows(by_index=wwtp_group[wwtp])\
					.filter_cols(by_tags=[otu["key"]])\
					.data.ravel()
				p = ax.boxplot([d], positions=[x + xpad + xw * (i + 0.5)],
					widths=[xw * 0.8], patch_artist=True, vert=True, zorder=3,
					boxprops={
						"edgecolor": color,
						"facecolor": color + "20",
					},
					flierprops={
						"marker": "x",
						"markersize": 5,
						"markeredgecolor": color,
					},
					medianprops={
						"color": color,
					},
					whiskerprops={
						"color": color,
					},
					capprops={
						"color": color,
					},
				)

			# wwtp misc
			p = matplotlib.patches.Rectangle((x, 0), 1, 1,
				facecolor=pylib.supp.CATE_COLOR[cate] + "40",
				edgecolor="none",
				zorder=1)
			ax.add_patch(p)
			ax.axvline(x, linewidth=0.5, color="#ffffff", zorder=2)

			x += 1

	# misc
	ax.set_xlim(0, x)
	ax.set_ylim(0, 0.25)
	ax.set_xticks(numpy.arange(x) + 0.5)
	ax.set_xticklabels([pylib.supp.WWTP_DISPLAY[i] for i in wwtp_full_list],
		rotation=90)
	for t, wwtp in zip(ax.get_xticklabels(), wwtp_full_list):
		cate = pylib.supp.WWTP_CATE[wwtp]
		t.set_color(pylib.supp.CATE_COLOR[cate])
	ax.yaxis.set_major_formatter(
		matplotlib.ticker.PercentFormatter(xmax=1, decimals=0)
	)
	ax.set_ylabel("relative abundance")

	# legend
	handles = list()
	for otu in otu_list:
		p = matplotlib.patches.Rectangle((0, 0), 0, 0, edgecolor=otu["color"],
			facecolor=otu["color"] + "20", label=otu["key"])
		handles.append(p)
	ax.legend(handles=handles, loc=1, bbox_to_anchor=[0.98, 0.98],
		frameon=False, handlelength=0.75)

	# savefig and clean-up
	png = cfg.get("file", sys.stdout.buffer)
	figure.savefig(png, dpi=cfg.get("dpi", 300))
	matplotlib.pyplot.close()
	return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("config", type=str, nargs="?", default="-",
		help="plot config file in json format")
	args = ap.parse_args()
	if args.config == "-":
		args.config = sys.stdin
	return args


def main():
	args = get_args()
	cfg = load_json(args.config)
	table = pylib.table.RelaAbundTable.from_file(cfg["data"]["file"])
	plot(cfg["figure"], table, otu_list=cfg.get("otu_list", None))
	return


if __name__ == "__main__":
	main()
