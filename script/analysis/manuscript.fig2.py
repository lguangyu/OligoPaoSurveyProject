#!/usr/bin/env python3

import numpy
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import matplotlib.patches
import matplotlib.patheffects
import mpllayout
import os
import scipy.cluster
import scipy.spatial
import scipy.stats
import skbio
import sklearn
import sys
# custom lib
import pylib


ORGN_PLOT_CFG = [
	{
		"key": "acc",
		"display": "Ca_Accumulibacter",
		"path_key": "candidatus_accumulibacter",
		"color": "#4040ff",
		"tsne_random_state": 675622,
	},
	{
		"key": "dech",
		"display": "Dechloromonas",
		"path_key": "dechloromonas",
		"color": "#f0691a",
		"tsne_random_state": 157224,
	},
	{
		"key": "tetr",
		"display": "Tetrasphaera",
		"path_key": "tetrasphaera",
		"color": "#0a6e16",
		"tsne_random_state": 822881,
	}
]


def load_list(fname):
	with pylib.util.get_fp(fname, "r") as fp:
		ret = fp.read().splitlines()
	return ret


def create_layout() -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.7,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=0.7,
	)

	tsne_size = 2.5
	dendro_dfrac = 0.3  # in relative to tsne_size
	hmap_hpadfrac = 0.15  # in relative to tsne_size
	cate_wfrac = 0.05  # in relative to tsne_size
	hmap_width = 3.0
	center_gap = 1.2

	for i, orgn_cfg in enumerate(ORGN_PLOT_CFG[::-1]):
		key = orgn_cfg["key"]
		# tsne axes
		tsne = lc.add_frame(key + "_tsne")
		if i == 0:
			tsne.set_anchor("bottomleft")
		else:
			prev_tsne = lc.get_frame(ORGN_PLOT_CFG[-i]["key"] + "_tsne")
			tsne.set_anchor("bottomleft", prev_tsne, "topleft",
				offsets=(0, 0.3))
		tsne.set_size(tsne_size, tsne_size)

		# top dendro axes
		tden = lc.add_frame(key + "_tden")
		tden.set_anchor("topleft", tsne, "topright", offsets=(center_gap, 0))
		tden.set_size(hmap_width, dendro_dfrac * tsne_size)

		# cate patch axes
		cate = lc.add_frame(key + "_cate")
		cate.set_anchor("topleft", tden, "bottomright")
		cate.set_anchor("bottomright", tsne, "bottomright",
			offsets=(center_gap + hmap_width + cate_wfrac * tsne_size,
				hmap_hpadfrac * tsne_size)
		)

		# right dendro axes
		rden = lc.add_frame(key + "_rden")
		rden.set_anchor("topleft", cate, "topright")
		rden.set_anchor("bottomright", cate, "bottomright",
			offsets=(dendro_dfrac * tsne_size, 0))

		# heatmap axes
		hmap = lc.add_frame(key + "_hmap")
		hmap.set_anchor("topleft", tden, "bottomleft")
		hmap.set_anchor("bottomright", cate, "bottomleft")

		# color bar axes
		if i == len(ORGN_PLOT_CFG) - 1:
			cbar = lc.add_frame("cbar")
			cbar.set_anchor("bottomleft", rden, "topleft",
				offsets=(0, dendro_dfrac * 0.2 * tsne_size))
			cbar.set_size(dendro_dfrac * tsne_size,
				dendro_dfrac * 0.2 * tsne_size)

	# create layout
	layout = lc.create_figure_layout()

	# apply axes styles
	for name, axes in layout.items():
		if name == "figure":
			continue

		axes.set_facecolor("#ffffff")

		if name.endswith("_tsne"):
			for sp in axes.spines.values():
				sp.set_visible(True)
			axes.tick_params(
				left=True, labelleft=True,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=True, labelbottom=True,
			)
		elif name.endswith("_tden") or name.endswith("_rden"):
			for sp in axes.spines.values():
				sp.set_visible(False)
			axes.tick_params(
				left=False, labelleft=False,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)
		elif name.endswith("_cate"):
			for sp in axes.spines.values():
				sp.set_visible(False)
			axes.tick_params(
				left=False, labelleft=False,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)
		elif name.endswith("_hmap"):
			for sp in axes.spines.values():
				sp.set_visible(True)
			axes.tick_params(
				left=False, labelleft=False,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)
		elif name == "cbar":
			for sp in axes.spines.values():
				sp.set_visible(False)
			axes.tick_params(
				direction="in",
				left=False, labelleft=False,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)

	return layout


def _adjust_color_brightness(c: str, multiplier: float) -> str:
	hsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
	hsv[2] = min(hsv[2] * multiplier, 1.0)
	return matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb(hsv))


def tsne_add_hull(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable,
		*, xys):
	for wwtp, cate in pylib.supp.WWTP_CATE.items():
		color = pylib.supp.CATE_COLOR[cate]

		# filter positions per wwtp
		wwtp_mask = [pylib.supp.SAMPLE_WWTP[i] == wwtp for i in data.rows]
		wwtp_xys = xys[wwtp_mask]
		if len(wwtp_xys) == 0:
			continue

		# use convex hull instead of elipse
		hull = scipy.spatial.ConvexHull(wwtp_xys)
		hull_xys = wwtp_xys[hull.vertices]
		patch = matplotlib.patches.Polygon(hull_xys, closed=True, alpha=0.5,
			linestyle="--", facecolor=color + "40", edgecolor=color + "80",
			zorder=2
		)
		axes.add_patch(patch)

		# add wwtp label text
		xy = numpy.median(wwtp_xys, axis=0)
		axes.text(*xy, pylib.supp.WWTP_DISPLAY_ABBR[wwtp],
			fontsize=12, color=_adjust_color_brightness(color, 0.75),
			horizontalalignment="center", verticalalignment="center",
			path_effects=[matplotlib.patheffects.withStroke(linewidth=3,
				foreground="#ffffff")],
		)

	return


def euclidean_distance_wrapper(x1, x2):
	return sklearn.metrics.pairwise.euclidean_distances(x1.reshape(1, -1),
		x2.reshape(1, -1))


def tsne_anosim(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable, *,
		permutations=100000):
	grouping = [pylib.supp.SAMPLE_WWTP[i] for i in data.rows]
	dist_mat = skbio.DistanceMatrix.from_iterable(data.data,
		euclidean_distance_wrapper, validate=False)
	anosim = skbio.stats.distance.anosim(dist_mat, grouping=grouping,
		permutations=permutations)

	# add anosim result as text
	pvalue = anosim["p-value"]
	text = ("\n").join([
		anosim["method name"],
		"%s=%.3f" % (anosim["test statistic name"], anosim["test statistic"]),
		"p-val%s" % ("<0.001" if pvalue < 1e-3 else "=%.3f" % pvalue),
	])

	color = "#ff0000" if anosim["p-value"] < 0.01 else "#000000"
	axes.text(1.02, 0.98, text, fontsize=10, color=color, clip_on=False,
		transform=axes.transAxes,
		horizontalalignment="left", verticalalignment="top",
	)

	return


def plot_tsne(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable, *,
		tsne_random_state=None, with_xlabel=True, with_ylabel=True,
		anosim_permutations=100000, title=None,
	):
	# decomposition
	tsne = pylib.decomposition.DECOMP_METHODS.get("t-sne",
		random_state=tsne_random_state,)
	xys = tsne.fit_transform(data.data)

	# plot scatter
	edgecolors = [
		pylib.supp.CATE_COLOR[pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[i]]]
		for i in data.rows
	]
	axes.scatter(xys[:, 0], xys[:, 1], s=20, marker="o",
		linewidth=1.0, facecolors="#ffffff80", edgecolors=edgecolors, zorder=3,
	)

	# add convex hull
	tsne_add_hull(axes, data, xys=xys)

	# add anosim results
	tsne_anosim(axes, data, permutations=anosim_permutations)

	# add title
	if title:
		axes.text(0.98, 0.98, title, fontsize=12, color="#000000", zorder=4,
			transform=axes.transAxes,
			horizontalalignment="right", verticalalignment="top",
			path_effects=[matplotlib.patheffects.withStroke(linewidth=3,
				foreground="#ffffff")],
		)

	# misc
	if with_xlabel:
		axes.set_xlabel(tsne.xlabel_str, fontsize=12, fontweight="semibold")
	if with_ylabel:
		axes.set_ylabel(tsne.ylabel_str, fontsize=12, fontweight="semibold")

	return


def hca_dendro(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable, *,
		orientation: str) -> numpy.ndarray:
	if orientation not in ["top", "right"]:
		raise ValueError("orientation must be either 'top' or 'right', "
			"got '%s'" % orientation)

	# hca
	d = data.data if orientation == "right" else data.data.T
	linkage = scipy.cluster.hierarchy.linkage(d, method="average",
		optimal_ordering=True)
	dendro = scipy.cluster.hierarchy.dendrogram(linkage,
		orientation=orientation, no_plot=True)
	leaves = numpy.asarray(dendro["leaves"])

	# plot dendrogram
	if orientation == "top":
		xys_iter = zip(dendro["icoord"], dendro["dcoord"])
	else:
		xys_iter = zip(dendro["dcoord"], dendro["icoord"])
	for x, y in xys_iter:
		axes.plot(x, y, clip_on=False, linestyle="-", linewidth=1.0,
			color="#200064")

	# misc
	imax = 10 * len(leaves)
	dmax = numpy.max(dendro["dcoord"])
	axes.set_xlim(0, imax if orientation == "top" else dmax)
	axes.set_ylim(0, dmax if orientation == "top" else imax)

	return leaves


def oligo_clr_ttest_with_fdr(data: pylib.table.RelaAbundTable,
		cate_grouping: list, sig_level: float = 0.05
	) -> numpy.ndarray:

	mask = list()
	for cate in cate_grouping:
		m = [pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[i]] == cate
			for i in data.rows
		]
		mask.append(m)

	ttest_res = scipy.stats.ttest_ind(
		data.filter_rows(by_index=mask[0]).data,
		data.filter_rows(by_index=mask[1]).data,
		axis=0, equal_var=False,
	)

	# fdr control
	pval = ttest_res.pvalue
	fdr_thres = numpy.linspace(0, sig_level, len(pval) + 1)[1:]
	argsort_idx = numpy.argsort(pval)
	sig_flag = numpy.zeros(len(pval), dtype=bool)
	sig_flag[argsort_idx] = (pval[argsort_idx] <= fdr_thres)

	return pval, sig_flag


def add_1v1_ttest_with_fdr(axes: matplotlib.axes.Axes,
		data: pylib.table.RelaAbundTable, *,
		sig_level: float = 0.05,
		rden_native=None, rden_leaves=None,
		tden_native=None, tden_leaves=None,
	):

	wwtp_list = numpy.asarray([pylib.supp.SAMPLE_WWTP[i] for i in data.rows],
		dtype=object)
	wwtp_unique = numpy.unique(wwtp_list)
	pval = numpy.empty((len(wwtp_unique), data.n_cols), dtype=float)

	# test for each otu, per wwtp at one vs rest basis
	for i, otu in enumerate(data.cols):
		for j, wwtp_j in enumerate(wwtp_unique):
			pval_temp = list()
			for k, wwtp_k in enumerate(wwtp_unique):
				if j == k:
					continue
				mask_j = (wwtp_list == wwtp_j)
				mask_k = (wwtp_list == wwtp_k)
				ttest_res = scipy.stats.ttest_ind(
					data.data[mask_j, i], data.data[mask_k, i],
					equal_var=False,
				)
				pval_temp.append(ttest_res.pvalue)
			pval[j, i] = max(pval_temp)

	# fdr control
	pval_ravel = pval.ravel()
	fdr_thres = numpy.linspace(0, sig_level, len(pval_ravel) + 1)[1:]
	argsort_idx = numpy.argsort(pval_ravel)
	sig_flag = numpy.zeros(len(pval_ravel), dtype=bool)
	sig_flag[argsort_idx] = (pval_ravel[argsort_idx] <= fdr_thres)
	sig_flag = sig_flag.reshape(pval.shape)

	# add boxed for significant wwtp-otu
	if tden_native is None:
		tden_native = data.cols
	if tden_leaves is None:
		tden_leaves = numpy.arange(data.n_cols)
	if rden_native is None:
		rden_native = wwtp_unique
	if rden_leaves is None:
		rden_leaves = numpy.arange(len(wwtp_unique))

	local_tmap = {v: i for i, v in enumerate(data.cols)}
	local_rmap = {v: i for i, v in enumerate(wwtp_unique)}

	xs, ys = list(), list()
	for x, c in enumerate(tden_leaves):
		for y, r in enumerate(rden_leaves):
			sig_flag_r = local_rmap[rden_native[r]]
			sig_flag_c = local_tmap[tden_native[c]]
			if sig_flag[sig_flag_r, sig_flag_c]:
				xs.append(x + 0.5)
				ys.append(y + 0.5)
				axes.scatter(xs, ys, s=20, marker="*", facecolors="#000000",
					linewidths=1.0, zorder=4
				)

	return


def plot_hmap(data: pylib.table.RelaAbundTable, *,
		hmap_axes: matplotlib.axes.Axes,
		tden_axes: matplotlib.axes.Axes,
		rden_axes: matplotlib.axes.Axes,
		cate_axes: matplotlib.axes.Axes,
		unified_vmin: float = None,
		unified_vmax: float = None,
		cbar_axes: matplotlib.axes.Axes = None,
		with_cate_legend: bool = None,
	):

	# t-test oligo abundances
	_, oligo_sig_flag = oligo_clr_ttest_with_fdr(data, ("EBPR", "S2EBPR"),
		sig_level=0.05,
	)

	# prepare data for plot
	plot_data = data.grouped_row_mean(grouping=[pylib.supp.SAMPLE_WWTP[i]
		for i in data.rows])

	# hca and plot
	tden_leaves = hca_dendro(tden_axes, plot_data, orientation="top")
	rden_leaves = hca_dendro(rden_axes, plot_data, orientation="right")

	# heatmap
	cmap = matplotlib.colormaps["bwr"]
	p = hmap_axes.pcolor(plot_data.data[numpy.ix_(rden_leaves, tden_leaves)],
		cmap=cmap, vmin=unified_vmin, vmax=unified_vmax
	)

	# add per wwtp significance testing results
	add_1v1_ttest_with_fdr(hmap_axes, data, sig_level=0.01,
		rden_native=plot_data.rows,
		rden_leaves=rden_leaves,
		tden_native=plot_data.cols,
		tden_leaves=tden_leaves,
	)

	# color bar if necessary
	if cbar_axes:
		hmap_axes.figure.colorbar(p, cax=cbar_axes, orientation="horizontal",
			ticks=numpy.linspace(*p.get_clim(), 3), drawedges=False,
			label="CLR",
		)
		cbar_axes.xaxis.tick_top()
		cbar_axes.xaxis.set_label_position("top")

	# heatmap xticklabels
	for x, i in enumerate(tden_leaves):
		s = " " + plot_data.cols[i]
		if oligo_sig_flag[i]:
			s += "*"
			color = "#ff0000"
		else:
			color = "#000000"
		hmap_axes.text(x + 0.5, 0, s, fontsize=10, color=color, rotation=270,
			horizontalalignment="center", verticalalignment="top",
		)

	# heatmap yticklabel and category patch
	for y, i in enumerate(rden_leaves):
		wwtp = plot_data.rows[i]
		color = pylib.supp.CATE_COLOR[pylib.supp.WWTP_CATE[wwtp]]
		# ticklabel
		hmap_axes.text(0, y + 0.5, pylib.supp.WWTP_DISPLAY[wwtp] + " ",
			fontsize=8, color=_adjust_color_brightness(color, 0.75),
			horizontalalignment="right", verticalalignment="center",
		)
		# patch
		p = matplotlib.patches.Rectangle((0, y), 1, 1, edgecolor="none",
			facecolor=color,
		)
		cate_axes.add_patch(p)
	# misc
	cate_axes.set_xlim(0, 1)
	cate_axes.set_ylim(0, plot_data.n_rows)

	if with_cate_legend:
		handles = list()
		for k in sorted(pylib.supp.CATE_COLOR.keys()):
			p = matplotlib.patches.Rectangle((0, 0), 1, 1,
				edgecolor="none", facecolor=pylib.supp.CATE_COLOR[k],
				label=pylib.supp.CATE_DISPLAY[k],
			)
			handles.append(p)
		hmap_axes.legend(handles=handles, loc=3, bbox_to_anchor=[0.98, 1.02],
			handlelength=0.75, fontsize=10, frameon=False,
		)

	return


def plot(png, dpi=300):
	# create layout
	layout = create_layout()
	figure = layout["figure"]

	for orgn_cfg in ORGN_PLOT_CFG:
		orgn_key = orgn_cfg["key"]
		data = pylib.table.RelaAbundTable.from_file(
			os.path.join("data", orgn_cfg["path_key"] + ".clr.tsv")
		)
		otu_list = load_list(
			os.path.join("data", orgn_cfg["path_key"] + ".abund_oligo.list")
		)

		plot_tsne(layout[orgn_key + "_tsne"], data,
			tsne_random_state=orgn_cfg.get("tsne_random_state", None),
			# only the bottom-most one shows x-axis label
			with_xlabel=orgn_key == "tetr", with_ylabel=True,
			anosim_permutations=10000, title=orgn_cfg["display"],
		)

		plot_hmap(data.filter_cols(by_tags=otu_list),
			hmap_axes=layout[orgn_key + "_hmap"],
			tden_axes=layout[orgn_key + "_tden"],
			rden_axes=layout[orgn_key + "_rden"],
			cate_axes=layout[orgn_key + "_cate"],
			with_cate_legend=(orgn_key == "tetr"),
			unified_vmin=-12,
			unified_vmax=12,
			cbar_axes=(layout["cbar"] if orgn_key == "acc" else None),
		)

	# savefig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def main():
	plot("manuscript.fig2.png", dpi=600)
	return


if __name__ == "__main__":
	main()
