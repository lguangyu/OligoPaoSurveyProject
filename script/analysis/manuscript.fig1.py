#!/usr/bin/env python3

import collections
import numpy
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import matplotlib.patches
import mpllayout
import scipy.cluster
import scipy.stats
import sys
# custom lib
import pylib


OTU_PLOT_CFG = [
	{
		"key": "Ca_Accumulibacter",
		"display": "candidatus_accumulibacter",
		"color": "#4040ff"
	},
	{
		"key": "Dechloromonas",
		"display": "dechloromonas",
		"color": "#f0691a"
	},
	{
		"key": "Tetrasphaera",
		"display": "tetrasphaera",
		"color": "#0a6e16"
	}
]

SIG_LEVEL_CFG = {
	"breaks": [
		{
			"mark": "**",
			"alpha": 0.01,
			"color": "#ff0000",
		},
		{
			"mark": "*",
			"alpha": 0.05,
			"color": "#4040ff",
		},
	],
	"default_color": "#606060",
}


def load_abund_table(fname):
	return pylib.table.RelaAbundTable.from_file(fname)


def create_layout() -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.7,
		right_margin=0.3,
		top_margin=0.3,
		bottom_margin=2.0,
	)

	bv_col_width = 0.20
	bv_height = 3.0
	heatmap_width = 2.5

	# bar axes
	bar_ebpr = lc.add_frame("bar_ebpr")
	bar_ebpr.set_anchor("bottomleft")
	bar_ebpr.set_size(bv_col_width * 8, bv_height)

	bar_pilotebpr = lc.add_frame("bar_pilotebpr")
	bar_pilotebpr.set_anchor("bottomleft", bar_ebpr, "bottomright",
		offsets=(0.1, 0))
	bar_pilotebpr.set_size(bv_col_width, bv_height)

	bar_pilots2ebpr = lc.add_frame("bar_pilots2ebpr")
	bar_pilots2ebpr.set_anchor("bottomleft", bar_pilotebpr, "bottomright",
		offsets=(0, 0))
	bar_pilots2ebpr.set_size(bv_col_width, bv_height)

	bar_s2ebpr = lc.add_frame("bar_s2ebpr")
	bar_s2ebpr.set_anchor("bottomleft", bar_pilots2ebpr, "bottomright",
		offsets=(0.1, 0))
	bar_s2ebpr.set_size(bv_col_width * 5, bv_height)

	# violin axes
	violin_ebpr = lc.add_frame("violin_ebpr")
	violin_ebpr.set_anchor("bottomleft", bar_ebpr, "topleft",
		offsets=(0, 0.2))
	violin_ebpr.set_size(bv_col_width * 8, bv_height)

	violin_pilotebpr = lc.add_frame("violin_pilotebpr")
	violin_pilotebpr.set_anchor("bottomleft", violin_ebpr, "bottomright",
		offsets=(0.1, 0))
	violin_pilotebpr.set_size(bv_col_width, bv_height)

	violin_pilots2ebpr = lc.add_frame("violin_pilots2ebpr")
	violin_pilots2ebpr.set_anchor("bottomleft", violin_pilotebpr, "bottomright",
		offsets=(0, 0))
	violin_pilots2ebpr.set_size(bv_col_width, bv_height)

	violin_s2ebpr = lc.add_frame("violin_s2ebpr")
	violin_s2ebpr.set_anchor("bottomleft", violin_pilots2ebpr, "bottomright",
		offsets=(0.1, 0))
	violin_s2ebpr.set_size(bv_col_width * 5, bv_height)

	# cbar axes
	cbar = lc.add_frame("cbar")
	cbar.set_anchor("topleft", violin_s2ebpr, "topright",
		offsets=(0.3, -0.4))
	cbar.set_size(2.5, 0.15)

	# dendro axes
	dendro = lc.add_frame("dendro")
	dendro.set_anchor("topleft", cbar, "bottomleft",
		offsets=(0, -0.3))
	dendro.set_anchor("bottomright", cbar, "bottomright",
		offsets=(0, -1.5))

	# heatmap axes
	heatmap = lc.add_frame("heatmap")
	heatmap.set_anchor("bottomleft", bar_s2ebpr, "bottomright",
		offsets=(0.3, 0))
	heatmap.set_anchor("topright", dendro, "bottomright",
		offsets=(0, 0))

	# fdr axes
	fdr = lc.add_frame("fdr")
	fdr.set_anchor("topleft", cbar, "topright",
		offsets=(0.8, 0))
	fdr.set_size(1.6, 1.2)

	# create layout
	layout = lc.create_figure_layout()
	for name, axes in layout.items():
		if name == "figure":
			continue

		axes.set_facecolor("#ffffff")

		if name == "fdr":
			for sp in axes.spines.values():
				sp.set_visible(True)
			axes.tick_params(
				left=True, labelleft=True,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)
		elif name == "dendro":
			for sp in axes.spines.values():
				sp.set_visible(False)
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
		elif name == "heatmap":
			for sp in axes.spines.values():
				sp.set_visible(True)
			axes.tick_params(
				left=False, labelleft=False,
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)
		else:
			# induce the position of axes
			border_vis = dict(
				left=("_ebpr" in name),
				right=("_s2ebpr" in name),
				top=True,
				bottom=True,
			)
			# axes spines
			for k, v in border_vis.items():
				axes.spines[k].set_visible(v)

			axes.tick_params(
				left=border_vis["left"], labelleft=border_vis["left"],
				right=False, labelright=False,
				top=False, labeltop=False,
				bottom=False, labelbottom=False,
			)

	# apply axes styles

	return layout


def _adjust_color_brightness(c: str, multiplier: float):
	hsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
	hsv[2] = min(hsv[2] * multiplier, 1.0)
	return matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb(hsv))


def _plot_styled_box_element(axes: matplotlib.axes.Axes, *, vals, position,
		width, linecolor, label, zorder
	) -> dict:
	p = axes.boxplot([vals], positions=[position], widths=[width],
		patch_artist=True, vert=True, labels=[label], zorder=zorder,
		boxprops={
			"edgecolor": linecolor,
			"facecolor": linecolor + "40",
		},
		flierprops={
			"marker": "x",
			"markersize": 5,
			"markeredgecolor": linecolor,
		},
		medianprops={
			"color": linecolor,
		},
		whiskerprops={
			"color": linecolor,
		},
		capprops={
			"color": linecolor,
		},
	)
	return p


def plot_box(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable, *,
		wwtp_cate: str, otu_cfg_list: list, ymin=0, ymax=1,
		with_each_wwtp=True, with_overall=True, with_ylabel=False,
		with_legend=False,
	):

	# extract otus
	data = data.filter_cols(by_tags=[i["key"] for i in otu_cfg_list])

	# sort by wwtp
	data_by_wwtp = dict()
	for k, v in pylib.supp.WWTP_CATE.items():
		if v != wwtp_cate:
			continue
		mask = [pylib.supp.SAMPLE_WWTP[i] == k for i in data.rows]
		data_by_wwtp[k] = data.filter_rows(by_index=mask)

	xi = 0  # left-boundry of current drawing column
	xpad = 0.05  # half-padding between elements of neigboring columns
	xw = (1 - xpad * 2) / len(otu_cfg_list)  # width of elements
	cate_color = pylib.supp.CATE_COLOR[wwtp_cate]
	darkened_cate_color = _adjust_color_brightness(cate_color, 0.70)

	if with_each_wwtp:
		for wwtp in sorted(data_by_wwtp.keys()):
			for i, otu_cfg in enumerate(otu_cfg_list):
				color = otu_cfg["color"]
				# plot box
				d = data_by_wwtp[wwtp].filter_cols(by_tags=[otu_cfg["key"]])\
					.data.ravel()
				p = _plot_styled_box_element(axes, vals=d, linecolor=color,
					position=xi + xpad + xw * (i + 0.5), width=xw * 0.9,
					label=otu_cfg["display"], zorder=3
				)

			# add wwtp background
			p = matplotlib.patches.Rectangle((xi, ymin), 1, ymax - ymin,
				edgecolor="none",
				facecolor=cate_color + ("45" if (xi & 0x1) else "30"), zorder=1)
			axes.add_patch(p)

			xi += 1

	# add overall column if necessary
	if with_overall:
		for i, otu_cfg in enumerate(otu_cfg_list):
			color = otu_cfg["color"]
			# plot box
			d = numpy.concatenate(
				[v.filter_cols(by_tags=[otu_cfg["key"]]).data.ravel()
					for v in data_by_wwtp.values()
				]
			)
			p = _plot_styled_box_element(axes, vals=d, linecolor=color,
				position=xi + xpad + xw * (i + 0.5), width=xw * 0.9,
				label=otu_cfg["display"], zorder=3
			)

		# add wwtp background
		p = matplotlib.patches.Rectangle((xi, ymin), 1, ymax - ymin,
			edgecolor=("none" if not with_each_wwtp else darkened_cate_color),
			facecolor=cate_color + ("45" if (xi & 0x1) else "30"),
			clip_on=False, zorder=1)
		axes.add_patch(p)

		xi += 1

	# legend if necessary
	if with_legend:
		handles = list()
		for otu_cfg in otu_cfg_list:
			color = otu_cfg["color"]
			handles.append(
				matplotlib.patches.Rectangle((0, ymin), 1, ymax,
					edgecolor=color, facecolor=color + "40",
					label=otu_cfg["key"]
			))
		axes.legend(handles=handles, loc=1, bbox_to_anchor=[1.02, 1.02],
			handlelength=0.75, fontsize=12, frameon=True, edgecolor="none",
			facecolor="#ffffff80",
		)

	# misc
	axes.set_xlim(0, xi)
	axes.set_ylim(ymin, ymax)
	if with_ylabel:
		axes.set_ylabel("Relative abundance", fontsize=12)
	# add xticklabels manually
	xticklabels = list()
	if with_each_wwtp:
		xticklabels.extend([pylib.supp.WWTP_DISPLAY[i]
			for i in sorted(data_by_wwtp.keys())])
	if with_overall:
		xticklabels.append(pylib.supp.CATE_DISPLAY[wwtp_cate]
			+ (" overall" if with_each_wwtp else ""))
	for i, label in enumerate(xticklabels):
		axes.text(i + 0.5, 0, label + " ", color=darkened_cate_color,
			fontsize=12, rotation=90,
			horizontalalignment="center", verticalalignment="top"
		)

	return


def _plot_styled_violin_element(axes: matplotlib.axes.Axes, *, vals, position,
		width, linecolor, label, zorder
	) -> dict:
	if not len(vals):
		return

	v = axes.violinplot([vals], positions=[position], widths=[width], vert=True)
	for p in v["bodies"]:
		p.set(facecolor=linecolor + "40", edgecolor=linecolor, alpha=None,
			zorder=zorder)
	for n in ["cmaxes", "cbars", "cmins"]:
		v[n].set(linewidth=1.0, color=linecolor, zorder=zorder + 1)
	return v


def plot_violin(axes: matplotlib.axes.Axes, data: pylib.table.RelaAbundTable, *,
		wwtp_cate: str, ymin: float, ymax: float, with_each_wwtp=True,
		with_overall=True, with_ylabel=False,
	):

	# sort by wwtp
	data_by_wwtp = dict()
	for k, v in pylib.supp.WWTP_CATE.items():
		if v != wwtp_cate:
			continue
		mask = [pylib.supp.SAMPLE_WWTP[i] == k for i in data.rows]
		data_by_wwtp[k] = scipy.stats.entropy(
			data.filter_rows(by_index=mask).data, axis=1,
		)
		assert len(data_by_wwtp[k]) == sum(mask)

	xi = 0  # left-boundry of current drawing column
	cate_color = pylib.supp.CATE_COLOR[wwtp_cate]
	darkened_cate_color = _adjust_color_brightness(cate_color, 0.70)

	if with_each_wwtp:
		for wwtp in sorted(data_by_wwtp.keys()):
			p = _plot_styled_violin_element(axes, vals=data_by_wwtp[wwtp],
				linecolor=cate_color, position=xi + 0.5, width=0.8,
				label=pylib.supp.WWTP_DISPLAY[wwtp], zorder=3
			)

			# add wwtp background
			p = matplotlib.patches.Rectangle((xi, ymin), 1, ymax - ymin,
				edgecolor="none",
				facecolor=cate_color + ("45" if (xi & 0x1) else "30"), zorder=1)
			axes.add_patch(p)

			xi += 1

	# add overall column if necessary
	if with_overall:
		d = numpy.concatenate(list(data_by_wwtp.values()))
		p = _plot_styled_violin_element(axes, vals=d,
			linecolor=cate_color, position=xi + 0.5, width=0.8,
			label=pylib.supp.CATE_DISPLAY[wwtp_cate]
				+ (" overall" if with_each_wwtp else ""),
			zorder=3
		)

		# add wwtp background
		p = matplotlib.patches.Rectangle((xi, ymin), 1, ymax - ymin,
			edgecolor=("none" if not with_each_wwtp else darkened_cate_color),
			facecolor=cate_color + ("45" if (xi & 0x1) else "30"),
			clip_on=False, zorder=1)
		axes.add_patch(p)

		xi += 1

	# misc
	axes.set_xlim(0, xi)
	axes.set_ylim(ymin, ymax)
	if with_ylabel:
		axes.set_ylabel("Genus alpha diversity (Shannon)", fontsize=12)

	return


def plot_dendro(axes: matplotlib.axes.Axes, data: numpy.ndarray, *,
		clr_prior=None, right_padding: int = 0) -> numpy.ndarray:
	# in this context, row=wwtp, col=otu
	# it might be transposed in caller's context
	n_row, n_col = data.shape

	if axes is None:
		return numpy.arange(n_row)  # the null order

	# clr transformation
	if clr_prior is None:
		clr_prior = 0.01 / n_col  # the total row-sum is 0.01 bigger
	log_data = numpy.log(data + clr_prior)
	clr = log_data - log_data.mean(axis=1, keepdims=True)
	# hca
	linkage = scipy.cluster.hierarchy.linkage(clr, method="average",
		optimal_ordering=True)
	dendro = scipy.cluster.hierarchy.dendrogram(linkage, orientation="top",
		no_plot=True)
	leaves = numpy.asarray(dendro["leaves"])

	# plot dendrogram
	for x, y in zip(dendro["icoord"], dendro["dcoord"]):
		axes.plot(x, y, clip_on=False, linestyle="-", linewidth=1.0,
			color="#200064")

	# misc
	axes.set_xlim(0, 10 * (len(leaves) + right_padding))
	axes.set_ylim(0, numpy.max(dendro["dcoord"]) * 1.1)

	return leaves


def abund_ttest_with_fdr(data, cate_grouping: list, max_n_otu: int,
		sig_level_cfg: dict,
	) -> (numpy.ndarray, numpy.ndarray):

	mask = list()
	for cate in cate_grouping:
		m = [pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[i]] == cate
			for i in data.rows
		]
		mask.append(m)

	# ttest for each row
	ttest_res = scipy.stats.ttest_ind(
		data.filter_rows(by_index=mask[0]).data[:, :max_n_otu],
		data.filter_rows(by_index=mask[1]).data[:, :max_n_otu],
		axis=0, equal_var=False,
	)

	# fdr thresholding
	pval = ttest_res.pvalue
	argsort_idx = numpy.argsort(pval)
	sorted_pval = pval[argsort_idx]
	n_pval = len(pval)
	sig_level = numpy.full(n_pval, len(sig_level_cfg["breaks"]), dtype=int)
	for cfg in sig_level_cfg["breaks"]:
		fdr_thres = numpy.linspace(0, cfg["alpha"], n_pval + 1)[1:]
		sig_level[sorted_pval <= fdr_thres] -= 1

	sig_flag = numpy.empty(n_pval, dtype=int)
	sig_flag[argsort_idx] = sig_level

	return pval, sig_flag


def plot_abund_ttest_fdr(axes: matplotlib.axes.Axes, pval: numpy.ndarray,
		sig_flag: numpy.ndarray, sig_level_cfg: dict,
	):
	if axes is None:
		return

	n_pval = len(pval)
	argsort_idx = numpy.argsort(pval)
	sorted_pval = pval[argsort_idx]
	n_sig_level = len(sig_level_cfg["breaks"])

	edgecolors = [
		sig_level_cfg["default_color"]
			if i == n_sig_level
			else sig_level_cfg["breaks"][i]["color"]
		for i in sig_flag[argsort_idx]
	]
	facecolors = [i + "40" for i in edgecolors]

	x = numpy.arange(1, n_pval + 1)
	axes.scatter(x, sorted_pval, clip_on=False, s=40, marker="o",
		edgecolors=edgecolors, facecolors=facecolors, zorder=2)

	# plot threshold lines and legend
	x = numpy.arange(n_pval + 1)
	handles = list()
	labelcolor = list()
	for i, cfg in enumerate(sig_level_cfg["breaks"]):
		a = cfg["alpha"]
		c = cfg["color"]
		y = numpy.linspace(0, a, n_pval + 1)
		p = axes.plot(x, y, linestyle="-", linewidth=1, color=c, zorder=3,
			label=(r"$\alpha_{%u}$=%.4f" % (i + 1, a)).rstrip("0"),
		)[0]
		handles.append(p)
		labelcolor.append(c)
	axes.legend(handles=handles, loc=2, bbox_to_anchor=[0.0, 1.0],
		fontsize=10, frameon=False, labelcolor=labelcolor,
	)

	# misc
	axes.grid(linewidth=1.0, color="#c0c0c0")
	axes.set_xlim(0, n_pval)
	axes.set_ylim(0, 1.0)

	axes.set_xlabel(r"$i$-th smallest p-value", fontsize=12)
	axes.set_ylabel("p-value", fontsize=12)

	return


def plot_abund(axes: matplotlib.axes.Axes, data, max_n_otu: int, *,
		cate_list: list, with_overall: list = None,
		dendrogram_axes: matplotlib.axes.Axes = None,
		cbar_axes: matplotlib.axes.Axes = None, cmap: str = "Spectral_r",
		vmin: float = None, vmax: float = None, zorder: int = 2,
		fdr_axes: matplotlib.axes.Axes = None, sig_level_cfg: dict = None,
	):
	wwtp_cate = pylib.supp.WWTP_CATE
	sig_level_cfg["breaks"].sort(key=lambda x: x["alpha"])

	# group data by category
	xlabels = list()
	xcolors = list()
	plot_data = list()
	for cate in cate_list:
		for wwtp in sorted(wwtp_cate.keys()):
			if wwtp_cate[wwtp] != cate:
				continue
			xlabels.append(pylib.supp.WWTP_DISPLAY[wwtp])
			xcolors.append(pylib.supp.CATE_COLOR[cate])
			# extract data and find the mean
			mask = [pylib.supp.SAMPLE_WWTP[i] == wwtp for i in data.rows]
			plot_data.append(data.filter_rows(by_index=mask)
				.data.mean(axis=0, keepdims=True))

	# add dendrogram
	sample_reorder = plot_dendro(dendrogram_axes, numpy.vstack(plot_data),
		right_padding=len(with_overall),
	)
	xlabels = [xlabels[i] for i in sample_reorder]
	xcolors = [xcolors[i] for i in sample_reorder]
	plot_data = [plot_data[i] for i in sample_reorder]

	# add overall if necessary
	for cate in with_overall:
		xlabels.append(pylib.supp.CATE_DISPLAY[cate] + " overall")
		xcolors.append(_adjust_color_brightness(pylib.supp.CATE_COLOR[cate],
			0.70))
		# extract data and find the mean
		mask = [pylib.supp.WWTP_CATE[pylib.supp.SAMPLE_WWTP[i]] == cate
			for i in data.rows]
		plot_data.append(data.filter_rows(by_index=mask)
			.data.mean(axis=0, keepdims=True))

	plot_data = numpy.vstack(plot_data)[:, :max_n_otu].T

	# ttest
	pval, sig_flag = abund_ttest_with_fdr(data, max_n_otu=max_n_otu,
		cate_grouping=("EBPR", "S2EBPR"), sig_level_cfg=sig_level_cfg,
	)
	plot_abund_ttest_fdr(fdr_axes, pval=pval, sig_flag=sig_flag,
		sig_level_cfg=sig_level_cfg,
	)

	# plot heatmap
	if not isinstance(cmap, matplotlib.colors.Colormap):
		cmap = matplotlib.colormaps[cmap]
	p = axes.pcolor(plot_data, cmap=cmap, vmin=vmin, vmax=vmax)

	# colorbar
	if cbar_axes:
		cbar = axes.figure.colorbar(p, cax=cbar_axes, orientation="horizontal",
			ticks=numpy.linspace(vmin, vmax, 5), drawedges=False,
			label="Genus relative abundance",
		)
		cbar.set_label("Genus relative abundance", fontsize=12)
		cbar_axes.xaxis.tick_top()
		cbar_axes.xaxis.set_label_position("top")

	# misc
	n_row, n_col = plot_data.shape
	axes.set_xlim(0, n_col)
	axes.set_ylim(n_row, 0)

	# add separation line for overalls
	axes.axvline(n_col - len(with_overall), linewidth=2.0, color="#ffffff",
		zorder=zorder + 1)

	# xticklabels
	for i, label in enumerate(xlabels):
		color = _adjust_color_brightness(xcolors[i], 0.70)
		axes.text(i + 0.5, n_row, label + " ", color=color, fontsize=12,
			rotation=90, horizontalalignment="center", verticalalignment="top",
		)

	# yticklabels
	n_sig_level = len(sig_level_cfg["breaks"])
	for i, label in enumerate(data.cols[:max_n_otu]):
		if sig_flag[i] == n_sig_level:
			color = sig_level_cfg["default_color"]
		else:
			cfg = sig_level_cfg["breaks"][sig_flag[i]]
			label += cfg["mark"]
			color = cfg["color"]
		axes.text(n_col, i + 0.5, " " + label, fontsize=10, color=color,
			horizontalalignment="left", verticalalignment="center",
		)

	return


def plot(png, data: pylib.table.RelaAbundTable, dpi=300):
	# create layout
	layout = create_layout()
	figure = layout["figure"]

	# plot boxes
	plot_box(layout["bar_ebpr"], data, wwtp_cate="EBPR",
		ymax=0.3, otu_cfg_list=OTU_PLOT_CFG,
		with_ylabel=True,
	)
	plot_box(layout["bar_pilotebpr"], data, wwtp_cate="RC_EBPR",
		ymax=0.3, otu_cfg_list=OTU_PLOT_CFG,
		with_each_wwtp=True, with_overall=False,
	)
	plot_box(layout["bar_pilots2ebpr"], data, wwtp_cate="RC_PILOT",
		ymax=0.3, otu_cfg_list=OTU_PLOT_CFG,
		with_each_wwtp=True, with_overall=False,
	)
	plot_box(layout["bar_s2ebpr"], data, wwtp_cate="S2EBPR",
		ymax=0.3,
		otu_cfg_list=OTU_PLOT_CFG, with_legend=True,
	)

	# plot violins
	plot_violin(layout["violin_ebpr"], data, wwtp_cate="EBPR",
		ymin=3.5, ymax=6.0, with_ylabel=True,
	)
	plot_violin(layout["violin_pilotebpr"], data, wwtp_cate="RC_EBPR",
		ymin=3.5, ymax=6.0, with_each_wwtp=True, with_overall=False,
	)
	plot_violin(layout["violin_pilots2ebpr"], data, wwtp_cate="RC_PILOT",
		ymin=3.5, ymax=6.0, with_each_wwtp=True, with_overall=False,
	)
	plot_violin(layout["violin_s2ebpr"], data, wwtp_cate="S2EBPR",
		ymin=3.5, ymax=6.0,
	)

	# plot abund heatmap
	plot_abund(layout["heatmap"], data, max_n_otu=16,
		cate_list=["EBPR", "RC_EBPR", "RC_PILOT", "S2EBPR"],
		with_overall=["EBPR", "S2EBPR"],
		dendrogram_axes=layout["dendro"],
		cbar_axes=layout["cbar"], vmin=0.0, vmax=0.2,
		fdr_axes=layout["fdr"], sig_level_cfg=SIG_LEVEL_CFG,
	)

	# savefig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def main():
	abund = load_abund_table("data/genus.abund.tsv")
	plot("manuscript.fig1.png", abund, dpi=600)
	return


if __name__ == "__main__":
	main()
