#!/usr/bin/env python3

import collections
import json
import numpy
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import matplotlib.patches
import matplotlib.patheffects
import mpllayout
import sklearn.cluster
# custom lib
import pylib


ORGN_PLOT_CFG = [
	{
		"key": "acc",
		"display": "Ca_Accumulibacter (CA)",
		"path_key": "candidatus_accumulibacter",
		"color": "#4040ff",
		"common_seqs": "data/candidatus_accumulibacter.common_seqs.json",
		"oligo_rename": "data/candidatus_accumulibacter.rename.json",
		"jplace": "pplacer/candidatus_accumulibacter.aligned.selected.seqs.jplace",
		"reroot": 595,
		"xspan": 1.1,
		"theta_offset": -2.1,
		"x_offset": 0.35,
		"y_offset": -0.06,
		"scale_len": 0.1,
	},
	{
		"key": "dech",
		"display": "Dechloromonas (DE)",
		"path_key": "dechloromonas",
		"color": "#f0691a",
		"common_seqs": "data/dechloromonas.common_seqs.json",
		"oligo_rename": "data/dechloromonas.rename.json",
		"jplace": "pplacer/dechloromonas.aligned.selected.seqs.jplace",
		"reroot": 508,
		"xspan": 0.95,
		"theta_offset": 0,
		"x_offset": 0.26,
		"y_offset": 0.01,
		"scale_len": 0.1,
	},
	{
		"key": "tetr",
		"display": "Tetrasphaera (TE)",
		"path_key": "tetrasphaera",
		"color": "#0a6e16",
		"common_seqs": "data/tetrasphaera.common_seqs.json",
		"oligo_rename": "data/tetrasphaera.rename.json",
		"jplace": "pplacer/tetrasphaera.aligned.selected.seqs.jplace",
		"reroot": 336,
		"xspan": 0.5,
		"theta_offset": -1.7,
		"x_offset": 0.16,
		"y_offset": 0.008,
		"scale_len": 0.05,
	}
]

TAX_FILE = "data/midas.4.8.1.tax"


def load_json(fname):
	with pylib.util.get_fp(fname, "r") as fp:
		ret = json.load(fp)
	return ret


def load_tax_map(fname: str = TAX_FILE, *, delimiter="\t") -> dict:
	ret = dict()
	with pylib.util.get_fp(fname, "r") as fp:
		for line in fp.read().splitlines():
			seq_id, tax = line.split(delimiter)
			ret[seq_id] = pylib.util.Taxonomy.parse(tax)
	return ret


def _adjust_color_brightness(c: str, multiplier: float) -> str:
	hsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
	hsv[2] = min(hsv[2] * multiplier, 1.0)
	return matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb(hsv))


def create_layout() -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.05,
		right_margin=0.05,
		top_margin=0.05,
		bottom_margin=0.05,
	)

	tree_height = 2.2
	tree_width = 5.0

	for i, orgn_cfg in enumerate(ORGN_PLOT_CFG[::-1]):
		key = orgn_cfg["key"]
		# tree axes
		tree = lc.add_frame(key + "_tree")
		if i == 0:
			tree.set_anchor("bottomleft")
		else:
			prev_tree = lc.get_frame(ORGN_PLOT_CFG[-i]["key"] + "_tree")
			tree.set_anchor("bottomleft", prev_tree, "topleft",
				offsets=(0, 0.05))
		tree.set_size(tree_width, tree_height)

	# create layout
	layout = lc.create_figure_layout()

	# apply axes styles
	for name, axes in layout.items():
		if name == "figure":
			continue

		axes.set_facecolor("#ffffff")
		for sp in axes.spines.values():
			sp.set_visible(True)

		axes.tick_params(
			left=False, labelleft=False,
			right=False, labelright=False,
			top=False, labeltop=False,
			bottom=False, labelbottom=False,
		)

	layout["w2h_ratio"] = tree_width / tree_height

	return layout


def place_tree_node_recursive(node: pylib.jplace_tree.JTree.JTreeNode, *,
		x=0, y=0, theta_0=0, theta_1=numpy.pi * 2, theta_offset=0,
	):
	# place the current node
	theta_c = (theta_0 + theta_1) / 2 + theta_offset
	x += numpy.cos(theta_c) * node.dist
	y += numpy.sin(theta_c) * node.dist
	node.x = x
	node.y = y
	node.theta_0 = theta_0 + theta_offset
	node.theta_1 = theta_1 + theta_offset
	node.theta = theta_c

	# place children
	if not node.is_leaf:
		leaf_part = [i.n_leaves for i in node.children]
		theta_part = numpy.cumsum([0] + leaf_part)
		theta_part = theta_part / theta_part[-1] * (theta_1 - theta_0) + theta_0
		for i, child in enumerate(node.children):
			place_tree_node_recursive(child, x=x, y=y, theta_0=theta_part[i],
				theta_1=theta_part[i + 1], theta_offset=theta_offset)

	return


def _plot_tree(axes: matplotlib.axes.Axes, cfg: dict,
		jtree: pylib.jplace_tree.JTree, *, tax_map: dict = None,
		with_leaf_label=True, show_node_id=False,
	):
	for node in jtree.traverse():
		for c in node.children:
			axes.plot([node.x, c.x], [node.y, c.y], linestyle="-", linewidth=1,
				color="#606060"
			)

		if show_node_id:
			axes.text(node.x, node.y, node.edge_id, fontsize=2,
				horizontalalignment="center", verticalalignment="center")

		# labels
		if with_leaf_label:
			xspan = cfg.get("xspan", 0.01)
			label_min_dist = xspan / 100
			gap = xspan / 50
			if (node.name) and (node.dist > label_min_dist):
				gap_x = numpy.cos(node.theta) * gap
				gap_y = numpy.sin(node.theta) * gap

				if (tax_map is not None) and (node.name in tax_map):
					s = tax_map[node.name].s.split("_", maxsplit=1)[-1]
				else:
					s = node.name

				r = numpy.rad2deg(node.theta) % 360
				if (r > 100) and (r < 260):
					r += 180
					ha = "right"
				else:
					ha = "left"

				axes.text(node.x + gap_x, node.y + gap_y, s,
					fontsize=8, color=getattr(node, "label_color", "#000000"),
					rotation=r, rotation_mode="anchor",
					horizontalalignment=ha, verticalalignment="center",
				)

	return


def _load_seq_to_oligo_map(cfg: dict) -> dict:
	oligo_rename = load_json(cfg["oligo_rename"])
	common_seqs = load_json(cfg["common_seqs"])
	oligo_rename_inv = {v: k for k, v in oligo_rename.items()}
	ret = dict()
	for k, v in common_seqs.items():
		if k not in oligo_rename_inv:
			continue
		for s in v:
			ret[s] = oligo_rename_inv[k]
	return ret


def _plot_placements_by_oligo(axes: matplotlib.axes.Axes, cfg: dict,
		jtree: pylib.jplace_tree.JTree, *, edge_place: dict, label: str = None,
		edgecolor="#000000", facecolor=None, tax_map: dict = None,
	):

	ret = dict()

	# find the root (lca of all placements)
	nodes = [jtree.get_node_by_edge_id(i) for i in edge_place.keys()]
	lca = jtree.lowest_common_ancestor(*nodes)
	r = max([lca.dist_to_child(jtree.get_node_by_edge_id(k)) + v
		for k, v in edge_place.items()
	])

	# find taxonomy statistics
	tax = collections.Counter()
	for i in nodes:
		if (tax_map) and (i.name in tax_map):
			tax[tax_map[i.name].s] += 1
		else:
			tax[i.name] += 1

	ret["tax"] = tax

	if facecolor is None:
		facecolor = edgecolor + "40"

	thetas = [i.theta for i in nodes]
	theta1, theta2 = min(thetas), max(thetas)
	p = matplotlib.patches.Wedge((lca.x, lca.y), r,
		theta1=numpy.rad2deg(theta1),
		theta2=numpy.rad2deg(theta2),
		edgecolor=edgecolor, facecolor=facecolor, linewidth=0.5,
		clip_on=False,
		zorder=3,
		label=label,
	)
	axes.add_patch(p)
	ret["wedge"] = p
	ret["color"] = edgecolor

	p = matplotlib.patches.Circle((lca.x, lca.y), cfg.get("xspan", 0.1) / 200,
		edgecolor=edgecolor, facecolor=facecolor,
		clip_on=False,
		zorder=3,
	)
	axes.add_patch(p)

	min_tr = cfg.get("xspan", 0.01) / 50
	theta = (theta1 + theta2) / 2
	if r <= min_tr:
		theta += (numpy.random.random() - 0.5) * numpy.pi
		tr = min_tr
	else:
		tr = r * 0.80

	axes.text(
		x=lca.x + numpy.cos(theta) * tr,
		y=lca.y + numpy.sin(theta) * tr,
		s=label.split("_")[-1],
		color=_adjust_color_brightness(edgecolor, 1.00),
		zorder=4,
		fontsize=8,
		verticalalignment="center",
		horizontalalignment="center",
		path_effects=[matplotlib.patheffects.withStroke(linewidth=1.5,
			foreground="#ffffff")],
	)

	return ret


def _get_text_width(axes: matplotlib.axes.Axes, text: str, text_props=None,
	) -> float:
	if text_props is None:
		text_props = dict()
	figure = axes.figure
	dummy_text = axes.text(0, 0, text, **text_props)
	bbox = dummy_text.get_window_extent(renderer=figure.canvas.get_renderer())
	w, h = bbox.width, bbox.height
	dummy_text.remove()
	return axes.transAxes.inverted().transform([w, h])


def _add_labelled_pies(axes: matplotlib.axes.Axes, cfg: dict, *, label: str,
		oligo_tax: list, place_plot_res: dict, x: float, y: float,
		label_w: float = 0.0, pie_w: float = 0.05, gap: float = 0.01,
		w2h_ratio: float = 1, label_text_props=None,
	):
	if label_text_props is None:
		label_text_props = dict()

	axes.text(x, y, label, **label_text_props, transform=axes.transAxes,
		horizontalalignment="left", verticalalignment="center")

	xspan = cfg.get("xspan", 0.01)
	yspan = xspan / w2h_ratio
	radius = pie_w / 2 * xspan

	for i, oligo in enumerate(oligo_tax):
		res = place_plot_res[oligo]
		total = res["tax"].total()
		top = res["tax"].most_common()[0][1]

		wx = (x + label_w + gap + pie_w / 2 + (pie_w + gap) * i - 0.5) * xspan
		wy = (y - 0.5) * yspan

		theta2 = 90
		for t, c in res["tax"].most_common():
			theta = c / total * 360
			if theta2 == 90:  # the first wedge is highlighted
				edgecolor = res["color"]
				facecolor = res["color"] + "80"
				zorder = 3
			else:
				edgecolor = res["color"] + "40"
				facecolor = res["color"] + "20"
				zorder = 2
			w = matplotlib.patches.Wedge((wx, wy), r=radius,
				theta1=theta2 - theta, theta2=theta2,
				hatch=("////" if t is None else None),
				edgecolor=edgecolor,
				facecolor=facecolor,
				zorder=zorder,
			)
			theta2 -= theta
			axes.add_patch(w)

		# add label
		axes.text(wx, wy, oligo.split("_")[-1], fontsize=8, color=res["color"],
			horizontalalignment="center", verticalalignment="center",
			zorder=4,
			path_effects=[matplotlib.patheffects.withStroke(linewidth=2,
				foreground="#ffffff")],
		)

	return


def _add_tax_pie(axes: matplotlib.axes.Axes, cfg: dict, place_plot_res: dict,
		w2h_ratio: float = 1,
	):
	# assort placement results by taxonomy
	tax_oligo = collections.defaultdict(list)
	for k, v in place_plot_res.items():
		if not v["tax"]:
			continue
		tax_oligo[v["tax"].most_common()[0][0]].append(k)

	# add rows of label+pies
	pie_w = 0.05
	gap = 0.01
	label_gap = 0.005
	post_gap = 0.02
	max_space = 0.60
	w_min = 0.02
	h_max = 0.82
	h_space = 0.12

	row_space_used = list()
	label_text_props = dict(fontsize=8, color="#000000")
	for tax in sorted(tax_oligo.keys(),
			key=lambda x: "\xff" if x is None else x):
		# use \xff so that None will appear all other tax names

		t = "[unknown]: " if tax is None else tax + ":"

		if not t.startswith("midas"):
			t = t.split("_")[-1]
		# estimate required spqce
		l_w, l_h = _get_text_width(axes, t, text_props=label_text_props)
		l_w += label_gap
		n_pies = len(tax_oligo[tax])
		total_w = l_w + pie_w * n_pies + gap * n_pies + post_gap
		for i, rsu in enumerate(row_space_used):
			if (tax is not None) and (rsu + total_w < max_space):
				_add_labelled_pies(axes, cfg, label=t, oligo_tax=tax_oligo[tax],
					place_plot_res=place_plot_res, label_w=l_w, pie_w=pie_w,
					x=rsu, y=h_max - (h_space * i), gap=gap,
					w2h_ratio=w2h_ratio, label_text_props=label_text_props,
				)
				row_space_used[i] += total_w
				break
		else:
			_add_labelled_pies(axes, cfg, label=t, oligo_tax=tax_oligo[tax],
				place_plot_res=place_plot_res, label_w=l_w,
				x=w_min, y=h_max - (h_space * len(row_space_used)),
				pie_w=pie_w, gap=gap, w2h_ratio=w2h_ratio,
				label_text_props=label_text_props,
			)
			row_space_used.append(total_w + w_min)
	return


def _plot_placements(axes: matplotlib.axes.Axes, cfg: dict,
		jtree: pylib.jplace_tree.JTree, placements: dict, tax_map: dict = None,
		w2h_ratio: float = 1,
	):
	seq_to_oligo = _load_seq_to_oligo_map(cfg)

	# count placements on each node
	place = collections.defaultdict(dict)
	for i in placements:
		for nm in i["nm"]:
			oligo = seq_to_oligo[nm[0]]
			for p in i["p"]:
				dist, edge = p[:2]
				place[oligo][edge] = max(dist, place[oligo].get(edge, 0))

	color_list = matplotlib.colormaps["tab10"].colors \
		+ matplotlib.colormaps["Dark2"].colors
	color_list = [matplotlib.colors.to_hex(i) for i in color_list]

	oligo = sorted(place.keys(), key=lambda x: int(x.split("_")[1]))
	place_plot_res = dict()
	for k, c in zip(oligo, color_list):
		res = _plot_placements_by_oligo(axes, cfg, jtree, edge_place=place[k],
			label=k, edgecolor=c, tax_map=tax_map,
		)
		place_plot_res[k] = res

	_add_tax_pie(axes, cfg, place_plot_res, w2h_ratio=w2h_ratio)

	# legend
	# handles = [place_plot_res[k]["wedge"] for k in sorted(place_plot_res.keys())
	# axes.legend(handles=handles, loc=2, bbox_to_anchor=[0.02, 0.90], ncols=2,
	# fontsize=8, frameon=False, handlelength=0.75,
	# )

	return


def plot_tree_complex(axes: matplotlib.axes.Axes, cfg: dict,
		tax_map: dict = None, w2h_ratio: float = 1, theta_offset: float = 0,
	):
	# load tree data from jplace file
	jplace = load_json(cfg["jplace"])
	jtree = pylib.jplace_tree.JTree.parse(jplace["tree"])
	jtree.reroot(jtree.get_node_by_edge_id(cfg["reroot"]))

	place_tree_node_recursive(jtree.root,
		x=cfg.get("x_offset", 0), y=cfg.get("y_offset", 0),
		theta_offset=cfg.get("theta_offset", 0),
	)

	# plot major components of the tree and placements
	_plot_tree(axes, cfg, jtree, tax_map=tax_map, with_leaf_label=False,
		show_node_id=False)
	_plot_placements(axes, cfg, jtree, jplace["placements"], tax_map=tax_map,
		w2h_ratio=w2h_ratio)

	# misc
	if "xspan" in cfg:
		xrad = cfg["xspan"] / 2
		yrad = xrad / w2h_ratio
		axes.set_xlim(-xrad, xrad)
		axes.set_ylim(-yrad, yrad)
	axes.text(0.01, 0.98, cfg["display"], fontsize=12, color="#000000",
		transform=axes.transAxes,
		horizontalalignment="left", verticalalignment="top",
	)

	# add scale
	scale_len = cfg.get("scale_len", 0.1)
	x, y = axes.transData.inverted().transform(
		axes.transAxes.transform((0.02, 0.02 * w2h_ratio))
	)
	axes.plot([x, x + scale_len], [y, y], linestyle="-", linewidth=1.0,
		marker="|", color="#000000",
	)

	_, y = axes.transData.inverted().transform(
		axes.transAxes.transform((0, 0.025 * w2h_ratio))
	)
	axes.text(x + scale_len / 2, y, str(scale_len),
		fontsize=10, color="#000000",
		horizontalalignment="center", verticalalignment="bottom",
	)

	return


def plot(png, *, tax_map: dict = None, dpi=300):

	# create layout
	layout = create_layout()
	figure = layout["figure"]

	for cfg in ORGN_PLOT_CFG:
		plot_tree_complex(layout[cfg["key"] + "_tree"], cfg, tax_map=tax_map,
			w2h_ratio=layout["w2h_ratio"])

	# savefig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def main():
	tax_map = load_tax_map()
	numpy.random.seed(123864)
	plot("manuscript.fig3.png", tax_map=tax_map, dpi=600)
	return


if __name__ == "__main__":
	main()
