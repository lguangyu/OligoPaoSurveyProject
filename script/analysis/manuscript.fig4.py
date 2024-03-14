#!/usr/bin/env python3

import abc
import collections
import dataclasses
import functools
import itertools
import numpy
import matplotlib
import matplotlib.artist
import matplotlib.colors
import matplotlib.legend
import matplotlib.pyplot
import matplotlib.patches
import matplotlib.patheffects
import mpllayout
import scipy.stats
# custom lib
import pylib


CHORD_PLOT_CFG = [
	{
		"key": "oligo-oligo",
		"min_r": 0.5,
		"alpha": [0.01, 0.05],
		"groups": [
			{
				"key": "CA_Accumulibacter",
				"display": "CA_Accumulibacter",
				"color": "#4040ff",
				"file": "data/candidatus_accumulibacter.clr.tsv",
				"col_select": "data/candidatus_accumulibacter.abund_oligo.list",
				"col_translate": "performance/candidatus_accumulibacter.col_translate.tsv",
			},
			{
				"key": "Dechloromonas",
				"display": "Dechloromonas",
				"color": "#f0691a",
				"file": "data/dechloromonas.clr.tsv",
				"col_select": "data/dechloromonas.abund_oligo.list",
				"col_translate": "performance/dechloromonas.col_translate.tsv",
			},
			{
				"key": "Tetrasphaera",
				"display": "Tetrasphaera",
				"color": "#0a6e16",
				"file": "data/tetrasphaera.clr.tsv",
				"col_select": "data/tetrasphaera.abund_oligo.list",
				"col_translate": "performance/tetrasphaera.col_translate.tsv",
			},
		],
		"corr_meth": scipy.stats.pearsonr,
		"corr_res_prefix": "misc_out/manuscript.fig4.oligo-oligo.corr",
		"anchor": "CA_Accumulibacter",
		"anchor_angle": 90,
	},
	{
		"key": "oligo-perf",
		"theta_offset": -18,
		"min_r": 0.5,
		"alpha": [0.01, 0.05],
		"groups": [
			{
				"key": "CA_Accumulibacter",
				"display": "CA_Accumulibacter",
				"color": "#4040ff",
				"file": "performance/candidatus_accumulibacter.rep_combined.tsv",
				"col_select": "data/candidatus_accumulibacter.abund_oligo.list",
				"col_translate": "performance/candidatus_accumulibacter.col_translate.tsv",
			},
			{
				"key": "Dechloromonas",
				"display": "Dechloromonas",
				"color": "#f0691a",
				"file": "performance/dechloromonas.rep_combined.tsv",
				"col_select": "data/dechloromonas.abund_oligo.list",
				"col_translate": "performance/dechloromonas.col_translate.tsv",
			},
			{
				"key": "Tetrasphaera",
				"display": "Tetrasphaera",
				"color": "#0a6e16",
				"file": "performance/tetrasphaera.rep_combined.tsv",
				"col_select": "data/tetrasphaera.abund_oligo.list",
				"col_translate": "performance/tetrasphaera.col_translate.tsv",
			},
			{
				"key": "Performance",
				"display": "Performance",
				"color": "#bf0aa1",
				"file": "performance/perf.sample_matched.tsv",
				"col_translate": "performance/perf.col_translate.tsv",
			},
		],
		"group_corr_ipairs": [[0, 3], [1, 3], [2, 3]],
		"corr_meth": scipy.stats.pearsonr,
		"corr_res_prefix": "misc_out/manuscript.fig4.oligo-perf.corr",
		"anchor": "Performance",
		"anchor_angle": 270,
	},
	{
		"enable": False,
		"key": "oligo-batch",
		"theta_offset": -20,
		"min_r": 0.5,
		"alpha": [0.01, 0.05],
		"groups": [
			{
				"key": "CA_Accumulibacter",
				"display": "CA_Accumulibacter",
				"color": "#4040ff",
				"file": "performance/candidatus_accumulibacter.date_combined.tsv",
				"col_select": "data/candidatus_accumulibacter.abund_oligo.list",
				"col_translate": "performance/candidatus_accumulibacter.col_translate.tsv",
			},
			{
				"key": "Dechloromonas",
				"display": "Dechloromonas",
				"color": "#f0691a",
				"file": "performance/dechloromonas.date_combined.tsv",
				"col_select": "data/dechloromonas.abund_oligo.list",
				"col_translate": "performance/dechloromonas.col_translate.tsv",
			},
			{
				"key": "Tetrasphaera",
				"display": "Tetrasphaera",
				"color": "#0a6e16",
				"file": "performance/tetrasphaera.date_combined.tsv",
				"col_select": "data/tetrasphaera.abund_oligo.list",
				"col_translate": "performance/tetrasphaera.col_translate.tsv",
			},
			{
				"key": "Batch",
				"display": "Kinetics/Stoichiometry",
				"color": "#821c05",
				"file": "performance/ac_batch_test.tooker.tsv",
				"col_translate": "performance/ac_batch_test.tooker.col_translate.tsv",
			},
		],
		"group_corr_ipairs": [[0, 3], [1, 3], [2, 3]],
		"corr_meth": scipy.stats.pearsonr,
		"corr_res_prefix": "misc_out/manuscript.fig4.oligo-batch.corr",
		"anchor": "Batch",
		"anchor_angle": 270,
	},
]


def load_map(fname: str, *, delimiter="\t") -> dict:
	ret = dict()
	with pylib.util.get_fp(fname, "r") as fp:
		for line in fp.read().splitlines():
			k, v = line.split(delimiter)
			ret[k] = v
	return ret


def load_col_translate_safe(group_cfg: dict, *ka, **kw) -> dict:
	fname = group_cfg.get("col_translate", None)
	return dict() if fname is None else load_map(fname, *ka, **kw)


def _adjust_color_brightness(c: str, multiplier: float) -> str:
	hsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
	hsv[2] = min(hsv[2] * multiplier, 1.0)
	return matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb(hsv))


@dataclasses.dataclass
class CorrRes(object):
	node1: str
	node2: str
	r: float
	pvalue: float
	sig_flag: float = numpy.nan


def iter_corr_ipairs(cfg: dict):
	if "group_corr_ipairs" in cfg:
		return iter(cfg["group_corr_ipairs"])
	else:
		return itertools.combinations(range(len(cfg["groups"])), 2)
	return


def iter_corr_pairs(cfg: dict):
	return map(lambda x: (cfg["groups"][x[0]], cfg["groups"][x[1]]),
		iter_corr_ipairs(cfg))


def corr_analysis_with_fdr_control(cfg: dict) -> list:
	res_list = list()

	# calculate
	meth = cfg.get("corr_meth", scipy.stats.pearsonr)
	n_groups = len(cfg["groups"])

	for g1, g2 in iter_corr_pairs(cfg):
		t1 = g1["table"]
		t2 = g2["table"]
		# check if the node (col) names are different
		col_collision = set.intersection(set(t1.cols), set(t2.cols))
		if col_collision:
			raise ValueError("column name collision: %s" % col_collision)
		for i1, c1 in enumerate(t1.cols):
			for i2, c2 in enumerate(t2.cols):
				# get the data column from table
				d1 = t1.data[:, i1]
				d2 = t2.data[:, i2]
				# mask out nan
				mask = (numpy.isnan(d1) | numpy.isnan(d2)) ^ True
				corr_res = meth(d1[mask], d2[mask])
				res_list.append(
					CorrRes(
						node1=c1,
						node2=c2,
						r=corr_res.statistic,
						pvalue=corr_res.pvalue
					)
				)

	# fdr control
	pvalue = numpy.asarray([i.pvalue for i in res_list], dtype=float)
	pvalue_argsort = pvalue.argsort()
	sig_flag = numpy.full(len(pvalue), numpy.nan, dtype=float)

	#
	for alpha in sorted(cfg.get("alpha", [0.05]), reverse=True):
		thres = numpy.linspace(0, alpha, len(pvalue) + 1)[1:]
		mask = numpy.full(len(pvalue_argsort), False, dtype=bool)
		mask[pvalue_argsort] = (pvalue[pvalue_argsort] <= thres)
		sig_flag[mask] = alpha

	for res, sig in zip(res_list, sig_flag):
		res.sig_flag = sig

	return res_list


def save_corr_res(cfg: dict, corr_res_list: list):
	prefix = cfg.get("corr_res_prefix", None)
	if not prefix:
		return

	# ris = sorted(set([i[0] for i in iter_corr_ipairs(cfg)]))
	# cis = sorted(set([i[1] for i in iter_corr_ipairs(cfg)]))
	# rows = list(itertools.chain(*[cfg["groups"][i]["table"].cols for i in ris]))
	# cols = list(itertools.chain(*[cfg["groups"][i]["table"].cols for i in cis]))
	rows = list(itertools.chain(*[g["table"].cols for g in cfg["groups"]]))
	cols = rows
	rows_key_map = {v: i for i, v in enumerate(rows)}
	cols_key_map = {v: i for i, v in enumerate(cols)}

	delimiter = cfg.get("delimiter", "\t")

	# with pylib.util.get_fp(fname, "w") as fp:
	for tag, meth in zip(
			["r", "p", "sig_flag"],
			[lambda x: x.r, lambda x: x.pvalue, lambda x: x.sig_flag]
		):
		with pylib.util.get_fp("%s.%s.tsv" % (prefix, tag), "w") as fp:
			# head line
			line = delimiter.join([tag] + cols)
			print(line, file=fp)
			# prepare data
			data = numpy.full((len(rows), len(cols)), numpy.nan, dtype=float)
			for res in corr_res_list:
				rid = rows_key_map[res.node1]
				cid = cols_key_map[res.node2]
				data[rid, cid] = data[cid, rid] = meth(res)  # fancy
			for t, d in zip(rows, data):
				line = delimiter.join([t] + [str(i) for i in d])
				print(line, file=fp)

	return


def create_layout(cfg_list: list) -> dict:
	legend_space = 0.6

	lc = mpllayout.LayoutCreator(
		origin="topleft",
		left_margin=0.1,
		right_margin=0.1,
		top_margin=0.1,
		bottom_margin=0.1 + (legend_space if len(cfg_list) <= 2 else 0),
	)

	chord_axes_size = 3

	for i, cfg in enumerate(cfg_list):
		axes = lc.add_frame(cfg["key"])
		axes.set_size(chord_axes_size, chord_axes_size)

		if i == 0:
			axes.set_anchor("topleft")
		else:
			if i & 0x1:
				offsets = (0.1, 0)
			elif i == 2:
				offsets = (-0.1, -0.1 - legend_space)
			else:
				offsets = (-0.1, -0.1)
			axes.set_anchor("topleft" if (i & 0x1) else "topright",
				last_axes,
				ref_anchor=("topright" if (i & 0x1) else "bottomleft"),
				offsets=offsets,
			)

		last_axes = axes

	# create layout
	layout = lc.create_figure_layout()

	# apply axes style
	for name, axes in layout.items():
		if name == "figure":
			continue

		axes.set_facecolor("#ffffff")
		for sp in axes.spines.values():
			sp.set_visible(False)
		axes.tick_params(
			left=False, labelleft=False,
			right=False, labelright=False,
			top=False, labeltop=False,
			bottom=False, labelbottom=False,
		)

	return layout


def load_list(fname):
	with pylib.util.get_fp(fname, "r") as fp:
		ret = fp.read().splitlines()
	return ret


def load_group_data_inplace(cfg: dict) -> None:
	for g in cfg["groups"]:
		table = pylib.table.RelaAbundTable.from_file(g["file"])
		if "col_select" in g:
			col_select = load_list(g["col_select"])
			table = table.filter_cols(by_tags=col_select)
		g["table"] = table
	return


def _get_arc_steps(theta1, theta2, step_deg=0.02):
	# assumes theta1 < theta2, in degree
	return numpy.linspace(theta1, theta2,
		int(numpy.ceil((theta2 - theta1) / step_deg)),
		endpoint=True,
	)


def _get_annular_sector_path(*, theta1, theta2, r1, r2):
	# assumes theta1 < theta2, in degree
	arc_steps = numpy.deg2rad(_get_arc_steps(theta1, theta2))
	upath = numpy.empty((len(arc_steps), 2), dtype=float)
	upath[:, 0] = numpy.cos(arc_steps) * r1
	upath[:, 1] = numpy.sin(arc_steps) * r1
	dpath = numpy.empty((len(arc_steps), 2), dtype=float)
	dpath[:, 0] = numpy.cos(arc_steps[::-1]) * r2
	dpath[:, 1] = numpy.sin(arc_steps[::-1]) * r2
	return numpy.vstack([upath, dpath])


class ChordWidgetBase(abc.ABC):
	@abc.abstractmethod
	def draw(self, axes: matplotlib.axes.Axes, *ka, **kw):
		pass


class ChordNodeWidgetBase(ChordWidgetBase):
	@abc.abstractmethod
	def place(self, r1, r2, theta1, theta2, **kw):
		pass


class Connectible(object):
	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self.connections = list()
		return

	@property
	def n_connections(self) -> int:
		return len(self.connections)

	@property
	def has_connection(self) -> bool:
		return bool(self.connections)


class Collectible(object):
	def __init__(self, key: str, *ka, label: str = None, **kw):
		super().__init__(*ka, **kw)
		self.key = key
		self.label = key if label is None else label
		self.collection = None
		return

	@property
	def id_in_collection(self) -> int:
		return self._id_in_collection

	def set_id_in_collection(self, v: int):
		self._id_in_collection = v
		return


class Connection(ChordWidgetBase, Collectible):
	def __init__(self, node1: Connectible, node2: Connectible,
			*ka, data: CorrRes = None, **kw):
		key = (":").join([node1.key, node2.key])
		super().__init__(key, *ka, **kw)
		for s, n in zip(["node1", "node2"], [node1, node2]):
			if not isinstance(n, Connectible):
				raise TypeError("%s must be an instance of %s"
					% (s, Connectible.__name__))
		self.node1 = node1
		self.node2 = node2
		self.data = data
		self._term_thetas = dict()
		return

	@staticmethod
	def is_insig_static(sig_flag: float) -> bool:
		return numpy.isnan(sig_flag)

	@property
	def is_insig(self) -> bool:
		return self.is_insig_static(self.data.sig_flag)

	def share_node_with(self, other) -> bool:
		# test if this connection has a common node with another connection
		return ((self.node1 is other.node1) or (self.node1 is other.node2)
			or (self.node2 is other.node1) or (self.node2 is other.node2))

	def connects(self, node) -> bool:
		# test if the argument node is related to this connection
		return ((node is self.node1) or (node is self.node2))

	def crosses(self, other) -> bool:
		# test if this connection crosses over another connection
		if self.share_node_with(other):
			return False
		nodes = sorted([self.node1, self.node2, other.node1, other.node2],
			key=lambda x: x.global_node_id)
		# there is a crossover if (0, 2) belongs to the same connection
		# otherwise no
		return ((self.connects(nodes[0]) and self.connects(nodes[2])) or
			(other.connects(nodes[0]) and other.connects(nodes[2])))

	def set_term_thetas(self, node: ChordNodeWidgetBase,
			*, theta1: float, theta2: float):
		self._term_thetas[node.key] = (min(theta1, theta2), max(theta1, theta2))
		return

	def get_term_thetas(self, node: ChordNodeWidgetBase) -> (float, float):
		return self._term_thetas[node.key]

	@property
	def span_reserve(self) -> float:
		return 0.5 + 9.5 * ((abs(self.data.r) - 0.25) / 0.75) ** 2

	def init_draw_theta(self):
		self.node1.draw_theta = self.node1.theta1
		self.node2.draw_theta = self.node2.theta1
		return

	def _get_term_path_segment(self, node, *, end_div: float = 0.03,
			neighbor_div: float = 0.1) -> numpy.ndarray:
		base_radius = node.r1 - end_div
		theta1, theta2 = self.get_term_thetas(node)
		# draw a division between chords
		theta1 += neighbor_div
		theta2 -= neighbor_div

		if self.data.r > 0:
			# if corr coef r > 0 draw full ends
			end_theta = _get_arc_steps(theta1, theta2)
			radius = base_radius
		else:
			# if corr coef r < 0 draw notched ends
			end_theta = [
				theta1,
				(theta1 + theta2) / 2,
				theta2,
			]
			radius = [base_radius, base_radius * 0.95, base_radius]

		# the xy
		xy = numpy.asarray(
			[
				numpy.cos(numpy.deg2rad(end_theta)) * radius,
				numpy.sin(numpy.deg2rad(end_theta)) * radius,
			], dtype=float
		).T
		return xy

	@staticmethod
	def sig_to_alpha(sig: float, damping=1.0) -> str:
		if Connection.is_insig_static(sig):
			return "60"
		a = int((int((1 - min(1, sig / 0.05)) * 0x40) + 0x40) * damping)
		return "%02x" % min(a, 0xff)

	# class variable
	insig_rgb = "#c0c0c0"
	highr_insig_rgb = "#404040"
	outline_minr = 0.80

	@property
	def edgecolor_rgb(self) -> str:
		if abs(self.data.r) > self.outline_minr:
			if self.is_insig:
				return self.highr_insig_rgb
			else:
				return self.node1.color
		return "none"

	@property
	def facecolor_rgba(self) -> str:
		if numpy.isnan(self.data.sig_flag):
			return self.insig_rgb + self.sig_to_alpha(self.data.sig_flag)
		return self.node1.color + self.sig_to_alpha(self.data.sig_flag)

	def draw(self, axes: matplotlib.axes.Axes, zorder=None):
		codes = list()

		# node1 arc
		node1_arc = self._get_term_path_segment(self.node1)
		codes += [matplotlib.path.Path.MOVETO] \
			+ ([matplotlib.path.Path.LINETO] * (len(node1_arc) - 1))

		# node1->node2 bezier
		bezier1 = numpy.array([[0, 0]], dtype=float)
		codes += [matplotlib.path.Path.CURVE3]

		# node2 arc
		node2_arc = self._get_term_path_segment(self.node2)
		codes += [matplotlib.path.Path.CURVE3] \
			+ ([matplotlib.path.Path.LINETO] * (len(node2_arc) - 1))

		# node2->node1 bezier
		bezier2 = numpy.array([[0, 0], node1_arc[0]], dtype=float)
		codes += [matplotlib.path.Path.CURVE3, matplotlib.path.Path.CURVE3]

		# make path
		vertices = numpy.vstack([node1_arc, bezier1, node2_arc, bezier2])
		path = matplotlib.path.Path(vertices, codes, closed=True)

		# color = (self.node1.color if self.data.sig_flag else "#c0c0c0")
		patch = matplotlib.patches.PathPatch(path, linewidth=0.5,
			edgecolor=self.edgecolor_rgb, facecolor=self.facecolor_rgba,
			zorder=zorder
				if ((zorder is None) or (numpy.isnan(self.data.sig_flag)))
				else zorder + 1
		)
		axes.add_patch(patch)

		return


class Collection(object):
	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self._ele_list = list()
		self._ele_dict = dict()
		return

	# when using subscription, treat as a list
	def __getitem__(self, *ka):
		return self._ele_list.__getitem__(*ka)

	def __setitem__(self, *ka):
		return self._ele_list.__setitem__(*ka)

	@functools.wraps(list.sort)
	def sort(self, **kw):
		self._ele_list.sort(**kw)
		# reapply ids
		for i, v in enumerate(self):
			v.set_id_in_collection(i)
		return

	def filter(self, cond_func):
		"""
		filter elements in current collection, keep those cond_func returns True
		and remove those cond_func returns False/None; then update the
		id_in_collection attributes of kept elements, and remove the collection
		binds with the removed elements

		cond_func must be a callable with signature cond_func(element)
		"""
		remove = [v for v in self if not cond_func(v)]
		keep = [v for v in self if cond_func(v)]
		# update element container
		self._ele_list = keep
		self._ele_dict = {v.key: v for v in keep}
		# update id_in_collection
		for i, v in enumerate(self):
			v.set_id_in_collection(i)
		for v in remove:
			v.collection = None
		return

	def swap_by_id(self, id1: int, id2: int):
		if id1 == id2:
			raise ValueError("two ids cannot be identical")
		self[id1], self[id2] = self[id2], self[id1]
		self[id1].set_id_in_collection(id1)
		self[id2].set_id_in_collection(id2)
		return

	def __len__(self) -> int:
		return len(self._ele_list)

	def __iter__(self):
		return iter(self._ele_list)

	def get(self, *, key: str = None, id: int = None):
		if not ((key is None) ^ (id is None)):
			raise ValueError("must provide one of key or id, but not both")
		if id is not None:
			return self._ele_list[id]
		elif key is not None:
			return self._ele_dict[key]

	def add(self, element: Collectible) -> Collectible:
		if not isinstance(element, Collectible):
			raise TypeError("added element must be an instance of "
				"Collectible")
		if element.key in self._ele_dict:
			raise ValueError("key '%s' alredy exists" % key)
		element.set_id_in_collection(len(self))
		element.collection = self
		self._ele_list.append(element)
		self._ele_dict[element.key] = element
		return element

	def has_key(self, key: str):
		return key in self._ele_dict


class ChordNode(ChordNodeWidgetBase, Collectible, Connectible):
	def __init__(self, key, *ka, color="#c0c0c0", **kw):
		super().__init__(key, *ka, **kw)
		self.color = color
		return

	@property
	def global_node_id(self) -> int:
		parent = self.collection
		node_id = 0
		# sum the num of nodes in previous groups
		for g in parent.collection[:parent.id_in_collection]:
			node_id += len(g)
		return node_id + self.id_in_collection

	def place(self, *, r1, r2, theta1, theta2, unit_span):
		self.r1 = min(r1, r2)
		self.r2 = max(r1, r2)
		while theta2 < theta1:
			theta2 += 360
		self.theta1 = theta1
		self.theta2 = theta2
		self.unit_span = unit_span
		return

	def get_conn_remote_node(self, conn: Connection) -> ChordNodeWidgetBase:
		if not conn.connects(self):
			raise ValueError("connection must be related to current node")
		return (conn.node2 if self is conn.node1 else conn.node1)

	def get_mean_conn_niche(self) -> float:
		niches = list()
		self_group_niche = self.collection.id_in_collection
		for conn in self.connections:
			niche = self.get_conn_remote_node(conn).collection.id_in_collection
			if niche < self_group_niche:
				niche = len(self.collection.collection) - niche
			niches.append(niche)
		if not niches:
			return self_group_niche + 1
		else:
			return sum(niches) / len(niches)

	def _apply_connection_coords(self):
		curr = self.theta1
		for i, conn in enumerate(self.connections):
			span = conn.span_reserve * self.unit_span
			conn.set_term_thetas(self, theta1=curr, theta2=curr + span)
			curr += span
		return

	def sort_connections(self):
		self.connections.sort(key=lambda x:
				(self.get_conn_remote_node(x).theta1 - self.theta1) % 360,
			reverse=True
		)
		self._apply_connection_coords()
		return

	@property
	def span_reserve(self) -> float:
		return max(0.5, sum([i.span_reserve for i in self.connections]))

	def draw(self, axes: matplotlib.axes.Axes, zorder=None):
		path = _get_annular_sector_path(theta1=self.theta1, theta2=self.theta2,
			r1=self.r1, r2=self.r2)
		polygon = matplotlib.patches.Polygon(path, closed=True,
			facecolor=self.color, zorder=zorder,
		)
		axes.add_patch(polygon)
		self.patch = polygon
		# draw label
		rotation = ((self.theta1 + self.theta2) / 2) % 360
		theta = numpy.deg2rad(rotation)
		if self.collection.text_orient_v == "top":
			r_offset = 0.10
			rotation -= 90
		else:
			r_offset = 0.14
			rotation += 90

		text = axes.text(
			x=numpy.cos(theta) * (self.r2 + r_offset),
			y=numpy.sin(theta) * (self.r2 + r_offset),
			s=" " + self.label + " ",
			fontsize=4,
			color=self.color,
			rotation=rotation,
			rotation_mode="anchor",
			horizontalalignment="center",
			verticalalignment="center",
			zorder=(None if zorder is None else zorder + 1),
		)

		return


class ChordNodeGroup(ChordNodeWidgetBase, Collectible, Collection):
	def __init__(self, key, *ka, color="#e0e0e0", **kw):
		super().__init__(key, *ka, **kw)
		self.color = color
		return

	@property
	def span_reserve(self) -> float:
		return sum([i.span_reserve for i in self])

	def add(self, key, *ka, **kw):
		element = ChordNode(key, *ka, **kw)
		return super().add(element)

	def sort_nodes_by_proximity(self):
		# preliminarily sort by connection prximity
		self.sort(key=lambda x: x.get_mean_conn_niche(), reverse=True)
		return

	def _count_cross_change_by_swap(self, id1: int, id2: int) -> int:
		if id1 == id2:
			raise ValueError("two ids cannot be identical")
		# ensure id1 < id2
		if id1 > id2:
			id1, id2 = id2, id1

		count = 0
		for i in (id1, id2):
			test_node = self[i]
			# checking only the nodes between the proposed two is enough
			for inter_node in self[id1 + 1:id2 + 1]:
				# checking {conn from test node} x {conn from nodes in between}
				for c1, c2 in itertools.product(test_node.connections,
						inter_node.connections):
					if c1.share_node_with(c2):
						continue
					elif c1.crosses(c2):
						# if currently crossed, will be uncrossed after swap
						count -= 1
					else:
						# if currently uncrossed, will be crossed after swap
						count += 1
		return count

	def sort_nodes_by_cross(self):
		should_continue = True
		while should_continue:
			# fine tuning by reducing connection crosses
			id1, id2, count = None, None, 0
			for i in range(len(self)):
				for j in range(i + 1, len(self)):
					test_count = self._count_cross_change_by_swap(i, j)
					if (test_count < count):
						id1, id2, count = i, j, test_count
			#
			should_continue = count != 0
			if should_continue:
				self.swap_by_id(id1, id2)
		return

	def sort_nodes(self):
		self.sort_nodes_by_proximity()
		self.sort_nodes_by_cross()
		return

	def place(self, *, r1, r2, theta1, theta2, group_padding=0, node_gap=0,
			r1_indent=0, r2_indent=0):
		self.r1 = min(r1, r2)
		self.r2 = max(r1, r2)

		# ensure theta2 > theta1
		while theta2 < theta1:
			theta2 += 360
		self.theta1 = theta1
		self.theta2 = theta2

		self.group_padding = group_padding
		self.node_gap = node_gap

		# resolve the placement of nodes
		unit_span = (self.theta2 - self.theta1 - 2 * group_padding -
			(len(self) - 1) * node_gap) / self.span_reserve
		curr = self.theta1 + group_padding
		for i, node in enumerate(self):
			node_span = unit_span * node.span_reserve
			node.place(
				r1=self.r1 + r1_indent,
				r2=self.r2 - r2_indent,
				theta1=curr,
				theta2=curr + node_span,
				unit_span=unit_span,
			)
			curr += node_span + node_gap
		return

	@property
	def theta(self) -> bool:
		return (self.theta1 + self.theta2) / 2 % 360

	def _draw_outer_arc_label(self, axes: matplotlib.axes.Axes, zorder=None):
		fontsize = 8
		char_deg_base = fontsize / 2.7
		label = self.label
		r_base = self.r2 * 1.15

		# proposed degree position of each character
		char_deg_arr = ((numpy.arange(len(label)) - (len(label) - 1) / 2)
			* char_deg_base + self.theta)

		# text orientation
		if self.text_orient_v == "top":
			r_offset = 0.01
			rotation_offset = -90
			char_deg_arr = char_deg_arr[::-1]
		else:
			r_offset = 0.06
			rotation_offset = +90

		# first draw base, then character by character overlay
		for layer in range(2):
			path_effects = ([matplotlib.patheffects.withStroke(linewidth=0,
				foreground="#ffffff")] if layer == 0 else None)
			for c, d in zip(label, char_deg_arr):
				axes.text(
					x=numpy.cos(numpy.deg2rad(d)) * (r_base + r_offset),
					y=numpy.sin(numpy.deg2rad(d)) * (r_base + r_offset),
					s=(" " + c + " "),  # add two spaces to align some chars
					fontfamily="monospace",
					fontsize=fontsize,
					color=("#ffffff" if layer == 0 else self.color),
					rotation=d + rotation_offset,
					rotation_mode="anchor",
					horizontalalignment="center",
					verticalalignment="center",
					path_effects=path_effects,
					zorder=(None if zorder is None else zorder + 1 + layer),
				)

		return

	def draw(self, axes: matplotlib.axes.Axes, zorder=None):
		# inner arc
		path = _get_annular_sector_path(theta1=self.theta1, theta2=self.theta2,
			r1=self.r1, r2=self.r2)
		polygon = matplotlib.patches.Polygon(path, closed=True,
			edgecolor="none", facecolor=self.color + "60",
			zorder=zorder,
		)
		axes.add_patch(polygon)
		self.patch = polygon

		# outer arc
		arc = matplotlib.patches.Arc((0, 0),
			width=self.r2 * 2.2,
			height=self.r2 * 2.2,
			theta1=self.theta1 + self.group_padding,
			theta2=self.theta2 - self.group_padding,
			edgecolor=self.color,
		)
		axes.add_patch(arc)

		self._draw_outer_arc_label(axes, zorder=zorder)

		# child nodes
		for n in self:
			n.draw(axes, zorder=(None if zorder is None else zorder + 2))
		return

	@property
	def text_orient_v(self) -> str:
		if (((self.theta >= 0) and (self.theta <= 200)) or (self.theta >= 340)):
			return "top"
		else:
			return "bottom"


class ChordDiagram(object):
	def __init__(self, *ka, **kw):
		self.node_groups = self.NodeGroups()
		self.chords = self.ChordConnections()
		return

	def add_chord(self, key1: str, key2: str, *, data=None) -> Connection:
		node1 = self.get_node_by_key(key1)
		node2 = self.get_node_by_key(key2)
		new = self.chords.add(node1, node2, data=data)
		node1.connections.append(new)
		node2.connections.append(new)
		return new

	def draw(self, axes: matplotlib.axes.Axes, zorder=None):
		self.node_groups.draw(axes, zorder=zorder)
		self.chords.draw(axes, zorder=zorder)
		return

	def get_node_by_key(self, key: str):
		for g in self.node_groups:
			if g.has_key(key):
				return g.get(key=key)
		raise KeyError("no node of key '%s' found" % key)
		return

	class NodeGroups(ChordNodeWidgetBase, Collection):
		def add(self, key, *ka, **kw):
			element = ChordNodeGroup(key, *ka, **kw)
			return super().add(element)

		def filter_nodes(self, cond_func):
			for g in self:
				g.filter(cond_func)
			return

		def place(self, *, span=360, theta_offset=0, r1=1, r2=1.1, group_gap=0,
				group_padding=0, node_gap=0, node_r1_indent=0,
				node_r2_indent=0):
			# arrange nodes in each group
			self.sort_nodes()

			# total gaps (white space)
			n_node = sum([len(i) for i in self])
			total_gap = group_gap * (0 if len(self) == 1 else len(self)) \
				+ group_padding * (0 if len(self) == 1 else len(self) * 2) \
				+ node_gap * (n_node - (0 if len(self) == 1 else len(self)))

			# total units of node span reserve
			node_span_reserve = sum([i.span_reserve for i in self])

			# unit size of span reserve
			unit_span = (span - total_gap) / node_span_reserve

			curr = theta_offset
			for g in self:
				if len(self) == 1:
					group_span = span
				else:
					group_span = unit_span * g.span_reserve \
						+ node_gap * (len(g) - 1) \
						+ group_padding * 2
				g.place(r1=r1, r2=r2, theta1=curr, theta2=curr + group_span,
					group_padding=group_padding, node_gap=node_gap,
					r1_indent=node_r1_indent, r2_indent=node_r2_indent)
				curr += group_span + group_gap

			return

		def sort_nodes(self):
			for g in self:
				g.sort_nodes()
			return

		def sort_connections(self):
			for g in self:
				for n in g:
					n.sort_connections()
			return

		def draw(self, axes: matplotlib.axes.Axes, zorder=None):
			for i in self:
				i.draw(axes, zorder=(None if zorder is None else zorder + 2))
			return

	class ChordConnections(ChordWidgetBase, Collection):
		def add(self, *ka, **kw):
			conn = Connection(*ka, **kw)
			return super().add(conn)

		def draw(self, axes: matplotlib.axes.Axes, zorder=None):
			for i in self:
				i.init_draw_theta()
			for i in self:
				i.draw(axes, zorder=(None if zorder is None else zorder + 2))
			return

	def place(self, *ka, anchor=None, anchor_angle=0, **kw):
		# try put the node groups to adjust the angle of anchor, if necessary
		self.node_groups.place(*ka, theta_offset=0, **kw)
		if anchor is not None:
			theta_offset = anchor_angle - self.node_groups.get(key=anchor).theta
			self.node_groups.place(*ka, theta_offset=theta_offset, **kw)
		self.node_groups.sort_connections()
		return


def plot_corr_chord(axes: matplotlib.axes.Axes, cfg: dict, corr_res_list: list):
	# make chord diagram
	chord = ChordDiagram()

	# add groups and nodes
	for g in cfg["groups"]:
		col_translate = load_col_translate_safe(g)
		group = chord.node_groups.add(g["key"], label=g["display"],
			color=g["color"])
		for oligo in g["table"].cols:
			group.add(oligo, label=col_translate.get(oligo, oligo),
				color=g["color"])

	# add chords
	min_r = cfg.get("min_r", 0.3)
	for res in corr_res_list:
		if abs(res.r) < min_r:
			continue
		chord.add_chord(res.node1, res.node2, data=res)

	# filter nodes
	chord.node_groups.filter_nodes(lambda node: node.has_connection)

	# rearrange ndoes and connections, and calculate the coordinates of widgets
	chord.place(r1=4.15, r2=4.20, anchor=cfg.get("anchor", None),
		anchor_angle=cfg.get("anchor_angle", 0),
		group_gap=0.0, group_padding=4.0, node_gap=2.0,
		node_r1_indent=-0.03, node_r2_indent=0.00)

	# draw the chord diagram
	chord.draw(axes, zorder=2)

	# misc
	lim_radius = 5.0
	axes.set_xlim(-lim_radius, lim_radius)
	axes.set_ylim(-lim_radius, lim_radius)

	return


def gen_legend_chord_handle_path(xdescent, ydescent, width, height, *,
		notched=False, dep=0.1):
	xy = list()
	ts = numpy.linspace(-1, 1, 128, endpoint=True)
	# lower arc
	for t in ts:
		x = -xdescent + width * (1 + t) * 0.5
		y = -ydescent + height * (1 - t ** 2) * dep
		xy.append((x, y))
	# right notch
	if notched:
		xy.append((-xdescent + width * (1 - dep), -ydescent + height * 0.5))
	# upper arc
	for t in ts:
		x = -xdescent + width * (1 - t) * 0.5
		y = -ydescent + height * (1 - dep + t ** 2 * dep)
		xy.append((x, y))
	# left notch
	if notched:
		xy.append((-xdescent + width * dep, -ydescent + height * 0.5))

	return numpy.asarray(xy, dtype=float)


class DummyArtistNormChordPatch(matplotlib.patches.Patch):
	class LegendHandler(matplotlib.legend_handler.HandlerPatch):
		def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
				height, fontsize, trans):
			path = gen_legend_chord_handle_path(xdescent, ydescent, width,
				height, notched=False)
			p = matplotlib.patches.Polygon(path, closed=True)
			self.update_prop(p, orig_handle, legend)
			p.set_transform(trans)
			return [p]


class DummyArtistNotchChordPatch(matplotlib.patches.Patch):
	class LegendHandler(matplotlib.legend_handler.HandlerPatch):
		def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
				height, fontsize, trans):
			path = gen_legend_chord_handle_path(xdescent, ydescent, width,
				height, notched=True)
			p = matplotlib.patches.Polygon(path, closed=True)
			self.update_prop(p, orig_handle, legend)
			p.set_transform(trans)
			return [p]


class DummyArtistMultiColorPatch(matplotlib.artist.Artist):
	def __init__(self, edgecolors: list, facecolors: list, *ka, label=None,
			**kw):
		if len(edgecolors) != len(facecolors):
			raise ValueError("edgecolors and facecolors must have matching "
				"length")
		super().__init__(*ka, **kw)
		self.edgecolors = edgecolors
		self.facecolors = facecolors
		self.set_label(label)
		return

	class LegendHandler(matplotlib.legend_handler.HandlerPatch):
		def __update_prop_func(self, legend_handler, orig_handle):
			return

		def __init__(self, *ka, update_func=None, **kw):
			super().__init__(*ka, update_func=self.__update_prop_func, **kw)
			return

		def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
				height, fontsize, trans):
			# break drawing area into segments
			rela_gap = 0.1
			n_segments = len(orig_handle.edgecolors)
			seg_width = width / (n_segments + (n_segments - 1) * rela_gap)
			seg_step = seg_width * (1 + rela_gap)

			patches = list()
			for i in range(n_segments):
				p = matplotlib.patches.Rectangle(
					(-xdescent + seg_step * i, -ydescent),
					seg_width, height,
					edgecolor=orig_handle.edgecolors[i],
					facecolor=orig_handle.facecolors[i],
				)
				self.update_prop(p, orig_handle, legend)
				p.set_transform(trans)
				patches.append(p)

			return patches


def draw_legend(cfg_list: list, axes: matplotlib.axes.Axes, **kw):
	# group colors
	group_colors = set()
	for cfg in cfg_list:
		for g, _ in iter_corr_pairs(cfg):
			group_colors.add(g["color"])
	group_colors = sorted(group_colors)

	# sig levels
	sig_levels = set()
	for cfg in cfg_list:
		sig_levels.update(cfg["alpha"])
	sig_levels = sorted(sig_levels)

	# draw legend
	handles = list()

	# add positive/negative corr patches
	chord_handle_props = dict(linestyle="-", linewidth=0.5, edgecolor="none",
		facecolor="#a0a0a0")
	handles.append(DummyArtistNormChordPatch(**chord_handle_props,
		label="Pos. correlation"))
	handles.append(DummyArtistNotchChordPatch(**chord_handle_props,
		label="Neg. correlation"))

	# add sig level color patches
	for i in sig_levels:
		artist = DummyArtistMultiColorPatch(
			edgecolors=["none"] * len(group_colors),
			facecolors=[c + Connection.sig_to_alpha(i, damping=1.5)
				for c in group_colors],
			label=r"$\alpha$=%.2f (FDR)" % i,
		)
		handles.append(artist)

	# add high r color patch
	handles.append(DummyArtistNormChordPatch(facecolor="none",
		linewidth=0.5, edgecolor=Connection.highr_insig_rgb,
		label="$|r|$ > %s" % str(Connection.outline_minr))
	)

	# add insig color patch
	handles.append(
		DummyArtistNormChordPatch(edgecolor="none",
			facecolor=Connection.insig_rgb + Connection.sig_to_alpha(numpy.nan),
			label="Insignificant")
	)

	handler_map_override = {
		DummyArtistNormChordPatch: DummyArtistNormChordPatch.LegendHandler(),
		DummyArtistNotchChordPatch: DummyArtistNotchChordPatch.LegendHandler(),
		DummyArtistMultiColorPatch: DummyArtistMultiColorPatch.LegendHandler(),
	}
	handler_map = collections.ChainMap(handler_map_override,
		matplotlib.legend.Legend.get_default_handler_map()
	)

	return axes.legend(handles=handles, handler_map=handler_map,
		ncols=3, **kw)


def plot(png, cfg_list: list, *, dpi=300):
	plot_cfg_list = [i for i in cfg_list if i.get("enable", True)]
	if not plot_cfg_list:
		return

	# create layout
	layout = create_layout(plot_cfg_list)
	figure = layout["figure"]

	for cfg in plot_cfg_list:
		load_group_data_inplace(cfg)
		corr_res_list = corr_analysis_with_fdr_control(cfg)
		save_corr_res(cfg, corr_res_list)

		axes = layout[cfg["key"]]
		plot_corr_chord(axes, cfg, corr_res_list)

	axes = layout[plot_cfg_list[0]["key"]]
	draw_legend(plot_cfg_list, axes=axes, frameon=False, fontsize=8, loc=9,
		bbox_to_anchor=(1.02, -0.02), handlelength=2.5)

	# savefig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close()
	return


def main():
	plot("manuscript.fig4.png", CHORD_PLOT_CFG, dpi=600)
	return


if __name__ == "__main__":
	main()
