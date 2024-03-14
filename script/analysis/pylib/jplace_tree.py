#!/usr/bin/env python3

import dataclasses
import re
import typing
from . import util


class JTree(object):
	@dataclasses.dataclass
	class JTreeNode(object):
		name: typing.Optional[str] = None
		dist: util.NonNegFloat = 0
		edge_id: typing.Optional[int] = None
		children: list = dataclasses.field(default_factory=lambda: list())
		_parent = None

		def __hash__(self) -> int:
			return id(self)

		@property
		def parent(self):
			return self._parent

		@property
		def is_root(self) -> bool:
			return self._parent is None

		def add_child(self, node):
			self.children.append(node)
			node._parent = self
			return

		@classmethod
		def parse(cls, s: str):
			new = cls()
			last_comma = 0
			stack = 0
			for i, c in enumerate(s):
				if c == "(":
					stack += 1
				elif ((c == ",") or c == ")") and (stack == 1):
					new.add_child(cls.parse(s[last_comma + 1: i]))
					last_comma = i

				if c == ")":
					stack -= 1
				if stack == 0:
					break

			if c == ")":
				i += 1

			m = re.match("(.*):(.*)\{(.*)\};?$", s[i:])
			name, dist, edge_id = m.groups()

			new.name = name if name else None
			new.dist = float(dist)
			new.edge_id = int(edge_id)

			return new

		def traverse(self):
			yield self
			for c in self.children:
				yield from c.traverse()
			return

		@property
		def n_leaves(self) -> int:
			return sum([node.is_leaf for node in self.traverse()])

		@property
		def is_leaf(self):
			return not (self.children)

		def path_to_child(self, node) -> typing.Optional[list]:
			path = [node]
			if node is self:
				return path
			while node.parent is not None:
				path.append(node.parent)
				if node.parent is self:
					return path[::-1]
				node = node.parent
			else:
				return None
			return

		def dist_to_child(self, node) -> typing.Optional[float]:
			path = self.path_to_child(node)
			if path is None:
				return None
			dist = 0.
			for n in path[1:]:
				dist += n.dist
			return dist

		def reroot(self):
			pre_parent = self.parent
			if pre_parent is None:
				return
			if pre_parent.parent is not None:
				pre_parent.reroot()
			pre_parent.children = [i for i in self.parent.children
				if i is not self]
			pre_parent._parent = self
			pre_parent.dist = self.dist
			self._parent = None
			self.children.append(pre_parent)
			self.dist = 0
			return

	def __init__(self, root: JTreeNode, *ka, **kw):
		super().__init__(*ka, **kw)
		self.root = root
		self.renew_map_name_to_node()
		self.renew_map_edge_to_node()
		return

	def renew_map_name_to_node(self):
		d = dict()
		for node in self.traverse():
			if node.name:
				d[node.name] = node
		self._name_to_node = d
		return

	def renew_map_edge_to_node(self):
		d = dict()
		for node in self.traverse():
			d[node.edge_id] = node
		self._edge_to_node = d
		return

	def get_node_by_name(self, name: str) -> JTreeNode:
		return self._name_to_node[name]

	def get_node_by_edge_id(self, edge_id: int) -> JTreeNode:
		return self._edge_to_node[edge_id]

	@classmethod
	def parse(cls, s: str, *, node_cls=JTreeNode):
		new = cls(root=node_cls.parse(s))
		return new

	def traverse(self):
		yield from self.root.traverse()

	def get_path_to_node(self, node: JTreeNode) -> typing.Optional[list]:
		return self.root.path_to_child(node)

	def lowest_common_ancestor(self, *nodes) -> JTreeNode:
		paths = [self.get_path_to_node(i) for i in nodes]
		if any([i is None for i in paths]):
			raise ValueError("all nodes must be from the same tree")
		ret = self.root
		for i in zip(*paths):
			uniques = set(i)
			if len(uniques) != 1:
				break
			ret = uniques.pop()
		return ret

	def reroot(self, node):
		node.reroot()
		self.root = node
		return
