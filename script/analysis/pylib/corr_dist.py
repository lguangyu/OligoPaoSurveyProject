#!/usr/bin/env python3

import abc
import functools
import numpy
import sklearn.metrics
from . import method_registry


class CorrDistMethod(object):
	@abc.abstractmethod
	def __call__(self, X, Y = None, *ka, **kw) -> numpy.ndarray:
		pass
	@property
	@abc.abstractmethod
	def name_str(self) -> str:
		pass
	def to_plot_data(self, dist: numpy.ndarray) -> numpy.ndarray:
		return dist
	@property
	def vmax(self):
		return None
	@property
	def vmin(self) -> float:
		return 0
	@property
	def cmap(self):
		return "BuPu"


CORR_DIST_FUNC = method_registry.MethodRegistry(value_type = CorrDistMethod)


@CORR_DIST_FUNC.register("euclidean")
class EuclideanDist(CorrDistMethod):
	@functools.wraps(sklearn.metrics.pairwise.euclidean_distances)
	def __call__(self, *ka, **kw):
		return sklearn.metrics.pairwise.euclidean_distances(*ka, **kw)
	@property
	def name_str(self):
		return "Euclidean distance"


@CORR_DIST_FUNC.register("cosine", as_default = True)
class CosineDist(CorrDistMethod):
	@functools.wraps(sklearn.metrics.pairwise.cosine_distances)
	def __call__(self, *ka, **kw):
		return sklearn.metrics.pairwise.cosine_distances(*ka, **kw)
	@property
	def name_str(self):
		return "cosine similarity"
	def to_plot_data(self, dist):
		return 1 - dist
	@property
	def vmin(self):
		return -1
	@property
	def vmax(self):
		return 1
	@property
	def cmap(self):
		return "RdYlBu_r"


@CORR_DIST_FUNC.register("sqrt_cosine")
class SqrtCosineDist(CorrDistMethod):
	@functools.wraps(sklearn.metrics.pairwise.cosine_distances)
	def __call__(self, *ka, **kw):
		return numpy.sqrt(sklearn.metrics.pairwise.cosine_distances(*ka, **kw))
	@property
	def name_str(self):
		return "Sqrt. cosine distance"
	@property
	def vmax(self):
		return 1.4145
	@property
	def cmap(self):
		return "RdYlBu_r"
