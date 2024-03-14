#!/usr/bin/env python3

import abc
import functools
import numpy
import sklearn.feature_selection
from . import method_registry


class RelevanceMethod(object):
	@abc.abstractmethod
	def __call__(self, X, y = None, *ka, **kw) -> numpy.ndarray:
		pass
	@property
	@abc.abstractmethod
	def name_str(self) -> str:
		pass


RELEVANCE_FUNC = method_registry.MethodRegistry(value_type = RelevanceMethod)


@RELEVANCE_FUNC.register("mi", as_default = True)
class NMIRelevance(RelevanceMethod):
	@functools.wraps(sklearn.feature_selection.mutual_info_classif)
	def __call__(self, *ka, **kw):
		all_kw = dict(discrete_features = False, n_neighbors = 10)
		all_kw.update(kw)
		return sklearn.feature_selection.mutual_info_classif(*ka, **all_kw)
	@property
	def name_str(self):
		return "mutual information"


@RELEVANCE_FUNC.register("anova-f")
class ANOVAF(RelevanceMethod):
	@functools.wraps(sklearn.feature_selection.f_classif)
	def __call__(self, *ka, **kw):
		f, p = sklearn.feature_selection.f_classif(*ka, **kw)
		# remove all non-significant one
		f[p >= 0.05] = 0
		return f
	@property
	def name_str(self):
		return "ANOVA F-value"
