#!/usr/bin/env python3

import abc
import functools
import numpy
from . import method_registry


class DiversityMethod(object):
	@abc.abstractmethod
	def __call__(self, X, *ka, **kw) -> numpy.ndarray:
		pass
	@property
	@abc.abstractmethod
	def name_str(self) -> str:
		pass


DIVERSITY_FUNC = method_registry.MethodRegistry(value_type = DiversityMethod)


@DIVERSITY_FUNC.register("shannon", as_default = True)
class ShannonIndexDiversity(DiversityMethod):
	def __call__(self, X, *ka, pseudo_count = 0.0001, keepdims = False, **kw):
		assert X.ndim == 2
		count = X + pseudo_count
		frac = count / count.sum(axis = 1, keepdims = True)
		return -(frac * numpy.log(frac)).sum(axis = 1, keepdims = keepdims)
	@property
	def name_str(self):
		return "shannon index"


@DIVERSITY_FUNC.register("simpson")
class SimpsonIndexDiversity(DiversityMethod):
	def __call__(self, X, *ka, keepdims = False, **kw):
		assert X.ndim == 2
		frac = X / X.sum(axis = 1, keepdims = True)
		return numpy.power(frac, 2).sum(axis = 1, keepdims = keepdims)
	@property
	def name_str(self):
		return "Simpson's index"
