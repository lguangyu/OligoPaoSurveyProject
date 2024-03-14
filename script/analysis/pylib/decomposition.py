#!/usr/bin/env python3

import abc
import functools
import sklearn.base
import sklearn.decomposition
import sklearn.manifold
import umap
from . import method_registry


class DecompositionMethod(sklearn.base.TransformerMixin,
		sklearn.base.BaseEstimator, abc.ABC):
	@property
	@abc.abstractmethod
	def xlabel_str(self) -> str:
		pass

	@property
	@abc.abstractmethod
	def ylabel_str(self) -> str:
		pass


DECOMP_METHODS = method_registry.MethodRegistry(
	value_type = DecompositionMethod)


@DECOMP_METHODS.register("pca")
class PCA(sklearn.decomposition.PCA, DecompositionMethod):
	@property
	def xlabel_str(self):
		return "PC1 (%.1f%%)" % (self.explained_variance_ratio_[0] * 100)

	@property
	def ylabel_str(self):
		return "PC1 (%.1f%%)" % (self.explained_variance_ratio_[1] * 100)


@DECOMP_METHODS.register("t-sne", as_default = True)
class TSNE(sklearn.manifold.TSNE, DecompositionMethod):
	@functools.wraps(sklearn.manifold.TSNE.__init__)
	def __init__(self, *, learning_rate = "auto", init = "pca", **kw):
		super(TSNE, self).__init__(learning_rate = learning_rate, init = init,
			**kw)
		return

	@property
	def xlabel_str(self):
		return "t-SNE 1"

	@property
	def ylabel_str(self):
		return "t-SNE 2"


@DECOMP_METHODS.register("mds")
class MDS(sklearn.manifold.MDS, DecompositionMethod):
	@property
	def xlabel_str(self):
		return "MDS1"

	@property
	def ylabel_str(self):
		return "MDS2"


@DECOMP_METHODS.register("umap")
class UMAP(umap.UMAP, DecompositionMethod):
	@property
	def xlabel_str(self):
		return "UMAP1"

	@property
	def ylabel_str(self):
		return "UMAP2"
