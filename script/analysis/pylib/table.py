#!/usr/bin/env python3

import collections
import numpy
import skbio
from . import util


class Table(object):
	def __init__(self, data, rows, cols):
		self.set_all_data(data, rows, cols)
		return

	def set_all_data(self, data, rows, cols):
		data = numpy.asarray(data, dtype=float)
		if data.ndim != 2:
			raise ValueError("data must be a 2-d array")
		if data.shape != (len(rows), len(cols)):
			raise ValueError("shape of data does not comply with provided rows "
				"and/or cols: data (%u x %u), rows (%u), cols (%u)" %
				(*data.shape, len(rows), len(cols))
			)
		self.data = data
		self.rows = numpy.asarray(rows, dtype=object)
		self.cols = numpy.asarray(cols, dtype=object)
		self.rows_id_map = {v: i for i, v in enumerate(rows)}
		self.cols_id_map = {v: i for i, v in enumerate(cols)}
		return

	@property
	def n_rows(self) -> int:
		return len(self.rows)

	@property
	def n_cols(self) -> int:
		return len(self.cols)

	@classmethod
	def from_file(cls, fname, *, data_type=float):
		raw = numpy.loadtxt(fname, delimiter="\t", dtype=object)
		new = cls(
			data=raw[1:, 1:].astype(data_type),
			rows=raw[1:, 0],
			cols=raw[0, 1:])
		return new

	def transpose(self):
		return type(self)(
			data=self.data.T,
			rows=self.cols.copy(),
			cols=self.rows.copy(),
		)

	@property
	def T(self):
		return self.transpose()

	def filter_rows(self, *, by_index=None, by_tags=None) -> None:
		if (by_index is None) and (by_tags is None):
			raise ValueError("must specify by_index or by_tags")
		if (by_index is not None) and (by_tags is not None):
			raise ValueError("cannot specify both by_index and by_tags")
		if by_tags is not None:
			index = [self.rows_id_map[i] for i in by_tags]
			return self.filter_rows(by_index=index)
		if by_index is not None:
			return type(self)(
				data=self.data[by_index, :],
				rows=self.rows[by_index],
				cols=self.cols.copy(),
			)
		return

	def filter_cols(self, *, by_index=None, by_tags=None):
		if (by_index is None) and (by_tags is None):
			raise ValueError("must specify by_index or by_tags")
		if (by_index is not None) and (by_tags is not None):
			raise ValueError("cannot specify both by_index and by_tags")
		if by_tags is not None:
			index = [self.cols_id_map[i] for i in by_tags]
			return self.filter_cols(by_index=index)
		if by_index is not None:
			return type(self)(
				data=self.data[:, by_index],
				rows=self.rows.copy(),
				cols=self.cols[by_index],
			)
		return

	def grouped_row_mean(self, *, grouping: list):
		if len(grouping) != self.n_rows:
			raise ValueError("length of grouping must be identical to number of"
				" rows")

		grouping = numpy.asarray(grouping, dtype=object)
		group_unique = numpy.unique(grouping)
		grouped_mean = dict()
		for g in group_unique:
			mask = (grouping == g)
			grouped_mean[g] = self.data[mask, :].mean(axis=0, keepdims=True)
		data = numpy.vstack([grouped_mean[g] for g in group_unique])

		return type(self)(
			data=data,
			rows=group_unique,
			cols=self.cols.copy(),
		)

	def grouped_col_mean(self, *, grouping: list):
		if len(grouping) != self.n_cols:
			raise ValueError("length of grouping must be identical to number of"
				" rows")

		grouping = numpy.asarray(grouping, dtype=object)
		group_unique = numpy.unique(grouping)
		grouped_mean = dict()
		for g in group_unique:
			mask = (grouping == g)
			grouped_mean[g] = self.data[:, mask].mean(axis=1, keepdims=True)
		data = numpy.hstack([grouped_mean[g] for g in group_unique])

		return type(self)(
			data=data,
			rows=self.rows.copy(),
			cols=group_unique,
		)

	def to_file(self, fname, delimiter="\t"):
		with util.get_fp(fname, "w") as fp:
			# header
			fp.write(delimiter + delimiter.join(self.cols) + "\n")
			# each line
			for r, data in zip(self.rows, self.data):
				fp.write(delimiter.join([r] + [str(i) for i in data]) + "\n")
		return


class CountTable(Table):
	@classmethod
	def from_file(cls, fname):
		return super(cls, CountTable).from_file(fname, data_type=int)


class RelaAbundTable(Table):
	@classmethod
	def from_file(cls, fname):
		return super(cls, RelaAbundTable).from_file(fname, data_type=float)

	def get_alpha_diversity(self, metric="shannon"):
		data = skbio.diversity.alpha_diversity(metric, self.data).to_numpy()\
			.reshape(-1, 1)
		ret = type(self)(
			data=data,
			rows=self.rows,
			cols=["alpha_%s" % metric],
		)
		return ret
