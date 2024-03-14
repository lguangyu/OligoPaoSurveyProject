#!/usr/bin/env python3

import argparse
import json
import numpy
import os
import scipy.stats
import sys
import typing
# custom lib
import pylib


class RelaAbundTableTransform(pylib.table.RelaAbundTable):
	# compositional data transformations
	@staticmethod
	def _gmean_1d(vector):
		if (vector <= 0).any():
			raise ValueError("data in vector must be positive")
		s = vector.sum()
		if numpy.isclose(s, 1):
			gmean = scipy.stats.gmean(vector)
		elif s < 1:
			n = len(vector)
			gmean = scipy.stats.gmean(numpy.append(vector, 1.0 - s))
			# checks numpy.append above not to edit vector
			assert n == len(vector)
		else:
			raise ValueError("compositional data cannot be with a sum > 1")
		return gmean

	@staticmethod
	def _clr_1d(vector):
		return numpy.log(vector / RelaAbundTableTransform._gmean_1d(vector))

	def _calc_trans_meth_clr(self) -> numpy.ndarray:
		# centered log ratio transformation
		return pylib.table.RelaAbundTable(
			data=numpy.array([self._clr_1d(i) for i in self.data]),
			rows=self.rows,
			cols=self.cols)

	@staticmethod
	def _log_clr_1d(vector):
		ret = RelaAbundTableTransform._clr_1d(vector)
		lmask, umask = (ret <= 0), (ret > 0)
		ret[lmask] = -(numpy.log(1 - ret[lmask]) ** 2)
		ret[umask] = ret[umask] ** 2
		return ret

	def _calc_trans_meth_log_clr(self) -> numpy.ndarray:
		# log centered log ratio
		return pylib.table.RelaAbundTable(
			data=numpy.array([self._log_clr_1d(i) for i in self.data]),
			rows=self.rows,
			cols=self.cols)

	def _calc_trans_meth_none(self):
		return self

	COMP_TRANS_METH = dict()
	COMP_TRANS_METH["none"] = "_calc_trans_meth_none"
	COMP_TRANS_METH["clr"] = "_calc_trans_meth_clr"
	COMP_TRANS_METH["log-clr"] = "_calc_trans_meth_log_clr"

	def transform(self, meth):
		return getattr(self, self.COMP_TRANS_METH[meth])()

	def mark_oligos_above_abund(self, *, inplace=False,
			abund_thres: float = 0,
			list_output: typing.Optional[str] = None,
			no_renumbering: bool = False,
			rename_prefix: typing.Optional[str] = None,
			rename_output: typing.Optional[str] = None,
		):
		if (abund_thres > 0) and (list_output is None):
			raise ValueError("must provide an output file name to list abundant"
				" oligos")
		# rename abundant oligots
		mask = self.data.max(axis=0) >= abund_thres
		rename_map = dict()
		new_cols = list()
		abund_oligo = list()

		# this accouts for empty string, treat as None
		rename_pattern = ("%s_%%u" % rename_prefix) if rename_prefix else "%u"

		for i, is_abund in enumerate(mask):
			if is_abund:
				oligo_id = i if no_renumbering else len(abund_oligo)
				new_name = rename_pattern % oligo_id
				rename_map[new_name] = self.cols[i]
				new_cols.append(new_name)
				abund_oligo.append(new_name)
			else:
				new_cols.append(self.cols[i])

		# output abund oligo list
		if list_output:
			with pylib.util.get_fp(list_output, "w") as fp:
				print(os.linesep.join(abund_oligo), file=fp)

		# outout oligo rename map
		rename_json_str = json.dumps(rename_map)
		if rename_output is None:
			print(rename_json_str, file=sys.stderr)
		else:
			with pylib.util.get_fp(rename_output, "w") as fp:
				print(rename_json_str, file=fp)

		# return a new table
		if inplace:
			ret = self.set_all_data(
				data=self.data,
				rows=self.rows,
				cols=new_cols,
			)
		else:
			ret = type(self)(
				data=self.data.copy(),
				rows=self.rows.copy(),
				cols=new_cols,
			)
		return ret


def calc_oligo_rela_abund(oligo_count, *,
		pseudo_count: pylib.util.NonNegFloat = None,
	) -> RelaAbundTableTransform:
	# pseudo_count must be positive
	if (pseudo_count is not None) and (pseudo_count < 0):
		raise ValueError("pseudo_count must be non-negative")
	rows = oligo_count.rows.copy()
	cols = oligo_count.cols.copy()
	# the default pseudo count makes the row sum 1% higher
	if pseudo_count is None:
		pseudo_count = 0.01 * oligo_count.data.sum(axis=1, keepdims=True)\
			/ oligo_count.n_cols
	data = oligo_count.data + pseudo_count
	data /= data.sum(axis=1, keepdims=True)
	return RelaAbundTableTransform(data, rows, cols)


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		metavar="oligo-count",
		help="oligo count matrix from oligotying (default: stdin)")
	ap.add_argument("--pseudo-count", "-c", type=pylib.util.NonNegFloat,
		default=None, metavar="float",
		help="pseudo-count used in CLR and logCLR transformation, can be any "
			"positive real number; leave empty to use the auto mode, which "
			"will add pseudo count sums to 1 percent of raw count sum "
			"(default: auto)")
	ap.add_argument("--transform", "-t", type=str, default="clr",
		choices=sorted(RelaAbundTableTransform.COMP_TRANS_METH.keys()),
		help="transformation of compositional data (default: clr)")

	ag = ap.add_argument_group("output")
	ag.add_argument("--abund-output", "-a", type=str, default=None,
		metavar="tsv",
		help="output abundance table before compositional transformation "
			"(default: no)")
	ag.add_argument("--output", "-o", type=str, default="-",
		metavar="tsv",
		help="output transformed compositional data table (default: stdout)")
	ag.add_argument("--output-abund-table", type=str,
		metavar="tsv",
		help="output abundance table (default: omitted)")
	ag.add_argument("--output-alpha-diversity", type=str,
		metavar="tsv",
		help="output alpha diversity (default: omitted)")

	ag = ap.add_argument_group("misc")
	ag.add_argument("--abund-oligo-thres", type=pylib.util.NonNegFloat,
		metavar="float",
		help="filter oligos with abundance above this threshold in at least one"
			"sample; this filtering will not effect the --output; therefore, "
			"this option as no effect if --abund-oligo-list is not set "
			"(default: omitted)")
	ag.add_argument("--abund-oligo-output", type=str, default=None,
		metavar="txt",
		help="list oligo names with abundance above --abund-oligo-thres in at "
			"least one sample; has no effect if --abund-oligo-thres is not set "
			"(default: omitted)")
	ag.add_argument("--oligo-rename-prefix", type=str,
		metavar="str",
		help="if set, rename oligos with this prefix, only applies to those "
			"passed the abundance filtering (default: omitted)")
	ag.add_argument("--no-renumbering", action="store_true",
		help="do not renumber oligos that pass the abundance filtering, but "
			"set numbering using the indices of the original oligo table "
			"(default: no)")
	ag.add_argument("--oligo-rename-result", type=str,
		metavar="json",
		help="if set, save the new/old name mapping in json; otherwise, show in"
			" stderr (default: <stderr>")

	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


def main():
	args = get_args()
	oligo_count = pylib.table.CountTable.from_file(args.input)

	# calculate abundance table from count matrix
	oligo_abund = calc_oligo_rela_abund(oligo_count,
		pseudo_count=args.pseudo_count)

	# list abundant oligo and rename those oligos if necessary
	if args.abund_oligo_thres:
		oligo_abund.mark_oligos_above_abund(
			inplace=True,
			abund_thres=args.abund_oligo_thres,
			list_output=args.abund_oligo_output,
			no_renumbering=args.no_renumbering,
			rename_prefix=args.oligo_rename_prefix,
			rename_output=args.oligo_rename_result,
		)

	if args.abund_output:
		oligo_abund.to_file(args.abund_output)

	# output abundance table and diversity
	if args.output_abund_table:
		oligo_abund.to_file(args.output_abund_table)
	if args.output_alpha_diversity:
		oligo_abund.get_alpha_diversity().to_file(args.output_alpha_diversity)

	# transform into compositional data
	comp_data = oligo_abund.transform(args.transform)
	# output
	comp_data.to_file(args.output)
	return


if __name__ == "__main__":
	main()
