#!/usr/bin/env python3

import argparse
import json
import numpy
import skbio
import sklearn
import sys
# custom lib
import mpllayout
import pylib


def euclidean_distance(v1, v2) -> float:
	return sklearn.metrics.pairwise.euclidean_distances(
		v1.reshape(1, -1), v2.reshape(1, -1)
	)


def calc_wwtp_anosim(table: pylib.table.RelaAbundTable) -> None:
	grouping = [pylib.supp.SAMPLE_WWTP[i] for i in table.rows]
	dist_mat = skbio.DistanceMatrix.from_iterable(table.data,
		euclidean_distance, validate=False)
	res = skbio.stats.distance.anosim(dist_mat, grouping=grouping,
		permutations=10000)
	return res


def dump_json(f, obj) -> None:
	with pylib.util.get_fp(f, "w") as fp:
		json.dump(obj, fp)
	return


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		metavar="oligo-abund",
		help="oligo relative abundance matrix (default: stdin)")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="json",
		help="output JSON file of ANOSIM results (default: stdout)")
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	return args


def main() -> None:
	args = get_args()
	table = pylib.table.RelaAbundTable.from_file(args.input)
	anosim_res = calc_wwtp_anosim(table)
	dump_json(args.output, anosim_res.to_json())
	return


if __name__ == "__main__":
	main()
