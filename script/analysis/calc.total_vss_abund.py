#!/usr/bin/env python3

import argparse
import sys
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		metavar="oligo-abund",
		help="oligo abundance matrix before compositional transformation "
			"(default: stdin)")
	ap.add_argument("--rela-abund", "-a", type=str, required=True,
		metavar="tsv",
		help="total relative abundance of the oligotyping genus in each VSS "
			"sample (required)"
	)
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="tsv",
		help="output oligotype abundance table in relative to VSS "
			"(default: stdout)")

	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


def main():
	args = get_args()
	oligo_abund = pylib.table.RelaAbundTable.from_file(args.input)
	vss_abund = pylib.table.RelaAbundTable.from_file(args.rela_abund)
	oligo_abund.data *= vss_abund.data
	oligo_abund.to_file(args.output)
	return


if __name__ == "__main__":
	main()
