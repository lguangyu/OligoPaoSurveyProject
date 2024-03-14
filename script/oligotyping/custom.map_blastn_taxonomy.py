#!/usr/bin/env python3

import argparse
import io
import sys


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("--hit-accs", "-i", type=str, required=True,
		metavar="txt",
		help="single-column hit accession list (required)")
	ap.add_argument("--accs-tax-map", "-t", type=str, required=True,
		metavar="table",
		help="2-column accession to taxonomy map list, with the 1st column as "
			"accession and the 2nd column as taxonomy (required)")
	ap.add_argument("--delimiter", "-d", type=str, default="\t",
		metavar="char",
		help="delimiter in tabular input and output [tab]")
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="table",
		help="output table [stdout]")

	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout

	return args


def get_fp(f, *ka, factory=open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TabError("first argument of get_fp() must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


def cleanup_midas_taxon_name(s: str, delimiter=";") -> str:
	fields = s.rstrip(delimiter).split(delimiter)
	if fields[-1].startswith("midas_"):
		ret = fields[-2] + "_" + fields[-1]
	else:
		ret = fields[-1]
	return ret


def load_accs_tax_map(f, delimiter="\t") -> dict:
	ret = dict()
	with get_fp(f, "r") as fp:
		for line in fp:
			accs, tax = line.rstrip().split(delimiter)
			# clean up tax, we only want the last section wichi is species name
			ret[accs] = cleanup_midas_taxon_name(tax)
	return ret


def map_taxonomy_to_hit_accs(ifile, ofile, *, accs_tax_map: dict,
		delimiter="\t") -> None:
	with get_fp(ifile, "r") as ifp:
		hit_accs = ifp.read().splitlines()

	with get_fp(ofile, "w") as ofp:
		for accs in hit_accs:
			print(delimiter.join([accs, accs_tax_map[accs]]), file=ofp)

	return


def main():
	args = get_args()
	accs_tax_map = load_accs_tax_map(args.accs_tax_map,
		delimiter=args.delimiter
	)
	map_taxonomy_to_hit_accs(args.hit_accs, args.output,
		accs_tax_map=accs_tax_map,
		delimiter=args.delimiter
	)
	return


if __name__ == "__main__":
	main()
