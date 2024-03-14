#!/usr/bin/env python3

import argparse
import io
import sys

import Bio.SeqIO


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("--full-tax", "-t", type=str, required=True,
		metavar="tax",
		help="original complete tax file (required)")
	ap.add_argument("--truncated-fasta", "-f", type=str, required=True,
		metavar="fasta",
		help="truncated alignment fasta, with its headers being a subset of "
			"the input tax file (required)")
	ap.add_argument("--delimiter", "-d", type=str, default="\t",
		metavar="char",
		help="delimiter in input/output tax files [<tab>]")
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="tax",
		help="output recuited tax file [stdout]")

	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout

	return args


def load_fasta_ids_set(fname: str) -> set:
	ret = set()
	for seq in Bio.SeqIO.parse(fname, format="fasta"):
		ret.add(seq.id)
	return ret


def get_fp(f, *ka, factory=open, **kw) -> io.IOBase:
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument of get_fp() must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


def recruit_tax(ofile, tax_file: str, fasta_ids: set, delimiter="\t") -> None:
	with get_fp(tax_file, "r") as ifp:
		with get_fp(ofile, "w") as ofp:
			for i, line in enumerate(ifp):
				seq_id, *other_fields = line.rstrip().split(delimiter)
				if not other_fields:
					raise ValueError("apparent wrong format in file '%s', "
						"line '%u'" % (tax_file, i + 1))
				if seq_id in fasta_ids:
					ofp.write(line)
	return


def main():
	args = get_args()
	fasta_ids = load_fasta_ids_set(args.truncated_fasta)
	recruit_tax(args.output, args.full_tax, fasta_ids, delimiter=args.delimiter)
	return


if __name__ == "__main__":
	main()
