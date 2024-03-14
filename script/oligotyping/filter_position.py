#!/usr/bin/env python3

import argparse
import io
import sys


def get_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		help="2-column position-entropy table [stdin]")
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="txt",
		help="filtered numerical positions for oligotyping, comma-separated "
			"[stdout]")
	ap.add_argument("--threshold", "-t", type=float, default=0.2,
		metavar="float",
		help="filter out all positions below this threshold [0.2]")

	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout

	return args


def get_fp(f, *ka, factory=open, **kw) -> io.IOBase:
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = open(f, *ka, **kw)
	else:
		raise TypeError("first argument must be str or io.IOBase, got '%s'" \
			% type(f).__name__)
	return ret


def load_postion_entropy(f, delimiter="\t") -> list:
	# implement without using numpy
	ret = list()
	with get_fp(f, "r") as fp:
		for ln, line in enumerate(fp):
			fields = line.rstrip().split(delimiter)
			if len(fields) != 2:
				raise ValueError("incorrect number of fields at line %u: "
					"expect 2" % (ln + 1))
			pos = int(fields[0])
			entropy = float(fields[1])
			ret.append((pos, entropy))
	return ret


def filter_positions(pos_entropy: list, threshold: float) -> list:
	## sort by entropy in descending order
	#sorted_pos_entropy = sorted(pos_entropy, key=lambda x: x[1], reverse=True)
	# filter
	ret = [i[0] for i in filter(lambda x: x[1] >= threshold, pos_entropy)]
	return ret


def save_filtered_positions(f, filtered_pos) -> None:
	with get_fp(f, "w") as fp:
		fp.write((",").join([str(i) for i in filtered_pos]))
	return


def main():
	args = get_args()
	pos_entropy = load_postion_entropy(args.input)
	filtered_pos = filter_positions(pos_entropy, threshold=args.threshold)
	save_filtered_positions(args.output, filtered_pos)
	return


if __name__ == "__main__":
	main()

