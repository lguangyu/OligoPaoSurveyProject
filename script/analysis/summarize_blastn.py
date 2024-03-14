#!/usr/bin/env python3

import argparse
import collections
import dataclasses
import json
import math
import re
import sys
import typing
import os
# custom lib
import pylib


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("blastn-dir", type=str,
		metavar="dir",
		help="the directory for per-oligotype blastn results")
	ap.add_argument("-x", "--extension", type=str, default="blastn",
		metavar="str",
		help="only files ending with this extension will be discovered and "
			"processed [.blastn]")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="json",
		help="the output json file [stdout]")
	ap.add_argument("-H", "--human-readable", action="store_true",
		help="make the output json in human-readable format [off]")
	ap.add_argument("-c", "--cumu-frac-thres", type=pylib.util.Fraction,
		default=pylib.util.Fraction(0.95),
		metavar="float",
		help="minor seqs encoutered after reaching this cumulative seq count "
			"threshold will be discarded, 1.0 means no discarding; this option "
			"is functionally in addition to --min-frac-thres/-m [0.95]"
		)
	ap.add_argument("-m", "--min-frac-thres", type=pylib.util.Fraction,
		default=pylib.util.Fraction(0.05),
		metavar="float",
		help="seqs must be higher than this fraction to be selected, 0 means "
			"don't apply this filter; this option is functionally in addition "
			"to --cumu-frac-thres/-c [0.0]"
		)
	ap.add_argument("-t", "--taxonomy", type=str,
		metavar="txt",
		help="2-column taxonomy file with the first column as the seq ids and "
			"the second column as taxonomy labels; taxonomy labels should "
			"follow the format kingdom;phylum;class;order;family;genus;species "
			"without omitting any fields; this is an optional argument, when "
			"provided the seq ids will be translated into taxonomy labels"
		)
	ap.add_argument("--save-seq-ids", type=str,
		metavar="json",
		help="if set, save the filtered seq ids into the specified json file "
			"[no]")

	# parse and refine args
	args = ap.parse_args()
	# force args.extension starting with separator, typically dot (.), if user
	# input does not include it
	if not args.extension.startswith(os.extsep):
		args.extension = os.extsep + args.extension
	if args.output == "-":
		args.output = sys.stdout

	return args


def load_tax_map(fname, delimiter="\t", taxsep=";") -> dict:
	ret = dict()
	with pylib.util.get_fp(fname, "r") as fp:
		for line in fp.read().splitlines():
			seq_id, tax = line.split(delimiter)
			ret[seq_id] = pylib.util.Taxonomy.parse(tax)
	return ret


class TaxonomyMap(dict):
	@ classmethod
	def load_tax_file(cls, f, *, delimiter="\t", taxsep=";"):
		new = cls()
		with pylib.util.get_fp(f, "r") as fp:
			for line in fp:
				seq_id, tax = line.rstrip("\r\n").split(delimiter)
				k, p, c, o, f, g, s, *_ = tax.split(taxsep)
				new[seq_id] = ("_").join([g, s])
		return new


@ dataclasses.dataclass
class BlastTable(object):

	@ dataclasses.dataclass
	class Record(object):
		seq: str
		freq: int
		hit_seq: str
		identity: float
		evalue: float
		bitscore: float

		@ classmethod
		def parse(cls, s: str, *, delimiter="\t"):
			splits = s.rstrip().split(delimiter)
			seq = splits[0]
			freq = int(seq.split("|freq:")[1])
			new = cls(
				seq=seq,
				freq=freq,
				hit_seq=splits[1],
				identity=float(splits[2]),
				evalue=float(splits[3]),
				bitscore=float(splits[4]),
			)
			return new

	hit_dict: collections.defaultdict = dataclasses.field(
		default_factory=lambda: collections.defaultdict(list))
	freq_counter: collections.Counter = dataclasses.field(
		default_factory=lambda: collections.Counter()
	)
	_oligotype: typing.Optional[str] = dataclasses.field(default=None)

	@ property
	def empty(self) -> bool:
		return len(self.hit_dict) == 0

	@ property
	def oligotype(self) -> str:
		return self._oligotype

	@ oligotype.setter
	def oligotype(self, s: str):
		if self._oligotype is None:
			self._oligotype = s
		elif self._oligotype != s:
			raise ValueError("cannot set conflicting oligotypes: already set: "
				"%s; offender: %s" % (self._oligotype, s))
		return

	def add_record(self, rec: Record):
		self.hit_dict[rec.seq].append(rec)
		self.freq_counter[rec.seq] = rec.freq
		self.oligotype = rec.seq.split("_")[0]
		return

	@ classmethod
	def parse_file(cls, f, delimiter="\t"):
		new = cls()
		with pylib.util.get_fp(f, "r") as fp:
			for line in fp:
				rec = cls.Record.parse(line)
				new.add_record(rec)
		return new


def discover_files_in_ext(dir_path, ext: str, *, recursive=False):
	for i in os.scandir(dir_path):
		if i.is_dir() and recursive:
			yield from discover_files_in_ext(dir_path, ext)
		elif i.is_file() and os.path.splitext(i.path)[1] == ext:
			yield i.path
	return


def load_oligo_blastn_dir(dir_path, ext) -> dict:
	ret = dict()
	for f in discover_files_in_ext(dir_path, ext=ext):
		table = BlastTable.parse_file(f)
		if table.empty:
			continue

		ret[table.oligotype] = table

	return ret


def list_seq_best_indentity_hits(rec_list: list,
		tax_map: typing.Optional[dict] = None,
	) -> list:
	best_identity = max(map(lambda rec: rec.identity, rec_list))
	ret = [rec.hit_seq for rec in rec_list if rec.identity == best_identity]
	if tax_map:
		ret = [tax_map[i].s for i in ret]
	return ret


def filter_oligo_blastn_result(table: BlastTable, *,
		cumu_frac_thres: float = 1.0,
		min_frac_thres: float = 0.0,
		tax_map: typing.Optional[dict] = None,
	) -> (list, list):
	total_freq = table.freq_counter.total()
	cumu_count_thres = math.ceil(total_freq * cumu_frac_thres)
	min_count_thres = math.ceil(total_freq * min_frac_thres)
	curr_count = 0
	hits, major_seqs = list(), list()

	for seq, freq in table.freq_counter.most_common():
		if freq < min_count_thres:
			continue
		hits.extend(list_seq_best_indentity_hits(table.hit_dict[seq],
			tax_map=tax_map))
		major_seqs.append(seq)
		# stop checking more seqs when reaching count threshold
		curr_count += freq
		if curr_count >= cumu_count_thres:
			break

	return hits, major_seqs


def filter_oligo_blastn_results(raw_res: dict, *,
		cumu_frac_thres: float = 1.0,
		min_frac_thres: float = 0.0,
		tax_map: typing.Optional[dict] = None,
	) -> (dict, dict):
	hit_dict, seq_dict = dict(), dict()
	for k, v in raw_res.items():
		hits, seqs = filter_oligo_blastn_result(v,
			cumu_frac_thres=cumu_frac_thres,
			min_frac_thres=min_frac_thres,
			tax_map=tax_map
		)
		hit_dict[k] = hits
		seq_dict[k] = seqs
	return hit_dict, seq_dict


def save_json(f, obj: dict, *, human_readable=False) -> None:
	if f is None:
		return
	format_kw = dict(sort_keys=True, indent="\t") if human_readable else dict()
	with pylib.util.get_fp(f, "w") as fp:
		json.dump(obj, fp, **format_kw)
	return


def main():
	args = get_args()
	tax_map = load_tax_map(args.taxonomy) if args.taxonomy else None
	raw_hit = load_oligo_blastn_dir(getattr(args, "blastn-dir"), args.extension)
	flt_hit, flt_seq = filter_oligo_blastn_results(raw_hit,
		cumu_frac_thres=args.cumu_frac_thres,
		min_frac_thres=args.min_frac_thres,
		tax_map=tax_map,
	)
	save_json(args.output, flt_hit, human_readable=args.human_readable)
	save_json(args.save_seq_ids, flt_seq, human_readable=args.human_readable)
	return


if __name__ == "__main__":
	main()
