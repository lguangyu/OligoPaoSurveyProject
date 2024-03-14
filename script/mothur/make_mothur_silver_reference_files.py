#!/usr/bin/env python3

import argparse
import multiprocessing
import os
import subprocess
import sys


MOTHUR_SILVA_REF_FILE_URL = {
	"v138": {
		"common": "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138.tgz",
		"seed": "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138.tgz",
	},
	"v132": {
		"common": "https://www.mothur.org/w/images/3/32/Silva.nr_v132.tgz",
		"seed": "https://www.mothur.org/w/images/7/71/Silva.seed_v132.tgz",
	},
	"v128": {
		"common": "https://www.mothur.org/w/images/b/b4/Silva.nr_v128.tgz",
		"seed": "https://www.mothur.org/w/images/a/a4/Silva.seed_v128.tgz",
	},
	"v123": {
		"common": "https://www.mothur.org/w/images/b/be/Silva.nr_v123.tgz",
		"seed": "https://www.mothur.org/w/images/1/15/Silva.seed_v123.tgz",
	},
}


class NonNegInt(int):
	def __new__(cls, *ka, **kw):
		new = super(NonNegInt, cls).__new__(cls, *ka, **kw)
		if new < 0:
			raise ValueError("NonNegInt must be positive, got %d" % new)
		return new


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("dbver", type = str, metavar = "silva-version",
		choices = list(MOTHUR_SILVA_REF_FILE_URL.keys()),
		help = "version of SILVA database; choices: %s"\
			% ((", ").join(sorted(MOTHUR_SILVA_REF_FILE_URL.keys()))))
	ap.add_argument("-O", "--output-dir", type = str, default = ".",
		metavar = "path",
		help = "output directory (default: .)")
	# number of available cpus
	n_cpus = multiprocessing.cpu_count()
	ap.add_argument("-t", "--num-threads", type = NonNegInt, default = n_cpus,
		metavar = "int",
		help = "number of threads used by mothur; 0 stands for using all "
			"(default: %d)" % n_cpus)
	ap.add_argument("-k", "--keep", action = "store_true",
		help = "do not delete downloaded archives (default: off)")
	ap.add_argument("--use-wget", type = str, default = "wget",
		metavar = "/path/to/wget", help = "use specified wget program")
	ap.add_argument("--use-tar", type = str, default = "tar",
		metavar = "/path/to/tar", help = "use specified tar program")
	ap.add_argument("--use-mothur", type = str, default = "mothur",
		metavar = "/path/to/mothur", help = "use specified mothur program")
	# parse and refine args
	args = ap.parse_args()
	if args.num_threads == 0:
		args.num_threads = n_cpus
	if args.num_threads > n_cpus:
		raise ValueError("-t/--num-threads cannot exceed physical number of "
			"cores (%d)" % n_cpus)
	return args


class MothurSilvaRefDB(object):
	def __init__(self, silva_ver, output_dir, *ka, **kw):
		if silva_ver not in MOTHUR_SILVA_REF_FILE_URL:
			raise ValueError("SILVA database version %s not supported" % silva_ver)
		self.silva_ver = silva_ver
		self.output_dir = output_dir
		self.url = MOTHUR_SILVA_REF_FILE_URL[silva_ver].copy()
		self._flist = list()
		return

	def log(self, msg: str, *, file = sys.stderr):
		print(msg, file = file)
		return self

	def call_external(self, cmd: list):
		retcode = subprocess.call(cmd)
		if retcode:
			raise SystemError("%s returned non-zero return code (%d)"\
				% (str(cmd), retcode))
		return self

	def download_files(self, wget = "wget"):
		os.makedirs(self.output_dir, exist_ok = True)
		for key, url in self.url.items():
			fname = os.path.basename(url) # also works for url
			output = os.path.join(self.output_dir, fname)
			self.log("downloading: %s -> %s" % (url, output))
			self.call_external([wget, "-O", output, url])
			# add file name to list
			self._flist.append(output)
		return self

	def extract_files(self, tar = "tar"):
		for fname in self._flist:
			self.log("extracting: %s" % fname)
			self.call_external([tar, "-zxC", self.output_dir, "-f", fname])
		return self

	def filter_align(self, threads, mothur = "mothur"):
		for c in ["nr", "seed"]:
			proto = "%s_%s" % (c, self.silva_ver)
			# file names
			ifname = os.path.join(self.output_dir, "silva.%s.align" % proto)
			pfname = os.path.join(self.output_dir, "silva.%s.pcr.align" % proto)
			ofname = ifname + ".v4"
			# mothur align
			self.log("mothur pcr: %s" % ifname)
			script = ("#pcr.seqs(fasta=%s, start=11894, end=25319, "\
				"keepdots=F, processors=%d)") % (ifname, threads)
			self.call_external([mothur, script])
			# rename output
			self.call_external(["mv", pfname, ofname])
		return self

	def finalize(self, remove = True):
		if remove:
			while self._flist:
				fname = self._flist.pop()
				self.call_external(["rm", "-f", fname])
		return


def main():
	args = get_args()
	db_maker = MothurSilvaRefDB(args.dbver, output_dir = args.output_dir)
	db_maker.download_files(args.use_wget)
	db_maker.extract_files(args.use_tar)
	db_maker.filter_align(args.num_threads, args.use_mothur)
	db_maker.finalize(remove = not args.keep)
	return


if __name__ == "__main__":
	main()
