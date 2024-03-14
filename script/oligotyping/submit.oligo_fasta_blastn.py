#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
import typing
import warnings

import Bio.SeqIO # for reading fasta


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("oligo_output", type=str,
		help="oligotyping output directory")
	ap.add_argument("--max-n-jobs", "-j", type=int, default=1,
		metavar="int",
		help="number of maximum SLURM jobs to submit [1]")
	ap.add_argument("--output-dir", "-O", type=str, default="blastn",
		metavar="dir",
		help="output directory [blastn]")
	ap.add_argument("--log-dir", type=str, default=".log",
		metavar="dir",
		help="log directory [.log]")
	ap.add_argument("--dry-run", "-N", action="store_true",
		help="do not submit any jobs or make any changes")

	# parse and refine args
	args = ap.parse_args()
	if args.max_n_jobs < 1:
		warnings.warn("--max-n-jobs got an invalid value: '%d' and is reset to "
			"1" % args.max_n_jobs)
		args.max_n_jobs = 1 # fix offending values to the default

	return args


class OligoRepUniqStats(object):
	def __init__(self, files: list, *ka, **kw):
		super().__init__(*ka, **kw)
		self.files = files
		self._stat_file_num_seqs()
		return

	@classmethod
	def scan_oligo_output(cls, path: str, scan_suffix="_unique"):
		if not os.path.isdir(path):
			raise IOError("'%s' is not a directory" % path)
		rep_dir = os.path.join(path, "OLIGO-REPRESENTATIVES")
		if not os.path.isdir(rep_dir):
			raise IOError("dir '%s' is missing, may be incomplete oligotyping "
				"output" % path)
		# scan dir for target files
		files = list()
		for i in os.scandir(rep_dir):
			if i.name.endswith(scan_suffix):
				files.append(i.path)

		new = cls(files=files)
		return new

	@property
	def num_seqs(self) -> tuple:
		return self._num_seqs

	def _stat_file_num_seqs(self):
		num_seqs = list()
		for i in self.files:
			try:
				# use a formal fasta parser can also help check the file format
				n = len((list(Bio.SeqIO.parse(i, format="fasta"))))
			except:
				print("fail to parse file: %s" % i, file=sys.stderr)
				sys.exit(-1)
			num_seqs.append(n)
		self._num_seqs = tuple(num_seqs)
		return


class OligoRepBlastJobSubmit(object):
	def __init__(self, *ka, oligo_output: str, output_dir: str, log_dir: str,
			max_n_jobs: int = 1, **kw):
		super().__init__(*ka, **kw)
		self.oligo_output = oligo_output
		self.output_dir = output_dir
		self.log_dir = log_dir
		self.max_n_jobs = max_n_jobs
		self.fasta_stats = OligoRepUniqStats.scan_oligo_output(oligo_output)
		return

	def submit_jobs(self, dry_run=False) -> int:
		# create output dirs
		if not dry_run:
			os.makedirs(self.output_dir, exist_ok=True)
			os.makedirs(self.log_dir, exist_ok=True)

		fasta_lists = self.split_job_fasta_lists()
		# submit individual worker jobs
		worker_input_files = list()
		for i, flist in enumerate(fasta_lists):
			# save flist for workers to read
			flist_file = os.path.join(self.output_dir, "split_%03u.tmp" % i)
			if not dry_run:
				with open(flist_file, "w") as fp:
					for f in flist:
						print(f, file=fp)
			worker_input_files.append(flist_file)

		# submit worker jobs 
		return self._submit_slurm_workers(worker_input_files, dry_run=dry_run)

	def _submit_slurm_workers(self, worker_input_files: list, *, dry_run=False
			) -> int:
		jobids = list()

		for f in worker_input_files:
			job_name = "oligo_blastn." + os.path.basename(f)
			log_file = os.path.join(self.log_dir, job_name + ".log")
			cmd = ["sbatch", "-J", job_name, "-o", log_file,
				"script/worker.oligo_fasta_blastn.sh", f]
			if not dry_run:
				sp = subprocess.run(cmd, stdout=subprocess.PIPE)
				print(sp.stdout.decode("utf-8"), file=sys.stdout, end="")
				if sp.returncode:
					print(str(cmd) + " exited with non-zero return code",
						file=sys.stderr)
					self._clean_up_submitted_slurm_jobs(jobids)
					print("aborting", file=sys.stderr)
					return -1
				else:
					m = re.search(r"(\d+)$", sp.stdout.decode("utf-8"))
					if m:
						jobids.append(m.group(1))
			else:
				print("submitting: " + str(cmd), file=sys.stderr)

		return 0

	@staticmethod
	def	_clean_up_submitted_slurm_jobs(jobids: list) -> None:
		if not jobids:
			return
		print("clean up submitted jobs: " + ((" ").join(jobids)),
			file=sys.stderr)
		cmd = ["scancel"] + jobids
		subprocess.run(cmd)
		return

	def split_job_fasta_lists(self) -> list:
		ret = list()
		for s in self._split_solve(self.fasta_stats.num_seqs, self.max_n_jobs):
			ret.append([self.fasta_stats.files[i] for i in s])
		return ret

	@staticmethod
	def _split_attempt(arr: list, s: int, k: int) -> typing.Optional[list]:
		"""
		try to split the <ilist> (list of integers) into <k> number of sublists,
		s.t. minimizing the max sum among all sublists; intending to split child
		jobs as even as possible;

		arr: an array of integers, representing the size of each atomic job
		s: the desired max of sublist sums
		k: the number of sublists

		return: a list of lists of indices of elements in ilist; can be None if
			the desired split given by parameter s and k is not feasible
		"""
		sub_sums = list()
		sub_idxs = list()

		for i, v in enumerate(arr):
			#if v > s:
			#	# not feasible since a single element is larger than s
			#	return None

			for j in range(len(sub_sums)):
				if sub_sums[j] + v <= s:
					# the sublist can accept the current element
					sub_sums[j] += v
					sub_idxs[j].append(i)
					break
			else:
				# if reach here, need to add a new sublist
				sub_sums.append(v)
				sub_idxs.append([i])

		# check if the number of sublists is within k
		return sub_idxs if (len(sub_idxs) <= k) else None

	@staticmethod
	def _split_solve(arr: list, k: int) -> list:
		if not arr:
			raise ValueError("arr cannot be empty")
		if k < 1:
			raise ValueError("k must be positive, got '%d'" % k)
		if k == 1:
			# no need to go further if split into 1 split...
			return [list(range(len(arr)))]

		# lower bound of searching range
		ran_l = min(arr)
		# higher bound of searching range
		ran_h = sum(arr)

		# using a binary-search-like approach
		solution = None
		while ran_l <= ran_h:
			s = (ran_l + ran_h) // 2
			attempt = OligoRepBlastJobSubmit._split_attempt(arr, s, k)
			if attempt is not None:
				solution = attempt
				ran_h = s - 1
			else:
				ran_l = s + 1

		return solution


def main():
	args = get_args()
	o = OligoRepBlastJobSubmit(
		oligo_output=args.oligo_output,
		output_dir=args.output_dir,
		log_dir=args.log_dir,
		max_n_jobs=args.max_n_jobs,
	)
	o.submit_jobs(dry_run=args.dry_run)
	return


if __name__ == "__main__":
	main()
