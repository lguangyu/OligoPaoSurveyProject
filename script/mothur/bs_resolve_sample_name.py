#!/usr/bin/env python3

import argparse
import io
import json
import subprocess
import sys
import time


def get_args():
	ap = argparse.ArgumentParser("resolves the downloaded fastq file name with "
		"biosample information, both related to a given project")
	ag = ap.add_mutually_exclusive_group(required = True)
	ag.add_argument("--id", "-i", type = int,
		metavar = "int",
		help = "project id in Illumina basespace (exlusive with --names/-n)")
	ag.add_argument("--name", "-n", type = str,
		metavar = "str",
		help = "project name in Illumina basespace (exlusive with --id/-i)")
	ap.add_argument("--basespace", type = str, default = "bs",
		metavar = "path",
		help = "path to basespace program [bs]")
	ap.add_argument("--config", "-c", type = str, default = None,
		metavar = "cfg",
		help = "basespace access token config file [default]")
	ap.add_argument("--output", "-o", type = str, default = "-",
		metavar = "json",
		help = "output file in json format [<stdout>]")
	ap.add_argument("--verbose", "-v", action = "store_true",
		help = "increase verbosity [off]")

	# parse and refine args	
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout
	return args


def get_fp(f, *ka, factory = open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first argument of get_fp must be str or io.IOBase, "
			"got '%s'" % type(f).__name__)
	return ret


def subprocess_run(cmd, *ka, capture_output = True, **kw):
	"""
	wrapped call to subprocess.run that captures stderr info when got non-zero
	return code
	"""
	ret = subprocess.run(cmd, capture_output = True)
	if ret.returncode:
		raise OSError("system call to '%s' had non-zero return status with "
			"message below:\n%s" % (str(cmd), ret.stderr.decode("utf-8")))
	return ret


def bs_get_project_content(proj, *, bs = "bs", config = None,
		verbose = None) -> list:
	if verbose:
		print("fetching project file content...", file = sys.stderr)
	cmd = [bs, "project", "content"] +\
			(["-c", config] if config else []) +\
			(["-i", str(proj)] if isinstance(proj, int) else ["-n", proj]) +\
			["-f", "json"]
	res = subprocess_run(cmd, capture_output = True)
	return json.loads(res.stdout)
	

def bs_get_file_linked_sample_map(proj_content: list, *, bs = "bs",
		config = None, sleep_sec = 1, verbose = None) -> dict:
	n_file = len(proj_content)
	ret = dict()
	for i, v in enumerate(proj_content):
		if verbose:
			print("fetching file sample name (%u out of %u)"\
				% (i + 1, n_file), file = sys.stderr)
		cmd = [bs, "file", "get"] +\
				(["-c", config] if config else []) +\
				["-i", v["Id"], "-f", "json"]
		res = subprocess_run(cmd, capture_output = True)
		ret[v["Path"]] = json.loads(res.stdout)["ParentDataSet"]["Name"]
		time.sleep(sleep_sec)
	return ret


def save_as_json(f, obj):
	with get_fp(f, "w") as fp:
		json.dump(obj, fp, indent = "\t", sort_keys = True)
	return


def main():
	args = get_args()
	proj_content = bs_get_project_content(args.id or args.name,
		bs = args.basespace,
		config = args.config,
		verbose = args.verbose)
	file_sample_map = bs_get_file_linked_sample_map(proj_content,
		bs = args.basespace,
		config = args.config,
		verbose = args.verbose)
	save_as_json(args.output, file_sample_map)
	return


if __name__ == "__main__":
	main()
