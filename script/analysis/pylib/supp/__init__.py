#!/usr/bin/env python3

import os

SUPP_DIR = os.path.dirname(__file__)


class UserDict(dict):
	pass


def read_map_file(fname, delimiter="\t") -> UserDict:
	ret = UserDict()
	ret.def_key_order = list()
	with open(fname, "r") as fp:
		for line in fp:
			line = line.rstrip("\r\n")
			if (not line) or (line.startswith("#")):
				continue
			# this also checks a 2-column format, other wise raises error
			k, v = line.rstrip("\r\n").split(delimiter)
			ret[k] = v
			ret.def_key_order.append(k)
	return ret


SAMPLE_WWTP = read_map_file(os.path.join(SUPP_DIR, "sample_wwtp.map"))
WWTP_CATE = read_map_file(os.path.join(SUPP_DIR, "wwtp_cate.map"))
WWTP_DISPLAY = read_map_file(os.path.join(SUPP_DIR, "wwtp_display.map"))
WWTP_DISPLAY_ABBR = read_map_file(
	os.path.join(SUPP_DIR, "wwtp_display_abbr.map"))
CATE_DISPLAY = read_map_file(os.path.join(SUPP_DIR, "cate_display.map"))
CATE_COLOR = read_map_file(os.path.join(SUPP_DIR, "cate_color.map"))
