#!/usr/bin/env python3

import networkx
import numpy
# custom lib
import pylib


NETWORK_CFG = [
	{
		"input": {
			"r": "misc_out/manuscript.fig4.oligo-oligo.corr.r.tsv",
			"p": "misc_out/manuscript.fig4.oligo-oligo.corr.p.tsv",
			"sig_flag": "misc_out/manuscript.fig4.oligo-oligo.corr.sig_flag.tsv",
		},
	},
]


def load_data_inplace(cfg: dict) -> dict:
	res = dict()
	for k, v in cfg["input"].items():
		res[k] = pylib.table.RelaAbundTable.from_file(v)

	cfg["table"] = res
	return res


def network_analysis(cfg: dict):
	table_dict = load_data_inplace(cfg)

	r = table_dict["r"].data
	p = table_dict["p"].data
	s = table_dict["sig_flag"].data

	adj_mat = (numpy.abs(r) > 0.5) & (p < 0.01)
	adj_mat[numpy.isnan(adj_mat)] = False

	G = networkx.from_numpy_array(adj_mat)
	comm_maxq = networkx.algorithms.community.modularity_max\
		.greedy_modularity_communities(G)

	for c in comm_maxq:
		print((",").join([table_dict["r"].rows[i] for i in c]))

	q = (networkx.algorithms.community.modularity(G, comm_maxq))
	print(q)
	return


def main():
	for cfg in NETWORK_CFG:
		network_analysis(cfg)
	return


if __name__ == "__main__":
	main()
