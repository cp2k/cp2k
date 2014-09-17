#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from os import path
from os.path import basename
from glob import glob
from itertools import product, chain
from optparse import OptionParser
import re
from pprint import pprint

re_mnk    = re.compile("tune_(\d+)x(\d+)x(\d+)_")
re_winner = re.compile("\nWINNER: \d+ (.+)\n")
re_gflops = re.compile("# ([0-9.]+) GFlop/?s")
re_errors = re.compile("Number of errors: (\d+)\n")

#===============================================================================
def main():
	winners = dict()

	for d in glob("tune_*"):
		if(not path.isdir(d)):
			continue
					            
		for exe_fn in glob(d+"/tune_*main.cu"):
			mnk = tuple([int(i) for i in re_mnk.search(exe_fn).groups()])
			log_fn = exe_fn.replace("_main.cu", ".log")
			if(not path.exists(log_fn)):
				winners[mnk] = "log missing: "+log_fn
				continue
			
			process_log(log_fn, mnk, winners)

	for w in winners.values():
		print w
	##pprint(winners)

#===============================================================================
def process_log(log_fn, mnk, winners):
	print("Reading: "+log_fn)

	f = open(log_fn)
	#f.seek(-1000, os.SEEK_END)
	content = f.read()
	f.close()
	
	m = re_errors.search(content)
	if(not m):
		winners[mnk] = "log incomplete: "+log_fn
		return

	n_errors = int(m.group(1))
	if(n_errors != 0):
		winners[mnk] = "errors: "+log_fn
		return
	
	old_gflops = 0.0
	if(winners.has_key(mnk)):
		m = re_gflops.search(winners[mnk])
		if(not m):
			return
		old_gflops = float(m.group(1))
	
	new_winner = re_winner.search(content).group(1)
	new_gflops = float(re_gflops.search(new_winner).group(1))

	if(new_gflops > old_gflops):
		winners[mnk] = new_winner


#===============================================================================

main()

#EOF
