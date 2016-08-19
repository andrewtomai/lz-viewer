#!/usr/bin/env python2

'''
This script takes two arguments:
- an input filename (the output of epacts with a single phenotype)
- an output filename (please end it with `.json`)

It creates a json file which can be used to render a QQ plot.
'''

from __future__ import print_function, division, absolute_import

import os.path
import sys

import re
import json
import math
import collections
from lz_assoc_viewer import Data_reader
import scipy.stats

data_dir = '/var/pheweb_data/'
NEGLOG10_PVAL_BIN_SIZE = 0.05 # Use 0.05, 0.1, 0.15, etc
NEGLOG10_PVAL_BIN_DIGITS = 2 # Then keep 2 digits after the decimal
NUM_BINS = 1000

global NUM_MAF_RANGES
NUM_MAF_RANGES = 4 # Split MAF into 4 equally-sized ranges


def round_sig(x, digits):
	return 0 if x==0 else round(x, digits-1-int(math.floor(math.log10(abs(x)))))
assert round_sig(0.00123, 2) == 0.0012
assert round_sig(1.59e-10, 2) == 1.6e-10


def approx_equal(a, b, tolerance=1e-4):
	return abs(a-b) <= max(abs(a), abs(b)) * tolerance
assert approx_equal(42, 42.0000001)
assert not approx_equal(42, 42.01)


def gc_value(median_pval):
     # This should be equivalent to this R: `qchisq(p, df=1, lower.tail=F) / qchisq(.5, df=1, lower.tail=F)`
	return scipy.stats.chi2.ppf(1 - median_pval, 1) / scipy.stats.chi2.ppf(0.5, 1)
assert approx_equal(gc_value(0.49), 1.047457) # I computed these using that R code.
assert approx_equal(gc_value(0.5), 1)
assert approx_equal(gc_value(0.50001), 0.9999533)
assert approx_equal(gc_value(0.6123), 0.5645607)


Variant = collections.namedtuple('Variant', ['neglog10_pval', 'maf'])
def parse_variant_line(file_reader):
	
	if file_reader.get_pval() == 'NA':
		return None
	elif file_reader.get_maf() is None:
		maf = float(1)
		pval = float(file_reader.get_pval())
		global NUM_MAF_RANGES
		NUM_MAF_RANGES = 1
		return Variant(-math.log10(pval), maf)
	else:
		maf = float(file_reader.get_maf())
		pval = float(file_reader.get_pval())

		return Variant(-math.log10(pval), maf)


def rounded(x):
	return round(x // NEGLOG10_PVAL_BIN_SIZE * NEGLOG10_PVAL_BIN_SIZE, NEGLOG10_PVAL_BIN_DIGITS)



##FOR EPACTS FILES
def make_qq_stratified(file_reader):
	variants = []
	
	#skip the header info
	file_reader.skip_header()

	while True: 
		#get a new line
		file_reader.get_line()

		#check if we are at the end
		if file_reader.is_end():
			break
		variant = parse_variant_line(file_reader)
		if variant is not None:
			variants.append(variant)
	variants = sorted(variants, key=lambda v: v.maf)
	

	##QQ
	
	num_variants_in_biggest_maf_range = int(math.ceil(len(variants) / NUM_MAF_RANGES))
	max_exp_neglog10_pval = -math.log10(0.5 / num_variants_in_biggest_maf_range) #expected
	max_obs_neglog10_pval = max(v.neglog10_pval for v in variants) #observed
	# TODO: should max_obs_neglog10_pval be at most 9?  That'd break our assertions.
	# print(max_exp_neglog10_pval, max_obs_neglog10_pval)

	qqs = [dict() for i in range(NUM_MAF_RANGES)]
	for qq_i in range(NUM_MAF_RANGES):
		slice_indices = (len(variants) * qq_i//4, len(variants) * (qq_i+1)//NUM_MAF_RANGES - 1)
		qqs[qq_i]['maf_range'] = (variants[slice_indices[0]].maf, variants[slice_indices[1]].maf)
		neglog10_pvals = sorted((v.neglog10_pval for v in variants[slice(*slice_indices)]), reverse=True)
		qqs[qq_i]['count'] = len(neglog10_pvals)

		occupied_bins = set()
		for i, obs_neglog10_pval in enumerate(neglog10_pvals):
			exp_neglog10_pval = -math.log10( (i+0.5) / len(neglog10_pvals))
			exp_bin = int(exp_neglog10_pval / max_exp_neglog10_pval * NUM_BINS)
			obs_bin = int(obs_neglog10_pval / max_obs_neglog10_pval * NUM_BINS)
			occupied_bins.add( (exp_bin,obs_bin) )
		# print(sorted(occupied_bins))

		qq = []
		for exp_bin, obs_bin in occupied_bins:
			assert 0 <= exp_bin <= NUM_BINS, exp_bin
			assert 0 <= obs_bin <= NUM_BINS, obs_bin
			qq.append((
				exp_bin / NUM_BINS * max_exp_neglog10_pval,
				obs_bin / NUM_BINS * max_obs_neglog10_pval
			))
		qq = sorted(qq)

		qqs[qq_i]['qq'] = qq

	return qqs



##FOR PLINK AND RAREMETAL FILES
def make_qq(file_reader):
	neglog10_pvals = []
	#skip the header
	file_reader.skip_header()


	while True:
		#get a new line
		file_reader.get_line()
		#if we are at the end of file
		if file_reader.is_end():
			break
		pval = file_reader.get_pval()
		if pval == 'NA':
			continue
		pval = float(pval)

		neglog10_pvals.append(-math.log10(pval))


	neglog10_pvals = sorted(neglog10_pvals)



	# QQ
	max_exp_neglog10_pval = -math.log10(1/len(neglog10_pvals)) #expected
	max_obs_neglog10_pval = max(neglog10_pvals) #observed
	# print(max_obs_neglog10_pval, max_exp_neglog10_pval)

	occupied_bins = set()
	for i, obs_neglog10_pval in enumerate(neglog10_pvals):
		exp_neglog10_pval = -math.log10((len(neglog10_pvals)-i)/len(neglog10_pvals))
		exp_bin = int(exp_neglog10_pval / max_exp_neglog10_pval * NUM_BINS)
		obs_bin = int(obs_neglog10_pval / max_obs_neglog10_pval * NUM_BINS)
		occupied_bins.add( (exp_bin,obs_bin) )
	# print(sorted(occupied_bins))

	qq = []
	for exp_bin, obs_bin in occupied_bins:
		assert exp_bin <= NUM_BINS, exp_bin
		assert obs_bin <= NUM_BINS, obs_bin
		qq.append((
			exp_bin / NUM_BINS * max_exp_neglog10_pval,
			obs_bin / NUM_BINS * max_obs_neglog10_pval
		))
	qq = sorted(qq)

	# GC_value lambda
	median_neglog10_pval = neglog10_pvals[len(neglog10_pvals)//2]
	median_pval = 10 ** -median_neglog10_pval # I know, `10 ** -(-math.log10(pval))` is gross.
	gc_value_lambda = gc_value(median_pval)

	rv = {
		'qq': qq,
		'median_pval': median_pval,
		'gc_value_lambda': round_sig(gc_value_lambda, 5),
	}
	return rv



if __name__ == '__main__':
	epacts_filename = sys.argv[1]
	assert os.path.exists(epacts_filename)
	out_filename = sys.argv[2]
	
	file_reader = Data_reader.Data_reader.factory(epacts_filename, None)
	header = file_reader.check_header()

	rv = make_qq(file_reader)

	

	with open(out_filename, 'w') as f:
		json.dump(rv, f, sort_keys=True, indent=0)
	print('{} -> {}'.format(epacts_filename, out_filename))