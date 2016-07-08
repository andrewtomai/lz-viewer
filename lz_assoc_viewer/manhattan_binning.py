#!/usr/bin/env python2

'''
This script takes two arguments:
- an input filename (the output of epacts with a single phenotype)
- an output filename (please end it with `.json`)
It creates a json file which can be used to render a Manhattan plot.
'''

from __future__ import print_function, division, absolute_import

import os.path
import sys
import gzip
import re
import json
import math
import collections
import Data_reader

BIN_LENGTH = int(3e6)
NEGLOG10_PVAL_BIN_SIZE = 0.05 # Use 0.05, 0.1, 0.15, etc
NEGLOG10_PVAL_BIN_DIGITS = 2 # Then round to this many digits
BIN_THRESHOLD = 1e-4 # pvals less than this threshold don't get binned.

def round_sig(x, digits):
	return 0 if x==0 else round(x, digits-1-int(math.floor(math.log10(abs(x)))))
assert round_sig(0.00123, 2) == 0.0012
assert round_sig(1.59e-10, 2) == 1.6e-10



Variant = collections.namedtuple('Variant', 'chrom pos ref alt maf pval beta sebeta'.split())
def parse_variant_line(file_reader):
	
	if file_reader.get_pval() == 'NA':
		return None
	else:
		chrom, pos, pval, ref, alt =  file_reader.get_chrom(), file_reader.get_pos(), file_reader.get_pval(), file_reader.get_ref(), file_reader.get_alt()
	
		return {'chrom' : chrom, 'pos': long(pos), 'ref' : ref, 'alt' : alt, 'pval' : float(pval)}

def rounded_neglog10(pval):
	return round(-math.log10(pval) // NEGLOG10_PVAL_BIN_SIZE * NEGLOG10_PVAL_BIN_SIZE, NEGLOG10_PVAL_BIN_DIGITS)

def get_pvals_and_pval_extents(pvals):
# expects that NEGLOG10_PVAL_BIN_SIZE is the distance between adjacent bins.
	pvals = sorted(pvals)
	extents = [[pvals[0], pvals[0]]]
	for p in pvals:
		if extents[-1][1] + NEGLOG10_PVAL_BIN_SIZE * 1.1 > p:
			extents[-1][1] = p
		else:
			extents.append([p,p])
	rv_pvals, rv_pval_extents = [], []
	for (start, end) in extents:
		if start == end:
			rv_pvals.append(start)
		else:
			rv_pval_extents.append([start,end])
	return (rv_pvals, rv_pval_extents)

def bin_variants(file_reader):
	bins = []
	unbinned_variants = []

	prev_chrom, prev_pos = -1, -1

	file_reader.skip_header()
	
	while True:
		file_reader.get_line()
		
		if file_reader.is_end():
			
			break
		
		
		variant = parse_variant_line(file_reader)
		if variant is None: continue
		prev_chrom, prev_pos = variant['chrom'], variant['pos']

		if variant['pval'] < BIN_THRESHOLD:
			unbinned_variants.append({
				'chrom': variant['chrom'],
				'pos': variant['pos'],
				'ref': variant['ref'],
				'alt': variant['alt'],
				'pval': round_sig(variant['pval'], 2),
			})
				
		else:
			if len(bins) == 0 or variant['chrom'] != bins[-1]['chrom']:
				#we need a new bin, starting with this variant
				 bins.append({
                    'chrom': variant['chrom'],
                    'startpos': variant['pos'],
                    'neglog10_pvals': set(),
                })
			elif variant['pos'] > bins[-1]['startpos'] + BIN_LENGTH:
                # We need a new bin following the last one.
				bins.append({
					'chrom': variant['chrom'],
					'startpos': bins[-1]['startpos'] + BIN_LENGTH,
					'neglog10_pvals': set(),
				})
			bins[-1]['neglog10_pvals'].add(rounded_neglog10(variant['pval']))
		
	bins = [b for b in bins if len(b['neglog10_pvals']) != 0]
	for b in bins:
		b['neglog10_pvals'], b['neglog10_pval_extents'] = get_pvals_and_pval_extents(b['neglog10_pvals'])
		b['pos'] = int(b['startpos'] + BIN_LENGTH/2)
		del b['startpos']
	
	return bins, unbinned_variants

if __name__ == '__main__':
	epacts_filename = sys.argv[1]
	assert os.path.exists(epacts_filename)
	out_filename = sys.argv[2]
	#assert os.path.exists(os.path.dirname(out_filename))

	file_reader = Data_reader.Data_reader.factory(epacts_filename, None)
	header = file_reader.check_header()
	
	variant_bins, unbinned_variants = bin_variants(file_reader)
	
	
	

	rv = {
		'variant_bins': variant_bins,
		'unbinned_variants': unbinned_variants,
	}

	# Avoid getting killed while writing dest_filename, to stay idempotent despite me frequently killing the program
	with open(out_filename, 'w') as f:
		json.dump(rv, f, sort_keys=True, indent=0)
	print('{} -> {}'.format(epacts_filename, out_filename))