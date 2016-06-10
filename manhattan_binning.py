from __future__ import division, absolute_import

import gzip
import re
import json
import math
import collections
import make_manhattan_json
import tabix_solution


def average_and_add_bin(dict, bin, header):
	for key in bin:
		
		
		#if the type is a number
		if type(bin[key][0]) == int or type(bin[key][0]) == long or type(bin[key][0]) ==  float:
			average = sum(bin[key])/len(bin[key])
			dict[key].append(average)
		
			#if the item is an id
		elif re.search('^\S+:\d+_[ACTG]/[ACTG]', bin[key]):
			#it is a id
			continue

		#if the key is the chromosome
		elif key == 'chr':
			dict[key].append(bin[key][0])

		id_string = str(bin['chr'][0]) + ':' + str((bin['position'][0] + bin['position'][-1])/2)
		dict['id'].append(id_string)



##REQUIRES variant_lines does not include the header
##MODIFIES nothing
##EFFECTS thins data data for manhattan plots
def bin_variants(variant_lines, header):
	#initialize the dictionary with proper keys
	data_dict = { }
	for name in names:
		data_dict[name] = []

	current_bin = {}

	#get relevant columns 
	pval_column = tabix_solution.get_column(header, "pvalue")
	chrom_column = tabix_solution.get_column(header, "chr")
	pos_column = tabix_solution.get_column(header, "pos")



	

	#iterate through all lines
	for variant_line in variant_lines:
		#split the variant line
		variant_line = variant_line.split()


		
		#if the pvalue is significant enough
		if variant_line[pval_column] < make_manhattan_json.BIN_THRESHOLD:
			
			#loop through the data and add it to the dictionary
			column = 0
			for data in variant_line:
				tabix_solution.add_datum(header, column, data, data_dict)
				column += 1

		else:
			if len(current_bin) == 0 or variant_line[chrom_column] != current_bin['chr']:
				#we need a new bin that starts with this variant
				#first we add the previous bin to the dictionary
				data_dict = average_and_add_bin(data_dict, current_bin, header)
				current_bin = {}
				
				#add data to the new bin
				column = 0 
				for data in variant_line:
					tabix_solution.add_datum(header, column, data, current_bin)
					column += 1


			elif long(variant_line[pos_column]) > current_bin['position'][0] + make_manhattan_json.BIN_LENGTH:
				#we beed a new bin following the last one
				#first we add the previous bin to the dictionary
				data_dict =  average_and_add_bin(data_dict, current_bin, header)
				current_bin = {}

				#add the data to the new bin
				column = 0
				for data in variant_line:
					tabix_solution.add_datum(header, column, data, current_bin)
					column += 1

	if len(current_bin) != 0:
		data_dict = average_and_add_bin(data_dict, current_bin, header)




		
