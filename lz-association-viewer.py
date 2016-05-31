################################################################################################
#---------------------------------------helper functions---------------------------------------#

import pysam
import math
import argparse
import re 
import gzip
import zlib
import time
from cStringIO import StringIO
from struct import *
from flask import Flask, jsonify, request, render_template, url_for




##REQUIRES: nothing
##MODIFIES: file_name, port_number
##EFFECTS: parses the input arguments to set the correct file_name and port_number
def check_options():
	#initialize the parser with description
	parser = argparse.ArgumentParser(description="Open input_file and view it on the specified or default port, with a default data range based off of the minimum p-value")
	#add the port argument
	parser.add_argument("-p", "--port", help="Specify a port at which to view results", type=int)
	#add the filename argument
	parser.add_argument("filename", type=str, help="Provide the name of the file to be graphed")
	
	#add the range argument
	parser.add_argument("-r", "--range", type=str, help="Provide the range of positions to grab of format: [CHROMOSOME #]:[START]-[END]")
	#parse the arguments
	args = parser.parse_args()
	#set the filename
	file = args.filename
	#if a port was specified
	if args.port:
		port_number = args.port

	else: #otherwise the default port is 5000
		port_number = 5000
		
	minimum = False 
	range = None
	#if a range was specified
	if args.range:
		range = args.range
	else: #otherwise the default range is set to be based off of the minimum p-value
		minimum = True


	#return a dictionary including the filename and port number
	return {'filename' : file, 'port_number' : port_number, 'range' : range, 'minimum' : minimum}


##REQUIRES filename is a tabix file
##MODIFIES stdout
##EFFECTS reads the header of filename, returns a list of names.
def check_header(filename):
	#open the gzip file
	with gzip.open(filename) as f:
		text = f.readline()
	#split this text into words
	words = text.split()
	#initialize a list of names of columns
	names = []
	#loops through the columns to create list of names
	for word in words:
		if word == "BEGIN":
			names.append("position")
		elif word == "MARKER_ID":
			names.append("id")
		elif word == "#CHROM":
			names.append("chr")
		else: 
			names.append(word.lower())
	return names




######################################################################################################
#########--------------------------MY SOLUTION FOR TABIX-----------------------------#################
##REQUIRES names is a header list
##MODIFIES position_column
##EFFECTS returns the column that denotes position
def get_column(names, key):
	column = 0
	for name in names:
		if name == key:
			break
		column += 1
	return column 


##REQUIRES names is header list, current_column is < len(names), dict has keys from names
##MODIFIES dict
##EFFECTS adds 'datum' to the dict
def add_datum(names, column, datum, dict):
	#the datum should be stored as NULL
	if datum == "NA":
		dict[names[column]].append(None)
		#then the datum is a scientifically notated float
	elif (re.search('[a-zA-Z]', datum.replace('e', '')) == None) & (datum.count('e') == 1):
		dict[names[column]].append(float(datum))
	elif re.search('[a-zA-Z]', datum) != None:
		#then its the name of a variant
		dict[names[column]].append(datum)
	elif "." in datum:
		#then its a float
		dict[names[column]].append(float(datum))
	else:
		#it must be a long int or an int
		dict[names[column]].append(long(datum))
	return dict






##REQUIRES filename has an associated tabix file
##MODIFIES nothing
##EFFECTS creates a file object from the tabix file
def open_tabix(filename):
	tabix = gzip.open(filename + '.tbi', 'rb')
	return tabix







##REQUIRES tabix is a .tbi file object
##MODIFIES nothing
##EFFECTS creates a dictionary of the tabix header fields
def get_tabix_fields(tabix):
	fields = ['magic', 'n_ref', 'format', 'col_seq', 'col_beg', 'col_end', 'meta', 'skip', 'l_nm', 'names']
	field_dict= { }

	for x in range(0, 9):
		#read in four bytes of data
		information = tabix.read(4)
		#if this is the first read, the information is of type char
		if x == 0:
			field_dict[fields[x]] = information
		#otherwise, the information is an integer
		else:
			information = unpack('i', information)[0];
			field_dict[fields[x]] = information

	#'names' depends on the number that represents 'l_nm'
	field_dict['names'] = tabix.read(field_dict['l_nm'])

	return field_dict
	




##REQUIRES filename is a tabix file
##MODIFIES IO
##EFFECTS finds the byte number of the best position location to start reading data
def find_block(filename, start, end, chrom_in):

	#open the tabix (BGZF compressed) file
	tabix = open_tabix(filename)
	#get the fields for the tabix file
	field_dict = get_tabix_fields(tabix)
	#get the list of chromosomes from the fields dictionary
	chromosomes = field_dict['names'].split('\x00')
		
	for chromosome in chromosomes:
			
		#number of bins
		n_bin = unpack('i', tabix.read(4))[0]

		#unpack all bins, using the n_bin variable
		for x in range(0, n_bin):
			bin = unpack('I', tabix.read(4))[0]
			
			#find the number of chunks in the bin
			n_chunk = unpack('i', tabix.read(4))[0]
			#use the number of chunks to unpack all chunks
			for x in range(0, n_chunk):
				#chunk begins:
				cnk_beg = unpack('Q', tabix.read(8))[0]
				#chunk ends:
				cnk_end = unpack('Q', tabix.read(8))[0]

		#find the number of intervals
		n_intv = unpack('i', tabix.read(4))[0]


		#initialize a counter to count the current number of 16kb intervals passed
		
		#find the correct 'ioff' index 
		index = int(math.floor(start/16384))
		#find the correct index of 'ioff' for the end position
		end_index = int(math.floor(end/16384))
			
		#unpack the first offset
		offset_within_block =  unpack('H', tabix.read(2))[0]
		block_offset = unpack('I', tabix.read(4))[0]

		#should have no value
		highbits = unpack('H', tabix.read(2))[0]
		
			
		#loop through the ioff's, until we get the index of the end position
		for x in range(1, end_index):
			#if we havent gotten to the start index yet
			if x < index:
				offset_within_block = unpack('<H', tabix.read(2))[0]
				block_offset = unpack('<I', tabix.read(4))[0]
				highbits = unpack('H', tabix.read(2))[0]

			#if we have gotten or surpassed the start index
			elif x >= index:
				lowbits = tabix.read(2)
				if lowbits == '':
					end_actual_offset = -1
					break
				else: 
					lowbits = unpack('<H', lowbits)[0]
					end_actual_offset = unpack('<I', tabix.read(4))[0]
					highbits = unpack('H', tabix.read(2))[0]
			
		#now find the offset of the block one after the end block
		end_offset = end_actual_offset	
				
		while end_offset == end_actual_offset:
			lowbits = tabix.read(2)
			if lowbits == '':
				end_offset = -1
				break
			else:
				lowbits = unpack('H', lowbits)[0]
				end_offset = unpack('<I', tabix.read(4))[0]
				highbits = unpack('H', tabix.read(2))[0]
			#if the end_offset is empty, then the last block is also the block containing the end
			if end_offset == '':
				end_offset = -1
				break

		#if we are at the correct chromosome
		if int(chromosome) == chrom_in:
			break
			
		
	return offset_within_block, block_offset, end_offset


	##REQUIRES stream is a gzip compressed string
	##MODIFIES 
	##EFFECTS decrompresses the gzipped string, and returns the decompressed string
def stream_gzip_decompress(stream):
	# skip the header of the compressed string
	dec = zlib.decompressobj(32 + zlib.MAX_WBITS)
	#decrompress the stream
	rv = dec.decompress(stream)
	#make sure the string isnt NULL
	if rv:
		return rv




##REQUIRES filename is a tabix file, names is the list of header info from the file
##MODIFIES nothing
##EFFECTS reads the data of filename, returns a dictionary of data
def gather_data_gzip(filename, names, start, end, chrom_in):
	
		
	#initialize the dictionary with proper keys
	data_dict = { }
	for name in names:
		data_dict[name] = []

	#get the column of positions
	position_column = get_column(names, "position")
	pval_column = get_column(names, "pvalue")
	#set the current column to 0
	current_column = 0
	#get the number of columns
	num_columns = len(names)
	#open the gzip file


	#Get the offset information for the range we want
	offset_within_block, block_offset, end_offset = find_block(filename, start, end, chrom_in)
	
	#skip to the block_offset of the compressed file
	with open(filename, 'rb') as compressed:
		compressed.seek(block_offset)
		#if the end position is in the last block
		if end_offset == -1:
			block = compressed.read()
		#otherwise, only read the least amount of data to decompress
		else: 
			block = compressed.read(end_offset - block_offset + 1)

	#create a decompressed StringIO from the compressed data
	decompressed = StringIO(stream_gzip_decompress(block))
	#skip to the relevant data within the block
	decompressed.seek(offset_within_block)
	#if there was no offset, we need to skip the header
	if offset_within_block == 0 & block_offset == 0:
		decompressed.next()	 

	#loop through the lines
	for line in decompressed:
			
		#split up the line into a list of datums
		data = line.split()
		#if the data is not yet in range
		if long(data[position_column]) < start:
			continue
		#if the data is past the range
		elif long(data[position_column]) >= end:
			break
		#if the pvalue is not available
		if data[pval_column] == 'NA':
			continue
		#if data is in range
		for datum in data:
			if current_column == num_columns:
				current_column = 0			
			if datum == '\t':
				continue
			elif datum == '\n':
				continue
			else:
				data_dict = add_datum(names, current_column, datum, data_dict)
				current_column += 1

		
	return data_dict








		

#############----------------------------------------------------------------------##############
#################################################################################################




#################################################################################################
########-------------------------------solution for minimums----------------------------#########




##REQUIRES filename is a tabix file, names is the header of the file
##MODIFIES nothing
##EFFECTS finds the position of the minimum pvalue
def find_min_pvals(filename, names, num_minimums, region_buffer):
	#get the column of the pvals
	pval_column = get_column(names, "pvalue")
	#get the column of positions
	position_column = get_column(names, "position")
	#get the column of the chromosome
	chrom_column = get_column(names, "chr")

	#open the file
	with gzip.open(filename, 'rb') as f:
		#skip the header
		next(f)
		#create the minimums dictionary
		minimums = create_baseline_minimums(num_minimums)
		#find the highest of the minimums

		highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
		#loops through the lines in the file
		for line in f:
			data = line.split()
			if data[pval_column] == 'NA':
				continue
			#if the current pvalue is greater or equal to the highest minimum, we do not add it to the dictionary
			elif float(data[pval_column]) >= highest_min:
				continue
			#if the current pvalue should (possibly) be added to the dictionary of mins
			else:
				shares_region, shared_index = index_of_shared_region(minimums, num_minimums, long(data[position_column]), region_buffer)
				if shares_region:
					if float(data[pval_column]) < minimums['value'][shared_index]:
						minimums = replace_minimum(minimums, long(data[position_column]), float(data[pval_column]), int(data[chrom_column]), shared_index)
						highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
					else:
						continue
				else:
					minimums = replace_minimum(minimums, long(data[position_column]), float(data[pval_column]), int(data[chrom_column]), highest_min_index)
					highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
	
	return minimums


##REQUIRES minimums has at least 1 minimum
##MODIFIES minimums
##EFFECTS returns an updated dictionary of minimums
def replace_minimum(minimums, position, pvalue, chromosome, index):
	minimums['position'][index] = position
	minimums['value'][index] = pvalue
	minimums['chromosome'][index] = chromosome
	return minimums


##REQUIRES minimums has at least 1 minimum
##MODIFIES nothing
##EFFECTS returns a bool and a index, denoting that the current position is within a certain buffer region of another minimum
def index_of_shared_region(minimums, num_minimums, position, region_buffer):
	for x in range(0, num_minimums):
		position_diff = abs( position - minimums['position'][x] )
		if position_diff < region_buffer:
			return True, x
	return False, -1


##REQUIRES minimums has a least one 'minimum' in it 
##MODIFIES
##EFFECTS returns the highest minimum and the index it is stored at
def find_highest_min(minimums, num_minimums):
	current_max = 0
	
	for x in range(0, num_minimums):
		if minimums['value'][x] > current_max:
			current_max = minimums['value'][x]
			current_position = x
	return current_max, current_position


##REQUIRES num_minimums is > 0
##MODIFIES nothing
##EFFECTS creates a minimums dictionary, including position, value and chromosome
def create_baseline_minimums(num_minimums):
	minimums = {'position' : [], 'value' : [], 'chromosome' : [] }
	for x in range( 0 , num_minimums ):
		minimums['position'].append(-1000000)
		minimums['value'].append(1)
		minimums['chromosome'].append(0)
	
	return minimums 


##REQUIRES minimums is a dictionary of minimums
##MODIFIES nothing
##EFFECTS finds the index of the minimum of the minimums
def find_min_of_mins(minimums):
	current_min = 1
	counter = 0
	for min in minimums['value']:

		if current_min > min:
			current_min = min
			current_position = counter

		counter += 1
	return current_position

########-----------------------------------------------------------------------------------------########
#########################################################################################################



##REQUIRES data is a dictionary 
##MODIFIES data
##EFFECTS adds necessary 'lastPage' key to the dictionary  
def format_data(data):
	new_dict = {'data' : data, 'lastPage' : None}
	return new_dict


##REQUIRES filter is a query string 
##MODIFIES nothing
##EFFECTS returns the value of the given parameter from the filter string
def parse_query(filter):
	#split the filter into a list
	filter = filter.split()
	#initialize a counter to 0
	counter = 0
	#loop through the filter list
	for piece in filter:
		#if we find the word chromosome, we know that the nubmer of the chromosome is 2 after that
		if piece == 'chromosome':
			counter_chrom = counter + 2
		#denotes a start position
		elif piece == 'ge':
			counter_start = counter + 1
		#denotes an end position
		elif piece == 'le':
			counter_end = counter + 1
		counter += 1

	
	chromosome = int(filter[counter_chrom].strip("'"))
	
	start = long(filter[counter_start])
	
	end = long(filter[counter_end])
	return chromosome, start, end





#########################################################################################################
#--------------------------------------Flask initialization---------------------------------------------#	
lz_app = Flask(__name__, static_url_path='')
@lz_app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
	response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
	return response

@lz_app.route('/')
def home():
	
	return render_template("lz-association-viewer.html", port=port_number, region=range_opt, hits=hits)
@lz_app.route('/api/results/', methods=['GET'])
##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
def api():
	
	#get the filter query string
	filter = request.args.get('filter')
	
	#if we have a query string, we need to update the chromosome and position information
	if filter:
		chromosome, start, end = parse_query(filter)
		chrom_pos_dict['chromosome'] = chromosome
		chrom_pos_dict['start'] = start
		chrom_pos_dict['end'] = end

	#if there is no filter, then we know the user is trying to access api endpoint	
	else:
		chromosome = chrom_pos_dict['chromosome']
		start = chrom_pos_dict['start']
		end = chrom_pos_dict['end']

	#gather the data
	data = gather_data_gzip(filename, header, start, end, chromosome)

	#format the dictionary according to the portal API
	object = format_data(data)

	return jsonify(object)

#----------------------------------------------------------------------------------------------------#
######################################################################################################

#(main)#
if __name__ == '__main__':
	

	#check the input arguments
	arguments = check_options()
	global filename
	filename = arguments["filename"]
	global port_number
	port_number = arguments["port_number"]
	global minimum 
	minimum = arguments["minimum"]
	global range_opt 
	range_opt = arguments["range"]


	#get the header of the file
	global header
	header = check_header(filename)


	#find the minimums
	minimums = find_min_pvals(filename, header, 10, 102400)
	
	#create a list of 'hits'
	global hits
	hits = []
	
	##create the hits list for flask
	for x in range(0, 10):
		chr = minimums['chromosome'][x]
		chr = str(chr)
		pos = minimums['position'][x]
		pos = str(pos)
		hits.append([chr + ":" + pos, chr + ":" + pos])
	
#find the minimum pvalue
	if minimum:
		min_index = find_min_of_mins(minimums)

		minimum_chromosome = minimums['chromosome'][min_index]
		minimum_position = minimums['position'][min_index]

		range_opt = str(minimum_chromosome) + ":" + str(minimum_position) + "+150kb" 
	#create a dictionary that contains the most recent called chromosome and positions
	global chrom_pos_dict
	chrom_pos_dict = {'chromosome': 11, 'start': 1, 'end': 2}
	
	
	
	#run the flask webserver
	lz_app.run(port = port_number, debug=False, threaded=True)