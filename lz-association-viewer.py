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


##REQUIRES filename is a tabix file
##MODIFIES IO
##EFFECTS finds the byte number of the best position location to start reading data
def find_block(filename, start, end, chrom_in):

	#open the tabix (BGZF compressed) file
	with gzip.open(filename + '.tbi', 'rb') as tabix:
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
			#unpack the first offset
			offset_within_block =  unpack('H', tabix.read(2))[0]
			block_offset = unpack('I', tabix.read(4))[0]

			#should have no value
			highbits = unpack('H', tabix.read(2))[0]
			#find the correct start offsets
			for x in range(1, index):
				offset_within_block = unpack('<H', tabix.read(2))[0]
				block_offset = unpack('<I', tabix.read(4))[0]
				highbits = unpack('H', tabix.read(2))[0]
			
			end_index = int(math.floor(end/16384))
			
				

			#if we are at the correct chromosome
			if int(chromosome) == chrom_in:
				break


		return offset_within_block, block_offset


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
	#set the current column to 0
	current_column = 0
	#get the number of columns
	num_columns = len(names)
	#open the gzip file


	#Get the offset information for the range we want
	offset_within_block, block_offset = find_block(filename, start, end, chrom_in)

	#skip to the block_offset of the compressed file
	with open(filename, 'rb') as compressed:
		compressed.seek(block_offset)
		block = compressed.read()

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









##REQUIRES filename is a tabix file, names is the header of the file
##MODIFIES nothing
##EFFECTS finds the position of the minimum pvalue
def find_min_pval(filename, names):
	#get the column of the pvals
	pval_column = get_column(names, "pvalue")
	#get the column of positions
	position_column = get_column(names, "position")
	#get the column of the chromosome
	chrom_column = get_column(names, "chr")

	#open the file
	with gzip.open(filename, 'rb') as f:
		next(f)
		#set baselines for current position and min
		current_min = 1
		current_position = -1
		current_chrom = -1
		#find the minimum element's position
		for line in f:
			data = line.split()
			
			if data[pval_column] == "NA":
				continue
			
			if float(data[pval_column]) < current_min:
				current_min = float(data[pval_column])
				current_position = long(data[position_column])
				current_chrom = int(data[chrom_column])
	return current_chrom, current_position











##REQUIRES data is a dictionary 
##MODIFIES data
##EFFECTS adds necessary 'lastPage' key to the dictionary  
def format_data(data):
	new_dict = {'data' : data, 'lastPage' : None}
	return new_dict







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

	return render_template("lz-association-viewer.html", port=port_number, region=range_opt)
@lz_app.route('/api/results/', methods=['GET'])
##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
def api():
	

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

	global chromosome
	global start
	global end

#find the minimum pvalue
	if minimum:

		minimum_chromosome, minimum_position = find_min_pval(filename, header)
		range_opt = str(minimum_chromosome) + ":" + str(minimum_position) + "+150kb" 
		chromosome = minimum_chromosome
		start = minimum_position-150000
		end = minimum_position+150000

	else: #if we have a custom range
		range_list = re.split('[:-]', range_opt)
		chromosome = int(range_list[0])
		start = long(range_list[1])
		end = long(range_list[2])
	
	
	
	#run the flask webserver
	lz_app.run(port = port_number)