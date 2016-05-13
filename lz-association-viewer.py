################################################################################################
#---------------------------------------helper functions---------------------------------------#

import pysam
import math
import argparse
import re 
import gzip
from struct import *
from flask import Flask, jsonify, request, render_template, url_for



##REQUIRES: nothing
##MODIFIES: file_name, port_number
##EFFECTS: parses the input arguments to set the correct file_name and port_number
def check_options():
	#initialize the parser with description
	parser = argparse.ArgumentParser(description="Open input_file and view it on the specified or default port")
	#add the port argument
	parser.add_argument("-p", "--port", help="Specify a port at which to view results", type=int)
	#add the filename argument
	parser.add_argument("filename", type=str, help="Provide the name of the file to be jsonified")
	#parse the arguments
	args = parser.parse_args()
	#set the filename
	file = args.filename
	#if a port was specified
	if args.port:
		port_number = args.port

	else: #otherwise the default port is 5000
		port_number = 5000
	#return a dictionary including the filename and port number
	return {'filename' : file, 'port_number' : port_number }


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
def find_block(filename, start):

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
		for x in range(0, field_dict['names'].count('\x00')):
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
			#set the current block number to 0
			current_block = 0
			#set the previous block_offset to block_offset
			prev_off = block_offset
			for x in range(1, index):
				offset_within_block = unpack('<H', tabix.read(2))[0]
				block_offset = unpack('<I', tabix.read(4))[0]
				highbits = unpack('H', tabix.read(2))[0]
				if(block_offset > prev_off):
					current_block += 1
				prev_off = block_offset

		return offset_within_block, current_block



def make_voffset( current_block, offset_within_block ):
	return (current_block << 16) | offset_within_block


##REQUIRES filename is a tabix file, names is the list of header info from the file
##MODIFIES nothing
##EFFECTS reads the data of filename, returns a dictionary of data
def gather_data_gzip(filename, names, start, end):
	
	offset_within_block, current_block = find_block(filename, start)
	offset = make_voffset(current_block, offset_within_block)
	
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
	
	with gzip.open(filename, 'rb') as f:
		#skip to the correct block 
		f.seek(offset)
		if offset == 0:
			f.next()
		
		#loop through the lines
		for line in f:
			
			#split up the line into a list of datums
			data = line.split()
			
			#if the data is not yet in range
			if long(data[position_column]) < start:
				print line
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







#################################################################################################
#############-------------------PYSAM SOLUTION FOR TABIX---------------------####################
##REQUIRES: file_name is a TABIX file
##MODIFIES: file
##EFFECTS: opens the TABIX formatted file using pysam 
def open_file(filename):
	#open the tabix file
	file = pysam.TabixFile(filename)
	return file




##REQUIRES file is a pysam tabix file
##MODIFIES data
##EFFECTS creates a data dictionary by iterating through the tabix file
def create_data(file):
	#setting values of keys
	variant = "variant"
	position = "position"
	pvalue = "pvalue"
	#initializing the data dictionary
	data = { pvalue : [], position : [], variant : [] }
	#looping through the specific region of the tabix file
	for row in file.fetch(11, 113850000, 114150000):
		#splitting the row up into individual words

		words = row.split()
		#adding the data to the dictionary
		data[variant].append(words[3])
		if words[8] == 'NA': data[pvalue].append(None)
		else: data[pvalue].append(float(words[8]))

		data[position].append(long(words[1]))
		
	return data
##################---------------------------------------------------------------#################
##################################################################################################

##REQUIRES filename is a tabix file, names is the header of the file
##MODIFIES nothing
##EFFECTS finds the position of the minimum pvalue
def find_min_pval(filename, names):
	#get the column of the pvals
	pval_column = get_column(names, "pvalue")
	#get the column of positions
	position_column = get_column(names, "position")
	#open the file
	with gzip.open(filename, 'rb') as f:
		next(f)
		#set baselines for current position and min
		current_min = 1
		current_position = -1
		#find the minimum element's position
		for line in f:
			data = line.split()
			
			if data[pval_column] == "NA":
				continue
			
			if float(data[pval_column]) < current_min:
				current_min = float(data[pval_column])
				current_position = long(data[position_column])
	return current_position


##REQUIRES data is a dictionary 
##MODIFIES data
##EFFECTS adds necessary 'lastPage' key to the dictionary  
def format_data(data):
	new_dict = {'data' : data, 'lastPage' : None}
	return new_dict



##Flask initialization	
lz_app = Flask(__name__, static_url_path='')
@lz_app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
	response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
	return response

@lz_app.route('/')
def home():

	return render_template("lz-association-viewer.html", port=port_number, region="11:114mb+150kb")
@lz_app.route('/api/results/', methods=['GET'])
##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
def api():
	
	start = 113850000
	end = 114150000
	#if we, for some reason, did not find a minimum
	if minimum_position == -1:
		start = 113850000
		end = 114150000
	#gather the data
	data = gather_data_gzip(filename, header, start, end)

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

	#get the header of the file
	global header
	header = check_header(filename)
	#find the minimum pvalue
	global minimum_position
	minimum_position = find_min_pval(filename, header)

	
	#run the flask webserver
	lz_app.run(port = port_number)