################################################################################################
#---------------------------------------helper functions---------------------------------------#

import pysam
import argparse
import json
import subprocess
import re 
import gzip
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
##EFFECTS reads the header of filename using terminal command "tabix -H", returns a list of names.
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
			names.append("variant")
		else: 
			names.append(word.lower())
	return names




######################################################################################################
#########--------------------------MY SOLUTION FOR TABIX-----------------------------#################
##REQUIRES filename is a tabix file, names is the list of header info from the file
##MODIFIES stdout
##EFFECTS reads the data of filename using the terminal command "tabix [range]", returns a dictionary of data
def gather_data_gzip(filename, names, start, end):
	#initialize the dictionary with proper keys
	data_dict = { }
	for name in names:
		data_dict[name] = []

	#get the column of positions
	position_column = get_position(names)
	#set the current column to 0
	current_column = 0
	#get the number of columns
	num_columns = len(names)
	#open the gzip file

	with gzip.open(filename, 'rb') as f:
		#skip the header
		next(f)
		#loop through the lines
		for line in f:
			#split up the line into a list of dats
			data = line.split()
			#if the data is out of range
			if data[position_column] <= start or data[position_column] > end:
				continue
			else: #if data is in range
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



	#get the number of columns
	#num_columns = len(names)
	#split the data into datum
	#data = text.split()
	#initialize the current column
	#current_column = 0

	#gather the data
	#for datum in data:
	#	if current_column == num_columns:
	#		current_column = 0
	#	if datum == '\t':
	#		continue
	#	elif datum == '\n':
	#		continue

		#add the datum to the correct key using 'add_datum'
	#	data_dict = add_datum(names, current_column, datum, data_dict)
	#	current_column += 1
	#return data_dict




##REQUIRES names is header list, current_column is < len(names), dict has keys from names
##MODIFIES dict
##EFFECTS adds 'datum' to the dict
def add_datum(names, column, datum, dict):
	
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

##REQUIRES names is a header list
##MODIFIES position_column
##EFFECTS returns the column that denotes position
def get_position(names):
	column = 0
	for name in names:
		if name == "position":
			break
		column += 1
	return column 





		

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
	for row in file.fetch(11, 193154, 500000):
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





##REQUIRES data is a dictionary 
##MODIFIES data
##EFFECTS adds necessary 'lastPage' key to the dictionary  
def format_data(data):
	new_dict = {'data' : data, 'lastPage' : None}
	return new_dict



##Flask initialization	
lz_app = Flask(__name__, static_url_path='')
@lz_app.route('/')
def home():
	return lz_app.send_static_file('lz-association-viewer.html')
@lz_app.route('/api', methods=['GET'])
##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
def api():
	#check the input arguments
	arguments = check_options()
	filename = arguments["filename"]
	
	#get the header of the file
	header = check_header(filename)
	data = gather_data(filename, header, 193154, 193797)

	#format the dictionary according to the portal API
	object = format_data(data)

	return json.dumps(object)

#----------------------------------------------------------------------------------------------------#
######################################################################################################

#(main)#
if __name__ == '__main__':

	#check the input arguments
	arguments = check_options()

	port_number = arguments["port_number"]

	#run the flask webserver
	lz_app.run(port = port_number)