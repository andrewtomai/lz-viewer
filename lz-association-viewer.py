################################################################################################
#---------------------------------------helper functions---------------------------------------#

import pysam
import argparse
import json
from flask import Flask, jsonify

##REQUIRES: nothing
##MODIFIES: file_name, port_number
##EFFECTS: parses the input arguments to set the correct file_name and port_number
def check_options():
	#initialize the parser with description
	parser = argparse.ArgumentParser(description="Open input_file and view it on the specified or default port")
	#add the port argument
	parser.add_argument("-p", "--port", help="Specify a port at which to view results", type=int)
	#add the filename argument
	parser.add_argument("filename", type=str, help="Provide name of the file to be jsonified")
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




##REQUIRES: file_name is a TABIX file
##MODIFIES: file
##EFFECTS: opens the TABIX formatted file using pysam 
def open_file(file_name):
	#open the tabix file
	file = pysam.TabixFile(file_name)
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
		data[position].append(int(words[1]))
		data[pvalue].append(int(words[8]))
	print data
	return data

##REQUIRES data is a dictionary 
##MODIFIES data
##EFFECTS adds necessary 'lastPage' key  
def format_data(data):
	new_dict = {'data' : data, 'lastPage' : None}
	return new_dict



##Flask initialization	
lz_app = Flask(__name__)
@lz_app.route('/')
def home():
	return "hello world!"
@lz_app.route('/api', methods=['GET'])
##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
def api():
	#check the input arguments
	arguments = check_options()
	filename = arguments["filename"]
	
	#open the specified file
	file = open_file(filename)

	#create the dictionary from file
	data = create_data(file)

	#format the dictionary according to the portal API
	object = format_data(data)

	return jsonify(object)

#----------------------------------------------------------------------------------------------------#
######################################################################################################

#(main)#
if __name__ == '__main__':

	#check the input arguments
	arguments = check_options()
	#filename = arguments["filename"]
	port_number = arguments["port_number"]
	#open the specified file
	file = open_file(filename)
	#create the dictionary from file
	data = create_data(file)
	#format the dictionary according to the portal API
	object = format_data(data)
	#run the flask webserver
	#lz_app.run(port = port_number)