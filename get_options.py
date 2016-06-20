import re 
import math
import argparse
import gzip

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

	filetype = parser.add_mutually_exclusive_group()
	#add the "RAREMETAL" argument
	filetype.add_argument("-R", "--RAREMETAL", help="Specifies if RAREMETAL results are being used", action="store_true")
	#add the "PLINK" argument
	filetype.add_argument("-P","--PLINK", help="Specifies if PLINK results are being used", action="store_true")
	#add the "EPACTS" argument
	filetype.add_argument("-E", "--EPACTS", help="Specifies if EPACTS results are being used", action="store_true")
	#add the range argument
	parser.add_argument("-r", "--range", type=str, help="Provide the range of positions to grab of format: [CHROMOSOME #]:[START]-[END]")
	#add the manhattan argument
	parser.add_argument("-m", "--manhattan", help="Specify if the plot should be a manhattan plot", action="store_true")
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
	return {'filename' : file, 'port_number' : port_number, 'range' : range, 'minimum' : minimum, 'EPACTS': args.EPACTS, 'RAREMETAL' : args.RAREMETAL, 'PLINK' : args.PLINK, 'manhattan' : args.manhattan}



##REQUIRES: arguments is a list of arguments given by check_arguments, filename is valid
##MODIFIES nothing
##EFFECTS detects the filetype of the file to be read
def get_filetype(arguments, filename):
	EPACTS = arguments["EPACTS"]
	RAREMETAL = arguments["RAREMETAL"]
	PLINK = arguments["PLINK"]

	if EPACTS:
		filetype = "EPACTS"
	elif RAREMETAL:
		filetype = "RAREMETAL"
	elif PLINK:
		filetype = "PLINK"
	else:
		with gzip.open(filename) as f:
			line = f.readline()
			line = line.split()
			if line[0] == "#CHROM":
				filetype = "EPACTS"
			elif line[0] == "CHR":
				filetype = "PLINK"
			else:
				filetype = "RAREMETAL"
	return filetype






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
