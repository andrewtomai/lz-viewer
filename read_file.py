##extracting data from a tabix file, using PYSAM
import pysam
import argparse

##REQUIRES: nothing
##MODIFIES: file_name, port_number
##EFFECTS: parses the input arguments to set the correct file_name and port_number
def check_options():
	parser = argparse.ArgumentParser(description="Open input_file and view it on the specified or default port")
	parser.add_argument("-p", "--port", help="Specify a port at which to view results", type=int)
	parser.add_argument("filename", type=str, help="Provide name of the file to be jsonified")
	args = parser.parse_args()
	filename = args.filename
	if args.port:
		port_number = args.port
	else:
		port_number = 5000
	return {'filename' : filename, 'port_number' : port_number }


##REQUIRES: file_name is a TABIX file
##MODIFIES: file
##EFFECTS: opens the TABIX formatted file using pysam 
def open_file(file_name):
	file = pysam.TabixFile(file_name)
	return file

def create_json(file):

arguments = check_options()



file = open_file(file_name)

#for row in file.fetch(11, 114550452, 115067678):
#	print (str(row))

