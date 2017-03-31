import ast
import get_options
import minimum_solution
import manhattan_binning
import Minimum_Thread
import qq_to_json
import Data_reader
import json
import sys
import os
import site
import threading
from werkzeug.contrib.cache import SimpleCache
from flask import Flask, jsonify, request, render_template, url_for, Response, send_from_directory, redirect


#########################################################################################################
#--------------------------------------Flask initialization---------------------------------------------#	
lz_app = Flask(__name__, static_url_path='')

#if the user wants to shutdown the server
def shutdown_server():
	func = request.environ.get('werkzeug.server.shutdown')
	if func is None:
		raise RuntimeError('Not running with Werkzeug Server')
	func()

@lz_app.route('/shutdown/', methods=['POST', 'GET'])
def shutdown():
	shutdown_server()
	return 'Server shutting down...'

@lz_app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
	response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
	return response

@lz_app.route('/static/<path:path>')
def send_js(path):
	return lz_app.send_static_file(path)

##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the manhattan plot template
@lz_app.route('/')
def home():
	if range_opt:
		return redirect("/lz/" + range_opt)
	return redirect(url_for("manhattan"))
	

##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the manhattan plot template
@lz_app.route('/manhattan/')
def manhattan():
	if range_opt:
		default_range = range_opt
	
	else:
		default_range = minimum_range
	return render_template("manhattan.html", default_range=default_range)


##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the lz plot template
@lz_app.route('/lz/<region>')
def lz_region(region):
	print(region)
	if region == "null":
		region = minimum_solution.get_basic_region(filename, filetype)
	return render_template("lz-association-viewer.html", region=region, filetype=filetype)


##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the qq plot template
@lz_app.route('/QQ/')
def QQ_plot():
	return redirect(url_for("QQ_plot_unbin"))



##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS make a route for the binned json for the QQ plot
@lz_app.route('/QQ-unbin/')
def QQ_plot_unbin():
	if range_opt:
		default_range = range_opt
	else: 
		default_range = minimum_range
	return render_template("qq-unbin.html", default_range=default_range)



##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS make a route for the binned json for the QQ plot
@lz_app.route('/QQ-bin/')
def QQ_plot_bin():
	if range_opt:
		default_range = range_opt
	else: 
		default_range = minimum_range

	#if the filetype isnt EPACTS, go to the unbinned version
	if filetype is not "EPACTS":
		return redirect(url_for("QQ_plot_unbin"))

	return render_template("qq-bin.html", default_range=default_range)




##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
@lz_app.route('/api/lz-results/', methods=['GET'])
def api_lz():
	#get the filter query string
	filter = request.args.get('filter')
	
	#if we have a query string, we need to update the chromosome and position information
	if filter:
		chromosome, start, end = get_options.parse_query(filter)
		chrom_pos_dict['chromosome'] = chromosome
		chrom_pos_dict['start'] = start
		chrom_pos_dict['end'] = end

	#if there is no filter, then we know the user is trying to access api endpoint	
	else:
		chromosome = chrom_pos_dict['chromosome']
		start = chrom_pos_dict['start']
		end = chrom_pos_dict['end']

	#gather the data
	file_reader = Data_reader.Data_reader.factory(filename, filetype)
	data = file_reader.create_dict_from_file(start, end, chromosome)
		

	#format the dictionary according to the portal API
	object = get_options.format_data(data)

	return jsonify(object)




##REQUIRES nothing
##MODIFIES lz_app
## returns the JSON api for the "top hits"
@lz_app.route('/api/hits/', methods=['GET'])
def api_hits():
	minimums = minimumThread.join()
	hits = minimum_solution.create_hits(minimums)

#find the minimum pvalue, and assign to the global
	if minimum:
		min_index = minimum_solution.find_min_of_mins(minimums)
		minimum_chromosome = minimums['chromosome'][min_index]
		minimum_position = minimums['position'][min_index]

		global minimum_range
		left = minimum_position-100000
		right = minimum_position+100000
		minimum_range = str(minimum_chromosome) + ":" + str(left) + "-" + str(right)

	return jsonify(hits)




##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS returns the JSON api for manhattan plots
@lz_app.route('/api/manhattan/', methods=['GET'])
def api_manhattan():
	#cache.clear()
	rv = cache.get('manhattan')
	if rv is None:
		#if not cached on server, check if its cached on local disk
		if os.path.isfile(filename + '.manhattan.json'):
			cache_fileobj = open(filename + '.manhattan.json', 'r')
			rv = cache_fileobj.read()
			rv = ast.literal_eval(rv)
		
		#if it is not cached on the local disc, create the manhattan plot
		else:
			file_reader = Data_reader.Data_reader.factory(filename, filetype)
			variant_bins, unbinned_variants = manhattan_binning.bin_variants(file_reader)

			rv = {
				'variant_bins': variant_bins,
				'unbinned_variants': unbinned_variants,
				}
			cache_fileobj = open(filename + '.manhattan.json', 'w')
			cache_fileobj.write(json.dumps(rv))

		cache.set('manhattan', rv)
	
	return jsonify(rv)
	



##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS returns the JSON api for QQ plots
@lz_app.route('/api/QQ/', methods=['GET'])
def api_qq():
	#check if we saved a cache from a previous run
	
	bin = request.args.get('bin')
	rv = cache.get('qq-'+bin)

	if rv is None:
		#check if data is cached on the local disk!
		if (bin == 'true' or bin == 'True') and os.path.isfile(filename + '.qq.bin.json'):
			cache_fileobj = open(filename + '.qq.bin.json', 'r')
			rv = cache_fileobj.read()

		elif (bin == 'false' or bin == 'False') and os.path.isfile(filename + '.qq.unbin.json'):
			cache_fileobj = open(filename + '.qq.unbin.json', 'r')
			rv = cache_fileobj.read()

		else:
			#we need to create the data!
			file_reader = Data_reader.Data_reader.factory(filename, filetype)
			#if we are looking for the binned solution:
			if bin == 'true' or bin == 'True':
				rv = qq_to_json.make_qq_stratified(file_reader, True)
				rv = json.dumps(rv);
				#write to file
				cache_fileobj = open(filename + '.qq.bin.json', 'w')
				cache_fileobj.write(rv)
			else:
				rv = qq_to_json.make_qq_stratified(file_reader, False)
				rv = json.dumps(rv);
				#write to file
				cache_fileobj = open(filename + '.qq.unbin.json', 'w')
				cache_fileobj.write(rv)

	#rv should have data! set the cache, and return the data
	cache.set('qq-'+bin, rv)
	resp = Response(response=rv, status=200, mimetype="application/json")
	return resp

	


#----------------------------------------------------------------------------------------------------#
######################################################################################################

	
def main():
	#check the input arguments
	arguments = get_options.check_options()
	global filename
	filename = arguments["filename"]
	global port_number
	port_number = arguments["port_number"]
	global minimum 
	minimum = arguments["minimum"]
	global range_opt 
	range_opt = arguments["range"]
	global host
	host = arguments["host"]
	
	############################################################################################
	
	#try to open up the file 
	try:
		trail = open(filename, 'r')

	except IOError as detail:
		print >> sys.stderr, detail
		print >> sys.stderr,  "The file should be in the current working directory: "
		print >> sys.stderr,  os.getcwd()
		sys.exit(1)

	try:
		trail = open(filename + '.tbi', 'r')

	except IOError as detail:
		print >> sys.stderr, detail
		print >> sys.stderr,  "This program requires a tabix-indexed file, with the tabix file located in the current working directory: "
		print >> sys.stderr,  os.getcwd()
		sys.exit(1)
	#########################################################################################
	
	#Create a file object from said file		
	global filetype
	filetype = get_options.get_filetype(arguments, filename)
	file_reader = Data_reader.Data_reader.factory(filename, filetype)

	#get the header of the file
	global header
	header = file_reader.header
	
	#create a thread to find the minimum solution
	global minimumThread 
	minimumThread = Minimum_Thread.Minimum_Thread(filename, filetype, 10, 102400)
	minimumThread.start()
	global minimums 
	minimums = None
	
	#create a default minimum range
	global minimum_range
	minimum_range = "null"

	#create a dictionary that contains the most recent called chromosome and positions
	global chrom_pos_dict
	chrom_pos_dict = {'chromosome': 11, 'start': 1, 'end': 2}
	
	#create a simple cache for the session
	global cache
	cache = SimpleCache()
	

	#run the flask webserver
	lz_app.run(host = host, port = port_number, debug=True, threaded=True)

	#(main)#
if __name__ == '__main__':
	main()