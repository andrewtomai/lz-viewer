
import get_options
import minimum_solution
import manhattan_binning
import qq_to_json
import Data_reader

import json

from flask import Flask, jsonify, request, render_template, url_for, Response, send_from_directory, redirect

#########################################################################################################
#--------------------------------------Flask initialization---------------------------------------------#	
lz_app = Flask(__name__, static_url_path='')
@lz_app.after_request
def after_request(response):
	response.headers.add('Access-Control-Allow-Origin', '*')
	response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
	response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
	return response

@lz_app.route('/static/<path:path>')
def send_js(path):
    return send_from_directory('static', path)

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
	return render_template("manhattan.html", port=port_number, default_range=default_range)


##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the lz plot template
@lz_app.route('/lz/<region>')
def lz_region(region):
	return render_template("lz-association-viewer.html", port=port_number, region=region, hits=hits, filetype=filetype)


##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS renders the qq plot template
@lz_app.route('/QQ/')
def QQ_plot():
	if range_opt:
		default_range = range_opt
	else: 
		default_range = minimum_range
	return render_template("qq.html", port=port_number, default_range=default_range)

@lz_app.route('/api/lz-results/', methods=['GET'])

##REQUIRES object is a dictionary
##MODIFIES lz_app
##EFFECTS displays a json objects at route '/api'
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
##EFFECTS returns the JSON api for manhattan plots
@lz_app.route('/api/manhattan/', methods=['GET'])
def api_manhattan():
	file_reader = Data_reader.Data_reader.factory(filename, filetype)
	variant_bins, unbinned_variants = manhattan_binning.bin_variants(file_reader)
	rv = {
		'variant_bins': variant_bins,
		'unbinned_variants': unbinned_variants,
		}
	return jsonify(rv)
	
##REQUIRES nothing
##MODIFIES lz_app
##EFFECTS returns the JSON api for QQ plots
@lz_app.route('/api/QQ/', methods=['GET'])
def api_qq():
	bin = request.args.get('bin')
	file_reader = Data_reader.Data_reader.factory(filename, filetype)
	
	if (bin == 'true' or bin == 'True') and filetype == 'EPACTS' :
		
		rv = qq_to_json.make_qq_stratified(file_reader)
	
	else:
		rv = qq_to_json.make_qq(file_reader)

	resp = Response(response=json.dumps(rv), status=200, mimetype="application/json")
	return resp


#----------------------------------------------------------------------------------------------------#
######################################################################################################

#(main)#
if __name__ == '__main__':
	

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
	
	
	global filetype
	filetype = get_options.get_filetype(arguments, filename)
	
	file_reader = Data_reader.Data_reader.factory(filename, filetype)

	#get the header of the file
	global header
	header = file_reader.header
	

	#find the minimums
	minimums = minimum_solution.find_min_pvals(filename, filetype, 10, 102400)
	
	
	#create a list of 'hits'
	global hits
	hits = minimum_solution.create_hits(minimums)
	
#find the minimum pvalue
	if minimum:
		min_index = minimum_solution.find_min_of_mins(minimums)

		minimum_chromosome = minimums['chromosome'][min_index]
		minimum_position = minimums['position'][min_index]

		global minimum_range
		left = minimum_position-100000
		right = minimum_position+100000
		minimum_range = str(minimum_chromosome) + ":" + str(left) + "-" + str(right)

		
	#create a dictionary that contains the most recent called chromosome and positions
	global chrom_pos_dict
	chrom_pos_dict = {'chromosome': 11, 'start': 1, 'end': 2}
	
	

	
	#run the flask webserver
	lz_app.run(port = port_number, debug=True, threaded=True)