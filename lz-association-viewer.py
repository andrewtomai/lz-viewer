import get_options
import tabix_solution
import minimum_solution
from flask import Flask, jsonify, request, render_template, url_for

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
	data = tabix_solution.gather_data_gzip(filename, filetype, header, start, end, chromosome)

	#format the dictionary according to the portal API
	object = get_options.format_data(data)

	return jsonify(object)

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
	
	
	#get the header of the file
	global header
	header = get_options.check_header(filename, filetype)
	

	#find the minimums
	minimums = minimum_solution.find_min_pvals(filename, header, 10, 102400)
	
	
	#create a list of 'hits'
	global hits
	hits = minimum_solution.create_hits(minimums)
	
#find the minimum pvalue
	if minimum:
		min_index = minimum_solution.find_min_of_mins(minimums)

		minimum_chromosome = minimums['chromosome'][min_index]
		minimum_position = minimums['position'][min_index]

		range_opt = str(minimum_chromosome) + ":" + str(minimum_position) + "+150kb" 
		#create a dictionary that contains the most recent called chromosome and positions
		global chrom_pos_dict
		chrom_pos_dict = {'chromosome': 11, 'start': 1, 'end': 2}
			

	
	#run the flask webserver
	lz_app.run(port = port_number, debug=True, threaded=True)