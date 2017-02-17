import gzip
import Data_reader
##REQUIRES filename is a tabix file, names is the header of the file
##MODIFIES nothing
##EFFECTS finds the position of the minimum pvalue
def find_min_pvals(filename, filetype, num_minimums, region_buffer):
	#create a file reader from the file
	file_reader = Data_reader.Data_reader.factory(filename, filetype)
	#skip the header
	file_reader.skip_header()

	#create the minimums dictionary
	minimums = create_baseline_minimums(num_minimums)
	#find the highest of the minimums

	highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
	#loops through the lines in the file
	
	line = file_reader.get_line()
	while line != '#Genomic' or line != '':
		
		if line == '' or line.split()[0] == '#Genomic':
			break
		#if the pvalue is not available
		if file_reader.get_pval() == 'NA':
			line = file_reader.get_line()
			continue
		#if the pvalue is equal to the highest minimum, we do not add it to dictionary
		elif float(file_reader.get_pval()) >= highest_min:
			line = file_reader.get_line()
			continue
		#lastly, we must check other attributes of this pval if we want to add it to the dictionary
		else:
			#determine if this pvalue shares a region with another minimum
			shares_region, shared_index = index_of_shared_region(minimums, num_minimums, long(file_reader.get_pos()), region_buffer)
			#if it does share a region:
			if shares_region:
				#determine which is smaller, and place the smaller minimum in the list
				if float(file_reader.get_pval()) < minimums['value'][shared_index]:
					minimums = replace_minimum(minimums, long(file_reader.get_pos()), float(file_reader.get_pval()), int(file_reader.get_chrom()), shared_index)
					highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
				else:
					line = file_reader.get_line()
					continue
			#if it does not share a region, place replace the previous highest minimum with the new minimum
			else:
				minimums = replace_minimum(minimums, long(file_reader.get_pos()), float(file_reader.get_pval()), int(file_reader.get_chrom()), highest_min_index)
				highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
		line = file_reader.get_line()
	
	minimums = sort_minimums(minimums, num_minimums)
	
	return minimums


##REQUIRES minimums has at least two minimums
##MODIFIES minimums
##EFFECTS sorts (decreasing order) the dictionary of minimums based on pvalue
def sort_minimums(minimums, num_minimums):
	new_minimums = create_baseline_minimums(num_minimums)
	index = 0 
	for min in minimums['value']:
		best = find_min_of_mins(minimums)
		new_minimums['position'][index] = minimums['position'][best]
		new_minimums['value'][index] = minimums['value'][best]
		new_minimums['chromosome'][index] = minimums['chromosome'][best]
		minimums['value'][best] = 1
		index += 1
	
	return new_minimums


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


##REQUIRES: minimums is a dictionary of minimums
##MODIFIES nothing
##EFFECTS creats a top hits list
def create_hits(minimums):
	hits = []
	
	##create the hits list for flask
	for x in range(0, 10):
		chr = minimums['chromosome'][x]
		chr = str(chr)
		pos = minimums['position'][x]
		pos = str(pos)
		hits.append([chr + ":" + pos, chr + ":" + pos])
	return hits


##REQUIRES
##MODIFIES
##EFFECTS
def get_basic_region(filename, filetype):
	#create a file reader from the file
	file_reader = Data_reader.Data_reader.factory(filename, filetype)
	#skip the header
	file_reader.skip_header()

	#get a line
	line = file_reader.get_line()
	chrom = file_reader.get_chrom()
	position = file_reader.get_pos()

	return str(chrom) + ":" + str(position) + "-" + str(int(position) + 200000)