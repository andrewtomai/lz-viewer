import gzip
import tabix_solution
##REQUIRES filename is a tabix file, names is the header of the file
##MODIFIES nothing
##EFFECTS finds the position of the minimum pvalue
def find_min_pvals(filename, names, num_minimums, region_buffer):
	#get the column of the pvals
	pval_column = tabix_solution.get_column(names, "pvalue")
	#get the column of positions
	position_column = tabix_solution.get_column(names, "position")
	#get the column of the chromosome
	chrom_column = tabix_solution.get_column(names, "chr")

	#open the file
	with gzip.open(filename, 'rb') as f:
		#skip the header
		for line in f:
			line = line.split()
			if (line[0] == "#CHROM") | (line[0] == "CHR"):
				break
		

		#create the minimums dictionary
		minimums = create_baseline_minimums(num_minimums)
		#find the highest of the minimums

		highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
		#loops through the lines in the file
		for line in f:
			
			data = line.split()
			#if we have hit the end of a RAREMETAL file
			if data[0] == '#Genomic':
				break
			
			#if the data is not availalbe
			if data[pval_column] == 'NA':
				continue
			
			#if the current pvalue is greater or equal to the highest minimum, we do not add it to the dictionary
			elif float(data[pval_column]) >= highest_min:
				continue
			#if the current pvalue should (possibly) be added to the dictionary of mins
			else:
				shares_region, shared_index = index_of_shared_region(minimums, num_minimums, long(data[position_column]), region_buffer)
				if shares_region:
					if float(data[pval_column]) < minimums['value'][shared_index]:
						minimums = replace_minimum(minimums, long(data[position_column]), float(data[pval_column]), int(data[chrom_column]), shared_index)
						highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
					else:
						continue
				else:
					minimums = replace_minimum(minimums, long(data[position_column]), float(data[pval_column]), int(data[chrom_column]), highest_min_index)
					highest_min, highest_min_index = find_highest_min(minimums, num_minimums)
	
	return minimums


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