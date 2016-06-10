import re
import math
import gzip
import zlib
from cStringIO import StringIO
from struct import *
import binascii

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

		#then the datum is a scientifically notated float
	elif (re.search('[a-zA-Z]', datum.replace('e', '')) == None) & (datum.count('e') == 1):
		dict[names[column]].append(float(datum))

	elif re.search('[a-zA-Z:]', datum) != None:
		#then its the name of a variant
		relevant = re.search('^\S+:\d+_[ACTG]/[ACTG]', datum)
		if relevant:
			dict[names[column]].append(relevant.group())
		else:
			dict[names[column]].append(str(datum))
	elif "." in datum:
		#then its a float
		dict[names[column]].append(float(datum))
	else:
		#it must be a long int or an int
		dict[names[column]].append(long(datum))
	return dict






##REQUIRES filename has an associated tabix file
##MODIFIES nothing
##EFFECTS creates a file object from the tabix file
def open_tabix(filename):
	tabix = gzip.open(filename + '.tbi', 'rb')
	return tabix







##REQUIRES tabix is a .tbi file object
##MODIFIES nothing
##EFFECTS creates a dictionary of the tabix header fields
def get_tabix_fields(tabix):
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

	return field_dict
	




##REQUIRES filename is a tabix file
##MODIFIES IO
##EFFECTS finds the byte number of the best position location to start reading data
def find_block(filename, start, end, chrom_in):

	#open the tabix (BGZF compressed) file
	tabix = open_tabix(filename)
	#get the fields for the tabix file
	field_dict = get_tabix_fields(tabix)
	#get the list of chromosomes from the fields dictionary
	chromosomes = field_dict['names'].split('\x00')
		
	for chromosome in chromosomes:
		if chromosome =='':
			continue
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


		


		#if we are at the correct chromosome
		if int(chromosome) == chrom_in:	
			#initialize a counter to count the current number of 16kb intervals passed
		
			#find the correct 'ioff' index 
			index = int(math.floor(start/16384))
			#find the correct index of 'ioff' for the end position
			end_index = int(math.floor(end/16384))
			
			#unpack the first offset
			offset_within_block =  unpack('H', tabix.read(2))[0]
			block_offset = unpack('I', tabix.read(4))[0]

			#should have no value
			highbits = unpack('H', tabix.read(2))[0]
			#loop through the ioff's, until we get the index of the end position
			for x in range(1, end_index):
				#if we havent gotten to the start index yet
				if x < index:
					offset_within_block = unpack('<H', tabix.read(2))[0]
					block_offset = unpack('<I', tabix.read(4))[0]
					highbits = unpack('H', tabix.read(2))[0]

				#if we have gotten or surpassed the start index
				elif x >= index:
					lowbits = tabix.read(2)
					if lowbits == '':
						end_actual_offset = -1
						break
					else: 
						lowbits = unpack('<H', lowbits)[0]
						end_actual_offset = unpack('<I', tabix.read(4))[0]
						highbits = unpack('H', tabix.read(2))[0]
			
			#now find the offset of the block one after the end block
			end_offset = end_actual_offset	
				
			while end_offset == end_actual_offset:
				lowbits = tabix.read(2)
				if lowbits == '':
					end_offset = -1
					break
				else:
					lowbits = unpack('H', lowbits)[0]
					end_offset = unpack('<I', tabix.read(4))[0]
					highbits = unpack('H', tabix.read(2))[0]
				#if the end_offset is empty, then the last block is also the block containing the end
				if end_offset == '':
					end_offset = -1
					break

		#if we are not at the correct chromosome
		else:
			for x in range(0, n_intv):
				tabix.read(8)

		#if we are at the correct chromosome
		if int(chromosome) == chrom_in:
			break
			
	
	return offset_within_block, block_offset, end_offset


	##REQUIRES stream is a gzip compressed string
	##MODIFIES 
	##EFFECTS decrompresses the gzipped string, and returns the decompressed string
def stream_gzip_decompress(stream):
	# skip the header of the compressed string
	dec = zlib.decompressobj(32 + zlib.MAX_WBITS)
	#decrompress the stream
	rv = dec.decompress(stream)
	#make sure the string isnt NULL
	if rv:
		yield rv




##REQUIRES filename is a tabix file, names is the list of header info from the file
##MODIFIES nothing
##EFFECTS reads the data of filename, returns a dictionary of data
def gather_data_gzip(filename, filetype, names, start, end, chrom_in):
	
		
	#initialize the dictionary with proper keys
	data_dict = { }
	for name in names:
		data_dict[name] = []

	#get the column of positions
	position_column = get_column(names, "position")
	pval_column = get_column(names, "pvalue")
	#if the filetype is a raremetal file, we need to get a few more columns
	if filetype == "RAREMETAL":
		chrom_column = get_column(names, "chr")
		ref_column = get_column(names, "ref")
		alt_column = get_column(names, "alt")
	if filetype == "PLINK":
		chrom_column = get_column(names, "chr")
	#set the current column to 0
	current_column = 0
	#get the number of columns
	num_columns = len(names)



	#open the gzip file

	#Get the offset information for the range we want
	offset_within_block, block_offset, end_offset = find_block(filename, start, end, chrom_in)
	
	#skip to the block_offset of the compressed file
	with open(filename, 'rb') as compressed:

		compressed.seek(block_offset)

		#if the end position is in the last block
		if end_offset == -1:
			block = compressed.read()
		#otherwise, only read the least amount of data to decompress
		else: 
			block = compressed.read(end_offset - block_offset + 1)
			
	#create a decompressed StringIO from the compressed data
	
	decompressed = gzip.GzipFile(filename, fileobj=StringIO(block))
	
	#skip to the relevant data within the block
	decompressed.seek(offset_within_block)
	
	#if there was no offset, we need to skip the header
	if offset_within_block == 0 & block_offset == 0:
		#skip the header
		for line in decompressed:
			line = line.split()
			if (line[0] == "#CHROM") | (line[0] == "CHR"):
				break

	#loop through the lines
	for line in decompressed:
		
		#split up the line into a list of datums
		data = line.split()
		
		#make sure that data is a complete line
		if (filetype == "RAREMETAL") | (filetype == "PLINK"):
			if len(data) != (len(names) - 1):
				break

		#if the data is not yet in range
		if long(data[position_column]) < start:
			continue
		#if the data is past the range
		elif long(data[position_column]) >= end:
			break
		#if the pvalue is not available
		if data[pval_column] == 'NA':
			continue
		#if data is in range
		if filetype == "RAREMETAL":
			id = data[chrom_column] + ":" + data[position_column] + "_" + data[ref_column] + "/" + data[alt_column]
			data.append(id)
		elif filetype == "PLINK":
			id = data[chrom_column] + ":" + data[position_column]
			data.append(id)

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