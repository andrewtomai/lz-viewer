import re
import math
import gzip
import zlib
from cStringIO import StringIO
from struct import *


class Data_reader(object):
	#create a data reader factory based on filetype
	def factory(file, type):
		#return the correct type of reader
		if type == None:
			type = get_filetype(file)
		if type == "EPACTS": return Epacts_reader(file)
		if type == "RAREMETAL": return Raremetal_reader(file)
		if type == "PLINK": return Plink_reader(file)
		assert 0, "Filetype not recognized: " + type
	factory = staticmethod(factory)

	##REQUIRES a line has been read
	##MODIFIES this
	##EFFECTS returns true if the end of file has been reached
	def is_end(self):
		if self.line == None:
			return True
		return False

	##REQUIRES nothing
	##MODIFIES nothing
	##EFFECTS default returns NULL because we do not know where the information is stored
	def get_ref(self):
		return None
	##REQUIRES nothing
	##MODIFIES nothing
	##EFFECTS default returns NULL because we do not know where the information is stored
	def get_alt(self):
		return None


	##REQUIRES file object has been initialized
	##MODIFIES nothing
	##EFFECTS reads a line from the fileobject
	def get_line(self):
		line = self.file_obj.readline()
		self.line = line
		if line == '' or line.split()[0] == "#Genomic":
			self.line = None
		return line

	##REQUIRES a line has been read
	##MODIFIES nothing
	##EFFECTS returns the pvalue from the line as string
	def get_pval(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		return line['pvalue']

	##REQUIRES a line has been read
	##MODIFIES nothing
	##EFFECTS returns the chromosome from the line as string
	def get_chrom(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		return line['chr']

	##REQUIRES a line has been read
	##MODIFIES nothing
	##EFFECTS returns the position from the line as string
	def get_pos(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		return line['position']


	##REQUIRES names is a header list
	##MODIFIES position_column
	##EFFECTS returns the column that denotes position
	def get_column(self, key):
		column = 0
		for name in self.header:
			if name == key:
				break
			column += 1
		return column 
	
	##REQUIRES filename has an associated tabix file
	##MODIFIES nothing
	##EFFECTS creates a file object from the tabix file
	def open_tabix(self, filename):
		tabix = gzip.open(filename + '.tbi', 'rb')
		return tabix

	##REQUIRES line is not the header line
	##MODIFIES nothing
	##EFFECTS returns a dictionary with the header and the given value
	def make_line_dict(self, line):
		line = line.split()
		data = {}
		for x in range(0, len(line)):
			data[self.header[x]] = line[x]
		return data



	##REQUIRES names is header list, current_column is < len(names), dict has keys from names
	##MODIFIES dict
	##EFFECTS adds 'datum' to the dict
	def format_datum(self, datum):

			#then the datum is a scientifically notated float
		if (re.search('[a-zA-Z]', datum.replace('e', '')) == None) and (datum.count('e') == 1):
			return (float(datum))

		elif re.search('[a-zA-Z:]', datum) != None:
			#then its the name of a variant
			relevant = re.search('^\S+:\d+_[ACTG]/[ACTG]', datum)
			if relevant:
				return relevant.group()
			else:
				return (str(datum))
		elif "." in datum:
			#then its a float
			return (float(datum))
		else:
			#it must be a long int or an int
			return (long(datum))
	

	##REQUIRES filename is a gzipped file
	##MODIFIES stdout
	##EFFECTS reads the header of filename, returns a list of names.
	def check_header(self):
		#open the gzip file
		with gzip.open(self.filename) as f:
			#if the filetype is raremetal
			if self.filetype == "RAREMETAL":
				#loop through the lines until we get to the header line
				for line in f:
					#split up the line and look for the word #CHROM
					words = line.split()
					if words[0] == "#CHROM":
						break
			elif (self.filetype == "EPACTS") or (self.filetype == "PLINK"):
				text = f.readline()
				#split this text into words
				words = text.split()
			else:
				text = f.readline()
				words = text.split()
			#initialize a list of names of columns
			names = []
			#loops through the columns to create list of names
			for word in words:
				if word == "BEGIN":
					names.append("position")
				elif word == "MARKER_ID":
					names.append("id")
				elif word == "#CHROM":
					names.append("chr")
				elif (word == "POS") or (word == "BP"):
					names.append("position")
				elif word == "P":
					names.append("pvalue")
				else: 
					names.append(word.lower())
			if (self.filetype == "RAREMETAL") or (self.filetype == "PLINK"):
				names.append("id")
		return names



	##REQUIRES tabix is a .tbi file object
	##MODIFIES nothing
	##EFFECTS creates a dictionary of the tabix header fields
	def get_tabix_fields(self, tabix):
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
	def find_block(self, start, end, chrom_in):

		#open the tabix (BGZF compressed) file
		tabix = self.open_tabix(self.filename)
		#get the fields for the tabix file
		field_dict = self.get_tabix_fields(tabix)
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




	##REQUIRES file_obj is freshly opened file
	##MODIFIES nothing
	##EFFECTS skips to the header of the file
	def skip_header(self, *args):
		if len(args) == 1:

			for line in file_obj:
					line = line.split()
					if (line[0] == "#CHROM") or (line[0] == "CHR"):
						break
			return file_obj
		else:
			for line in self.file_obj:
				line = line.split()
				if (line[0] == "#CHROM") or (line[0] == "CHR"):
					self.line = ' '.join(line)
					break
			
				
			return None

	##REQUIRES start and end are proper numbers
	##MODIFIES nothing
	##EFFECTS skips to the correct block of info, and decompressed and returns a file object
	def skip_to_relevant(self, start, end, chrom_in):
		assert start >= 0 and end >= start, "Invalid start or end positions"
		
		##OPEN tHE BGZIPPED FILE
		#get offset information
		offset_within_block, block_offset, end_offset = self.find_block(start, end, chrom_in)

		#skip to the block offset of compressed file
		with open(self.filename, 'rb') as compressed:
			compressed.seek(block_offset)

			#if the end position is in the last block
			if end_offset == -1:
				block = compressed.read()
				
			#otherwise, only read the least amount of data to decompress
			else:
				block = compressed.read(end_offset - block_offset + 1)

		#create a decompressed StringIO from teh compressed data
		decompressed = gzip.GzipFile(self.filename, fileobj=StringIO(block))

		#skip to the relevant data within the block
		decompressed.seek(offset_within_block)

		#if there was no offset, we need to skip the header
		if(offset_within_block == 0 and block_offset == 0):
			#skip the header, 
			decompressed = self.skip_header(decompressed)
			

		return decompressed



class Epacts_reader(Data_reader):
	#initialize the reader with filename.
	def __init__(self, filename):
		self.filename = filename
		self.filetype = "EPACTS"
		self.file_obj = gzip.open(filename)
		self.line = None
		self.header = self.check_header()

	##REQUIRES a variant line has been read using get_line
	##MODIFIES self
	##EFFECTS returns the ref allele
	def get_ref(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		chr1, pos1, ref, alt, opt_info = re.match(r'([^:]+):([0-9]+)_([-ATCG]+)/([-ATCG]+)(?:_(.+))?', line['id']).groups()
		return ref

	##REQUIRES a variant line has been read using get_line
	##MODIFIES self
	##EFFECTS returns the alt allele
	def get_alt(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		chr1, pos1, ref, alt, opt_info = re.match(r'([^:]+):([0-9]+)_([-ATCG]+)/([-ATCG]+)(?:_(.+))?', line['id']).groups()
		return alt


	##REQUIRES start, end, and chrom_in are all valid
	##MODIFIES nothing
	##EFFECTS creates a dictionary based on the range positions
	def create_dict_from_file(self, start, end, chrom_in):
		
		data_dict = {}
		for name in self.header:
			data_dict[name] = []

		decompressed = self.skip_to_relevant(start, end, chrom_in)
		for line in decompressed:
			data = self.make_line_dict(line)
			#if the data is not yet in range
			if long(data['position']) < start:
				continue
			#if the data is past the range
			elif long(data['position']) >= end:
				break

			#if pvalue isnt available
			if data['pvalue'] == 'NA':
				continue

			#if data is in range
			for key in data:
				datum = self.format_datum(data[key])

				data_dict[key].append(datum)
		return data_dict


class Raremetal_reader(Data_reader):
	#initialize the reader with filename
	def __init__(self, filename):
		self.filename = filename
		self.filetype = "RAREMETAL"
		self.file_obj = gzip.open(filename)
		self.line = None
		self.header = self.check_header()


	##REQUIRES a variant line has been read using get_line
	##MODIFIES self
	##EFFECTS returns the ref allele
	def get_ref(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		
		return line['ref']

	##REQUIRES a variant line has been read using get_line
	##MODIFIES self
	##EFFECTS returns the alt allele
	def get_alt(self):
		assert self.line != None, "A line hasn't been read yet."
		line = self.make_line_dict(self.line)
		
		return line['alt']


	##REQUIRES start, end, and chrom_in are all valid
	##MODIFIES nothing
	##EFFECTS creates a dictionary based on the range positions
	def create_dict_from_file(self, start, end, chrom_in):
		data_dict = {}
		for name in self.header:
			data_dict[name] = []

		decompressed = self.skip_to_relevant(start, end, chrom_in)
		for line in decompressed:
			data = self.make_line_dict(line)

			#if the data is not yet in range
			if long(data['position']) < start:
				continue
			#if the data is past the range
			elif long(data['position']) >= end:
				break

			#if pvalue isnt available
			if data['pvalue'] == 'NA':
				continue

			#if the data is in range
			id = data['chr'] + ':' + data['position'] + '_' + data['ref'] + '/' + data['alt']
			data['id'] = id

			for key in data:
				datum = self.format_datum(data[key])

				data_dict[key].append(datum)

		return data_dict


class Plink_reader(Data_reader):
	#initialize the reader with the filename and type
	def __init__(self, filename):
		self.filename = filename
		self.filetype = "PLINK"
		self.file_obj = gzip.open(filename)
		self.line = None
		self.header = self.check_header()


	##REQUIRES start, end, and chrom_in are all valid
	##MODIFIES nothing
	##EFFECTS creates a dictionary based on the range positions
	def create_dict_from_file(self, start, end, chrom_in):
		data_dict = {}
		for name in self.header:
			data_dict[name] = []

		decompressed = self.skip_to_relevant(start, end, chrom_in)
		for line in decompressed:
			
			data = self.make_line_dict(line)

			#if the data is not yet in range
			if long(data['position']) < start:
				continue
			#if the data is past the range
			elif long(data['position']) >= end:
				
				break

			#if pvalue isnt available
			if data['pvalue'] == 'NA':
				continue

			id = data['chr'] + ':' + data['position']
			data['id'] = id

			for key in data:
				datum = self.format_datum(data[key])

				data_dict[key].append(datum)

		return data_dict


##REQUIRES the file is EPACTS, PLINK, or RAREMETAL
##MODIFIES nothing
##EFFECTS returns the filetype name based on header, or NULL if unknown
def get_filetype(filename):
	with gzip.open(filename) as f:
		line = f.readline()
		line = line.split()
		if line[0] == "#CHROM":
			filetype = "EPACTS"
		elif line[0] == "CHR":
			filetype = "PLINK"
		elif line[0] == "##ProgramName=RareMetalWorker":
			filetype = "RAREMETAL"
		else:
			filetype = None

	return filetype




