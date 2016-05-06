import gzip
from lz import *

header = check_header('assoc.q.lm.epacts.gz')
print get_position(header)
with gzip.open('assoc.q.lm.epacts.gz', 'rb') as f:
	next(f)
	for line in f:
		words = line.split()

print words[1]
position_column = get_position(header)
with gzip.open('assoc.q.lm.epacts.gz', 'rb') as f:
	#skip the header
	next(f)
	#loop through the lines
	for line in f:
		#split up the line into a list of dats
		data = line.split()
		#if the data is out of range
	
		if long(data[position_column]) < 1 or long(data[position_column]) >= 193797:
			continue
		#if data is in range
		print data[position_column]
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
