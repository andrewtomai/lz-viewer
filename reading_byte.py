from struct import *
import gzip
import time
with gzip.open("assoc.q.lm.epacts.gz.tbi", "rb") as tabix:
	fields = ['magic', 'n_ref', 'format', 'col_seq', 'col_beg', 'col_end', 'meta', 'skip', 'l_nm', 'names'];
	field_dict = { };
	outfile = open("data.txt", 'w')


	for x in range(0, 9):
		information = tabix.read(4);
		if x == 0:
			field_dict[fields[x]] = information;
		else:
			information = unpack('i', information)[0];
			field_dict[fields[x]] = information;
	field_dict['names'] = tabix.read(field_dict['l_nm']);
	outfile.write( str(field_dict) + '\n')
	outfile.write('\n')
	n_bin = unpack('i', tabix.read(4))[0];
	outfile.write( "n_bin: " + str(n_bin) + "\n") 
	for x in range(0, n_bin):

		bin = unpack('I', tabix.read(4))[0];
		n_chunk = unpack('i', tabix.read(4))[0];
		tabix.read(16*n_chunk)
			
	n_intv =  unpack('i', tabix.read(4))[0];
	outfile.write("n_intv: %d" % (n_intv))
	for x in range(0, n_intv):
		voffset =  unpack('H', tabix.read(2))[0]
		offset = unpack('I', tabix.read(4))[0]
		highbits = unpack('H', tabix.read(2))[0]

		outfile.write("voffset: %d   offset: %d, \n " % (voffset, offset))


	print tabix.read()
	