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
		outfile.write("bin: " + str(bin) + "\n")
		outfile.write("n_chunk: " + str(n_chunk) + "\n")
		for x in range(0, n_chunk):
			
			cnk_beg = unpack('l', tabix.read(8))[0];
			cnk_end = unpack('Q', tabix.read(8))[0];
			outfile.write("cnk_beg: %d \n" % (cnk_beg))
			outfile.write("cnk_end: %d \n" % (cnk_end))
			
	n_intv =  unpack('i', tabix.read(4))[0];
	outfile.write("n_intv: %d" % (n_intv))
	for x in range(0, n_intv):
		ioff =  unpack('l', tabix.read(8))[0];
		outfile.write("ioff: %d \n" % (ioff))


	print tabix.read()
	