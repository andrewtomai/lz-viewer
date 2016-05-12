from struct import *
import gzip
with gzip.open("assoc.q.lm.epacts.gz.tbi", "rb") as tabix:
	fields = ['magic', 'n_ref', 'format', 'col_seq', 'col_beg', 'col_end', 'meta', 'skip', 'l_nm', 'names'];
	field_dict = { };
	

	for x in range(0, 9):
		information = tabix.read(4);
		if x == 0:
			field_dict[fields[x]] = information;
		else:
			information = unpack('i', information)[0];
			field_dict[fields[x]] = information;
	field_dict['names'] = tabix.read(field_dict['l_nm']);
	
	n_bin = unpack('i', tabix.read(4))[0];
	print "n_bin: %d" % (n_bin)
	for x in range(0, n_bin):

		bin = unpack('I', tabix.read(4))[0];
		n_chunk = unpack('i', tabix.read(4))[0];
		#print "bin: %d" % (bin)
		#print "n_chunk: %d" % (n_chunk)
		for x in range(0, n_chunk):
			
			cnk_beg = unpack('L', tabix.read(8))[0];
			cnk_end = unpack('L', tabix.read(8))[0];
		#	print "cnk_beg: %d" % (cnk_beg)
		#	print "cnk_end: %d" % (cnk_end)
		
	n_intv =  unpack('i', tabix.read(4))[0];
	print "n_intv: %d" % (n_intv)
	for x in range(0, n_intv):
		ioff =  unpack('<L', tabix.read(8))[0];
		print "ioff: %d" % (ioff)
	