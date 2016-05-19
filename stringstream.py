from cStringIO import StringIO
import zlib
import io
import re
import gzip

filename = 'assoc.q.lm.epacts.gz'
with gzip.open(filename, 'rb') as f:
	outfile = open("assoc.txt", 'w')
	for line in f:
		outfile.write(line + '\n')


