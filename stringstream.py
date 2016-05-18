from cStringIO import StringIO
import zlib
import io
import re

filename = 'test_file.txt'
with open(filename, 'rb') as f:

	file = f.read()
	hype = f.read(1)
	if hype == '':
		print "HYPESWEG"


