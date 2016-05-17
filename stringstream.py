from cStringIO import StringIO
import zlib
import io
import re

filename = 'assoc.q.lm.epacts.gz'
with open(filename, 'rb') as f:
	f.seek(316547)
	file = f.read()

def stream_gzip_decompress(stream):
	dec = zlib.decompressobj(32 + zlib.MAX_WBITS)
	rv = dec.decompress(stream)
	if rv:
		return rv

data = StringIO(stream_gzip_decompress(file))
data.next()
print data.readline()