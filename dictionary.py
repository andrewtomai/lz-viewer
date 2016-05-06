import gzip
from lz import *

header = check_header('assoc.q.lm.epacts.gz')
print get_column(header, "position")
with gzip.open('assoc.q.lm.epacts.gz', 'rb') as f:
	next(f)
	for line in f:
		words = line.split()

position = find_min_pval('assoc.q.lm.epacts.gz', header)
print position