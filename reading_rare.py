import re
import gzip
with open("cleanfile.gz", 'rb') as f:
	data = f.read()
	outfile_clean = open("clean_out.txt", 'w')
	outfile_clean.write(data)

#with open("plink.qassoc.gz") as f:
#	data = f.read()

