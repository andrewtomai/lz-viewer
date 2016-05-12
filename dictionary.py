import os
import gzip
import time
f = open("test_file.txt")
f.seek(15)
for line in f:
	print line

f =  gzip.open("assoc.q.lm.epacts.gz")
print f.length()

middle = 323503/2
f.seek(middle)
f.next()
line = f.readline()
data = line.split()
position = long(data[1])

byte_location = middle
max_byte = 323503
min_byte = 0

print position
prev_position = -5
while position != 199411:
	if prev_position == position:
		break

	elif( position < 199411 ):
		min_byte = byte_location
		byte_location = (max_byte - byte_location)/2 + byte_location
		f.seek(byte_location)
		f.next()
		line = f.readline()
		data = line.split()
		prev_position = position
		position = long(data[1])
		print position
		print position < 73078586
		print byte_location
		time.sleep(1)


	elif( position > 199411 ):
		max_byte = byte_location
		byte_location = (byte_location - min_byte)/2 + min_byte
		f.seek(byte_location)
		f.next()
		line = f.readline()
		data = line.split()
		prev_position = position
		position = long(data[1])

		print position
		time.sleep(1)

print f.readline()


f2 = gzip.open('assoc.q.lm.epacts.gz.tbi')
for line in f2:
	print f2.readline()
