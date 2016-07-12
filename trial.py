import os

print os.getcwd()
file = open(os.path.join(os.getcwd(), "lz_assoc_viewer/teach.txt"), 'r')

print file.read()