from sys import argv

filename = argv[1]

with open(filename) as f:
	for line in f.readlines():
		line = line.strip()
		print line.replace('bv', 'bv.roar')