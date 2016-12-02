from sys import argv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

def read_level_occupancy(filename):
	level_list = []
	occupancy_list = []
	with open(filename) as f:
		for line in f.readlines():
			line = line.strip()
			columns = line.split()
			occupancy = float(columns[2])
			level = columns[4].count('*')
			if(level > 30):
				continue
			level_list.append(level)
			occupancy_list.append(occupancy)

	return level_list, occupancy_list

(cluster_level_list, cluster_occupancy_list) = read_level_occupancy(argv[1])
(solomon_level_list, solomon_occupancy_list) = read_level_occupancy(argv[2])

plt.plot(cluster_level_list, cluster_occupancy_list, linestyle='None', color='red', marker='*', label='Clustered Tree')
plt.plot(solomon_level_list, solomon_occupancy_list, linestyle='None', color='green', marker='.', label='Solomon Tree')

plt.xlabel('Node Depth')
plt.ylabel('Bloom Filter Occupancy')

plt.legend(loc='best')

plt.savefig(argv[3])