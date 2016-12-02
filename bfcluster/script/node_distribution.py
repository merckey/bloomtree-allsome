from sys import argv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

def read_node_dist(filename):
	level_num = {}
	with open(filename) as f:
		for line in f.readlines():
			line = line.strip()
			columns = line.split()
			occupancy = float(columns[2])
			level = columns[4].count('*')
			if level in level_num:
				level_num[level] += 1
			else:
				level_num[level] = 1
	level_list = []
	node_num_list = []
	
	for level in sorted(level_num):
		level_list.append(level)
		node_num_list.append(level_num[level])

	return level_list, node_num_list

(cluster_level_list, cluster_node_list) = read_node_dist(argv[1])
(solomon_level_list, solomon_node_list) = read_node_dist(argv[2])

plt.plot(cluster_level_list, cluster_node_list, color='red', marker='.', label='Clustered Tree')
plt.plot(solomon_level_list, solomon_node_list, color='green', marker='.', label='Solomon Tree')

plt.xlabel('Tree Depth')
plt.ylabel('Number of Nodes')

plt.legend(loc='best')

plt.savefig(argv[3])