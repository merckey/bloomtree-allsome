from sys import argv
import os
import matplotlib
# use agg to output plot to png
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# this script take three paramters:
# 	query_result_of_clustered_tree  	query output file on clustered tree
#	query_result_of_solomon_tree 		query output file on solomon tree
#	picture_name.png					output picture path
if (len(argv) < 3):
	print 'paramters: query_result_of_clustered_tree query_result_of_solomon_tree picture_name.png'

alpha_val = 0.3
def read_leaf_hits(filename):
	'''
		read query output file
		return two lists:
			1. number of leaves hits by each query
			2. number of nodes visited by corresponding each query
	'''
	leaf_hits_list = []
	nodes_visited_list = []
	with open(filename) as f:
		for line in f.readlines():
			line = line.strip()
			if line.startswith('*'):
				leaf_hits = float(line.split()[-1])
				leaf_hits_list.append(leaf_hits)
			if line.startswith('#'):
				nodes_visited = float(line.split()[1])
				nodes_visited_list.append(nodes_visited)
	assert len(leaf_hits_list) == len(nodes_visited_list)
	return leaf_hits_list, nodes_visited_list

(cluster_leaf_hits, cluster_nodes_visited) = read_leaf_hits(argv[3])
(solomon_leaf_hits, solomon_nodes_visited) = read_leaf_hits(argv[1])
(allsome_leaf_hits, allsome_nodes_visited) = read_leaf_hits(argv[2])
(allsome_clustered_leaf_hits, allsome_clustered_nodes_visited) = read_leaf_hits(argv[4])

def calculate_average_distance(l1, l2):
	assert len(l1) == len(l2)
	l1_visited = 0
	l2_visited = 0
	for c in l1:
		l1_visited += c
	for c in l2:
		l2_visited += c
	l1_avg = l1_visited / len(l1)
	l2_avg = l2_visited / len(l2)

	return  str(l1_avg) + ',' + str(l2_avg) + ',' + str((l1_avg-l2_avg)/l1_avg)

visited_list = [solomon_nodes_visited, allsome_nodes_visited, cluster_nodes_visited, allsome_clustered_nodes_visited]
visited_label_list = ['sbt-sk', 'allsome', 'clustering', 'sbt-also']

for i in range(4):
	for j in range(i+1,4,1):
		print 'average distance bwtween ' + visited_label_list[i] + ' and ' + visited_label_list[j] + ': ', calculate_average_distance(visited_list[i], visited_list[j])

def calculate_average_large_distance(x1, y1, x2, y2):
	assert len(y1) == len(y2)
	l1_visited = 0
	l2_visited = 0
	count = 0
	for i in range(len(x1)):
		if x1[i] > 800:
			l1_visited += y1[i]
			l2_visited += y2[i]
			count += 1

	l1_avg = l1_visited / count
	l2_avg = l2_visited / count

	return  str(l1_avg) + ',' + str(l2_avg) + ',' + str((l1_avg-l2_avg)/l1_avg)

print 'average large distance bwtween solomon and allsome: ', calculate_average_large_distance(solomon_leaf_hits, solomon_nodes_visited, allsome_leaf_hits, allsome_nodes_visited)

def visits_less_than_hits(x, y):
	assert len(x) == len(y)
	num = 0
	for i in range(len(x)):
		if x[i] >= y[i]:
			num += 1
	return str(num) + ',' + str(len(x)) + ',' + str(num / len(x))

print 'allsome, % #visits less than hits: ', visits_less_than_hits(allsome_leaf_hits, allsome_nodes_visited)
print 'SBT-ALSO, % #visits less than hits: ', visits_less_than_hits(allsome_clustered_leaf_hits, allsome_clustered_nodes_visited)

ax = plt.gca()
# plot points
sk_dot, = ax.plot(solomon_leaf_hits, solomon_nodes_visited, linestyle='None', color='green', marker='o', alpha=alpha_val)#, label='SBT-SK(base)')
as_dot, = ax.plot(allsome_leaf_hits, allsome_nodes_visited, linestyle='None', color='blue', marker='o', alpha=alpha_val)#, label='+ AllSome nodes')
cl_dot, = ax.plot(cluster_leaf_hits, cluster_nodes_visited, linestyle='None', color='red', marker='s', alpha=alpha_val)#, label='+ Clustering')
psu_dot, = ax.plot(allsome_clustered_leaf_hits, allsome_clustered_nodes_visited, linestyle='None', color='orange', marker='s', alpha=alpha_val)# label='+ AllSome & Clustering')

from scipy.optimize import curve_fit
import numpy as np
from numpy.polynomial import polynomial as P


#x = np.array(solomon_leaf_hits, dtype=float)
#y = np.array(solomon_nodes_visited, dtype=float)

def poly_curve_fit(leaf_hits, nodes_visited, color_sign):
	d = {}
	for i in range(len(leaf_hits)):
		x_element = leaf_hits[i]
		if x_element in d:
			d[x_element] = (d[x_element] + nodes_visited[i])/2
		else:
			d[x_element] = nodes_visited[i]



	#x_tmp = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600]

	x_tmp = [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1600]

	x = []
	y = []

	for x_element in x_tmp:
		total = 0
		count = 0
		for i in range(-20, 20, 1):
			k = x_element + i
			if k in d:
				total += d[k]
				count += 1
		if count == 0:
			continue
		x.append(x_element)
		y.append(total/count)

	last_x = sorted(d)[-1]
	last_y = d[last_x]

	x.append(last_x)
	y.append(last_y)

	z1 = np.polyfit(x, y, 3)
	p1 = np.poly1d(z1)

	yvals=p1(x)

	z1 = np.polyfit(x, y, 2)
	p1 = np.poly1d(z1)
	yvals=p1(x)
	#plot1=plt.plot(x, y, '^', color = color_sign)
	plot2, =ax.plot(x, yvals, color = color_sign, lw = 4)
	return plot2

sk_line = poly_curve_fit(solomon_leaf_hits, solomon_nodes_visited, 'green');
as_line = poly_curve_fit(allsome_leaf_hits, allsome_nodes_visited, 'blue');
cl_line = poly_curve_fit(cluster_leaf_hits, cluster_nodes_visited, 'red');
psu_line = poly_curve_fit(allsome_clustered_leaf_hits, allsome_clustered_nodes_visited, 'orange');

#print sk_dot, as_dot, cl_dot, psu_dot
#print sk_line, as_line, cl_line, psu_line

ax.legend([(sk_dot, sk_line), (as_dot, as_line), (cl_dot, cl_line), (psu_dot, psu_line)], ['SBT-SK(base)', '+ AllSome nodes', '+ Clustering', '+ AllSome & Clustering'], loc='best')
#ax.legend([(as_dot, as_line)], ['+ AllSome nodes'])
#ax.legend([(cl_dot, cl_line)], ['+ Clustering'])
#ax.legend([(psu_dot, psu_line)], ['+ AllSome & Clustering'])

# set legend
#plt.legend(loc='best')

#plt.plot([0, 2652], [1600, 2652], '--')

# plot line
plt.axhline(y=2652, linestyle='--', color='black')

# set x and y label
plt.xlabel('Number of Leaf Hits')
plt.ylabel('Nulber of Nodes Examined')

# save figure into file
plt.savefig(argv[5])
