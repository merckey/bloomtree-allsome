from sys import argv
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

def nodes_visited_per_query(filename):
	query_num = 0
	nodes_visited_list = []
	with open(filename) as f:
		for line in f.readlines():
			line = line.strip()
			if line.startswith('#'):
				query_num += 1
				nodes_visited = line.split()[1]
				nodes_visited_list.append(float(nodes_visited))
	nodes_visited_num = sum(nodes_visited_list)
	result = float(nodes_visited_num)/float(query_num)
	return float(result)

query_result_dir = '/storage/home/cxs1031/standard/sbt_private/query_stat/'

dataset = []

cluster_ten_data = []

for i in range(10):
	query_out_file = query_result_dir + 'clustered.ten.' + str(i) + '.seq.query.out'
	cluster_ten_data.append(nodes_visited_per_query(query_out_file))

cluster_hun_data = []

for i in range(3):
	query_out_file = query_result_dir + 'clustered.hundred.' + str(i) + '.seq.query.out'
	cluster_hun_data.append(nodes_visited_per_query(query_out_file))

cluster_thu_data = []

query_out_file = query_result_dir + 'clustered.thousand.seq.query.out'
cluster_thu_data.append(nodes_visited_per_query(query_out_file))

solomon_ten_data = []

for i in range(10):
	query_out_file = query_result_dir + 'solomon.ten.' + str(i) + '.seq.query.out'
	solomon_ten_data.append(nodes_visited_per_query(query_out_file))

solomon_hun_data = []

for i in range(3):
	query_out_file = query_result_dir + 'solomon.hundred.' + str(i) + '.seq.query.out'
	solomon_hun_data.append(nodes_visited_per_query(query_out_file))

solomon_thu_data = []

query_out_file = query_result_dir + 'solomon.thousand.seq.query.out'
solomon_thu_data.append(nodes_visited_per_query(query_out_file))

dataset = [cluster_ten_data, solomon_ten_data, cluster_hun_data, solomon_hun_data, cluster_thu_data, solomon_thu_data]

fig = plt.figure(1)

ax = fig.add_subplot(111)

bp = ax.boxplot(dataset)

ax.set_xticklabels(['C_10', 'S_10', 'C_100', 'S_100', 'C_1000', 'S_1000'])

plt.ylabel('Nulber of Nodes Visited per Query')

fig.savefig('fig.png', bbox_inches='tight')
