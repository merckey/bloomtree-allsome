import numpy
from collections import defaultdict
mode="mean"
#mode="all"
methods =["unclust-original","unclust-allsome","clust-original","clust-allsome"]
print "      "+'   '.join(methods)
times=defaultdict(lambda : defaultdict(list))
for d in ["individual_queries/individual","ten","hundred", "thousand"]:
	r = 10 if (d=="ten") else 3
	if d=="thousand": r = 1
	if d=="individual_queries/individual": r = 50
	for i in xrange(0,r):
		dataset=d+"."+str(i)
		if d=="thousand": dataset=d
		if mode=="all":
			print dataset+"  ",
		for method in methods:
			with open("log/"+dataset+"."+method) as f:
				time= f.readlines()[-2].split()[2].strip('elapsed')[:-3].replace(':','m')
				t = time.split('m')
				seconds = int(t[0])*60 +int(t[1])
				times[method][d].append(seconds)
				if mode=="all":
					print time+"s                ",
					
		if mode=="all":
			print 

	if mode == "mean":
		print "%11s" % (d if d != "individual_queries/individual" else "individual"), 
	for method in methods:
		if mode=="mean":
			#print method,d,times[method][d]
			mean, std= numpy.mean(times[method][d]), numpy.std(times[method][d])
			if d != "thousand":			
				stdpart =  "(%0.1fs)" % std
			else:
				stdpart = ""
			print   "%2dm%2d" % ( mean / 60,mean%60) + "s " + stdpart + "        ",
	if mode=="mean":
		print
