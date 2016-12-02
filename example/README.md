# Tutorial, Creating an Allsome Bloomtree

Initialize the hash function.

```bash  
bt hashes --k 20 example.hashfile
```

Estimate the best bloom filter size ... note that kmerSize 20 is hardwired in
get_bfsize.sh.

```bash  
ls experiment*.fa \
  | while read f ; do
	  echo "=== gzipping ${f} ==="
	  cat ${f} | gzip > ${f}.gz
	  done
ls experiment*.fa.gz > experimentfiles
../src/get_bfsize.sh experimentfiles
```

Output from get_bfsize.sh:  
```bash  
    Unique canonical kmer counts by number of occurences:
    1 4906564
    2 963379
    3 499717
```

So we want to use a bf_size of about 500,000.

Count kmers in each experiment.

```bash  
bf_size=500000
ls experiment*.fa \
  | while read f ; do
	  bv=`echo ${f} | sed "s/\.fa/.bf.bv/"`
	  echo "=== converting ${f} to ${bv} ==="
	  bt count --cutoff 3 example.hashfile ${bf_size} ${f} ${bv}
	  done
```

Build the SBT.

```bash  
ls experiment*.bf.bv | sed "s/\.bf\.bv//" > bitvectornames
../bfcluster/sbuild -f example.hashfile -p . -l bitvectornames -o . -b sbt.txt
```

Split the SBT nodes to create an SBT-AS.

```bash  
bt split sbt.txt sbt-as.txt
```

Compress the SBT-AS bit vectors, creating the node files for a compressed ABT-AS.

```bash  
cat sbt.txt \
  | tr -d "*" \
  | sed "s/,.*//" \
  | sed "s/\.bf\.bv//" \
  | while read node ; do
	  echo "=== compressing split ${node} to allsome ==="
	  bt compress-rrr-double example.hashfile \
		${node}.bf-all.bv \
		${node}.bf-some.bv \
		${node}.bf-allsome.bv
	  done
```

Create the topology file for the compressed ABT-AS.

```bash  
cat sbt.txt \
  | sed "s/\.bf\.bv/.bf-allsome.bv.rrr/" \
  > sbt-rrr-allsome.txt
```

Run a batch of queries.

```bash  
bt query -t 0.5 sbt-rrr-allsome.txt queries queryresults
```
