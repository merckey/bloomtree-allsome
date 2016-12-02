# BFCluster


BFCluster is a toolkit of SBT-AS. It contains the following functional modules: 


- `sbuild`: build sequence bloom tree by hierarchical clustering leaf bloom filters
- `squery`: querying bloom filter converted from all sequences in query file, instead of querying kmers in query file
- `pcompress`: parallel compress sequence bloom tree to RRR or Roaring compressed sequence bloom tree (This module requires Intel Threading Building Blocks installed) 


## Detailed Usage


1.   `sbuild -f hashfile -p leaf_bloom_filter_directory -l bloom_filter_name.list  -o output_directory -b bloom_filter_tree_topology_file`


- `-f` hash file used to generate leaf bloom filters
- `-p` directory of leaf bloom filters
- `-l` bloom filter name list file, each line in the file is a leaf bloom filter name, prefix excluded
    For instance, if one of your leaf bloom filters is "/example/leaves/SAMPLE.bf.bv", then the corresponding line in list is "SAMPLE"
- `-o` output directory of sequence bloom tree
- `-b` filename of tree topology




2. `squery -b tree_topology_file -f hashfile -q query.fasta -o output_directory -t threshold`


- `-b` sequence bloom tree topology file
- `-f` hash file used to generate bloom filters
- `-q` query sequences in FASTA file format
- `-o` query output directory
- `-t` query threshold




3. `pcompress -l bloom_filter_file_list -c compress_type -t thread`
- `-l` bloom filter file list, each line in the file is a bloom filter file. A valid sequence bloom tree topology file is also OK.
- `-c` compress_type, should be either "rrr" or "roar", defaultly use "rrr" if not indicated.
- `-t` number of threads used


