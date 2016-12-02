#!/bin/sh
rm individual_queries/*.fa
seqtk sample gencode.v25.transcripts.fa 50 > individual_queries/50_queries.fa
cd individual_queries
split -l 2 50_queries.fa
num=0
for file in x*; do
       mv "$file" "individual.$(printf "%u" $num).fa"
       let num=$num+1
done
