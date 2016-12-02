#!/bin/bash

special_thousand=0
special_gencode=0

for dataset in individual_queries/individual 
#ten hundred
do


if [ $dataset == "ten" ]
then
lasti=9
else
	if [ $dataset == "individual_queries/individual" ]
	then
		lasti=50
	else
		lasti=2
	fi
fi


for i in `seq 0 $lasti`
do 

if [ $special_thousand == "1" ]
then
	query=thousand
else
	if [ $special_gencode == "1" ]
	then
		query=gencode.v25.transcripts
	else
		query=$dataset.$i
fi
fi

echo "query: $query"


# query prep

cat /mnt/assemblage/sbt/$query.fa \
  | grep -v "^>" \
  > temp/$query.query

# original

do_original=1

if [ $do_original == "1" ]
then

	bash drop_cache > /dev/null
	echo "original SK"

	/usr/bin/time ../src/bt query-original \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/compressedSBT/SBT_list.txt \
	  temp/$query.query \
	  dat/$query.unclust-rrr-original.dat \
	   > log/$query.unclust-original  2>&1

	tail -n 2 log/$query.unclust-original |head -n 1
fi

# original + opt

do_original_opt=1

if [ $do_original_opt == "1" ]
then

	bash drop_cache > /dev/null
	echo "original SK + opt"

	/usr/bin/time ../src/bt query-redux \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/compressedSBT/SBT_list.txt \
	  temp/$query.query \
	  dat/$query.unclust-rrr-original.dat \
	   > log/$query.unclust-original-opt  2>&1

	tail -n 2 log/$query.unclust-original-opt |head -n 1
fi

# clust alone

do_clust_alone=1

if [ $do_clust_alone == "1" ]
then
	bash drop_cache > /dev/null
	echo "clust"

	/usr/bin/time ../src/bt query-original \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/clusteringSBT.conversion/sbt-rrr.txt \
	  temp/$query.query \
	  dat/$query.clust-rrr-original.dat \
	  > log/$query.clust-original 2>&1

	tail -n 2 log/$query.clust-original |head -n 1

fi

# old clust+allsome

#echo tableId      = SBT+clust+allsome
#echo queryId      = clust-rrr-split
#echo queryProcess = query-redux
#echo topologyFile = /gpfs/cyberstar/pzm11/backup/sbt/clusteringSBT.compressed/sbt-rrr-split.txt

#/usr/bin/time ../src/bt query-redux \
#  --query-threshold 0.9 \
#  /mnt/assemblage/sbt/clusteringSBT.compressed/sbt-rrr-split.txt \
#  temp/$query.query \
#  clust-rrr-split.$query.dat > log/$query.split  2>&1

#tail -n 2 log/$query.split | head -n 1

# original SK + allsome

do_original_allsome=1

if [ $do_original_allsome == "1" ]
then

	bash drop_cache > /dev/null
	echo "original+allsome"

	/usr/bin/time ../src/bt query-redux \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/compressedSBT.allsome/sbt-rrr-allsome.txt \
	  temp/$query.query \
	  dat/$query.unclust-rrr-allsome.dat \
	   > log/$query.unclust-allsome 2>&1

	tail -n 2 log/$query.unclust-allsome | head -n 1
fi

# original SK + allsplit-split(2 files)

do_original_allsome_split=1

if [ $do_original_allsome_split == "1" ]
then

	bash drop_cache > /dev/null
	echo "original+allsome-split"

	/usr/bin/time ../src/bt query-redux \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/compressedSBT.reconstructed/sbt-rrr-split.txt \
	  temp/$query.query \
	  dat/$query.unclust-rrr-split.dat \
	   > log/$query.unclust-split 2>&1

	tail -n 2 log/$query.unclust-split | head -n 1
fi



# new clust+allsome

do_clust_allsome=1

if [ $do_clust_allsome == "1" ]
then

	bash drop_cache > /dev/null
	echo "clust+allsome (1-file)"

	/usr/bin/time ../src/bt query-redux \
	  --query-threshold 0.9 \
	  /mnt/assemblage/sbt/clusteringSBT.compressed/sbt-rrr-allsome.txt \
	  temp/$query.query \
	  dat/$query.clust-rrr-allsome.dat \
	   > log/$query.clust-allsome 2>&1

	tail -n 2 log/$query.clust-allsome | head -n 1
fi




if [ $special_thousand == "1" ] || [ $special_gencode == "1" ]
then
exit 0
fi

done
done


