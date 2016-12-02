grep "^*" clust-rrr-original.ten.6.dat |sort > clust-rrr-original.ten.6.sequences.dat
grep "^*" clust-rrr-split.ten.6.dat |sort > clust-rrr-split.ten.6.sequences.dat
grep "^*" clust-rrr-allsome.ten.6.dat |sort > clust-rrr-allsome.ten.6.sequences.dat
diff clust-rrr-original.ten.6.sequences.dat clust-rrr-split.ten.6.sequences.dat
diff clust-rrr-split.ten.6.sequences.dat clust-rrr-allsome.ten.6.sequences.dat
