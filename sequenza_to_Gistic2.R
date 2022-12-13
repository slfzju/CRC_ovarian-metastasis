#!/bin/bash
##Num markers = N.BAF; Seg.CN = depth.ratio (but we will log transform this)
dirname="./segments";
sam=all54;
ls *.tsv >config1;
cat config1 | while read id;
     do
   k=${id};
  echo $k
cat $k | sed "1d" >> 202202_${sam}_allseg.tsv;
done

cat 202202_${sam}_allseg.tsv | awk 'BEGIN {print "Sample\tChromosome\tStart Position\tEnd Position\tNum markers\tSeg.CN"}{print $0}'>202202_${sam}_allseg2.tsv;
mkdir 202202_${sam}_allseg2;
mv *.tsv 202202_${sam}_allseg2;
rm config1

(base):/opt/GISTIC# ./gistic2 -b /mydata/2022_${sam}_seg_tsv -seg /mydata/202202_${sam}_allseg2/202202_dqf_allseg2.tsv -refgene ./refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -conf 0.9
