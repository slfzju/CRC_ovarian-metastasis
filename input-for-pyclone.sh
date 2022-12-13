#!bash/pyclone-input-tsv
ls /vep-maf/*.maf> config1;
ls /segment/*.cr.igv.seg >config2;
cat config1|while read id;
do 
i=${id};
echo $i 
cat $i | sed '1,2d' | awk -F '\t' '{print $5"\t"$6"\t"$7"\t"$5":"$6":"$1"\t"$41"\t"$42"\t"2"\t"0"\t"2"\t"}' ${i%%.maf}.tmp.tsv;
done

cat config1 | while read id;
i=${id};
echo $i
do
cat $i | sed '1,2d' | awk -F '\t' 'BEGIN{print "Chromosome\tStart_Position\tEnd_Position\tmutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn"}{print $5"\t"$6"\t"$7"\t"$5":"$6":"$1"\t"$41"\t"$42"\t"2"\t"0"\t"2"\t"}' ${i%%.maf}.tmp.tsv;
done

cat config2 | while read id;
j=${id};
echo $j
do 
cat $j | sed '1d' | awk 'BEGIN{OFS="\t"}{print $0"\t"int((2^$6)*2+0.5)}'| awk 'BEGIN{OFS="\t"}{if ($7!=0)print $0}' | cut -f 2-6  >/Users/sunlf/202006mac-yuli-GatkCNV/yl-segment/${j%%.crv*}.bed;
done

ls *.tmp.tsv >gz1;
ls *.bed>gz2;
paste gz1 gz2>config3;
cat config2 | while read id;	
do
      arr=${id};
      fq1=${arr[0]};
      fq2=${arr[1]};
echo $fq1 $fq2
do 
bedtools window -a $fq1 -b fq2 | cut -f 4-8,15 | awk 'BEGIN{OFS="\t";print "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn"}{print $0}' >${fq2%%.bed}.tsv;
done

#!/bin/bash/pyclone
cd /mydata/pyclone_tsv

ls *.tsv >config6;
cat config6 | while read id;
     do
   k=${id};
  echo $k
PyClone run_analysis_pipeline --prior major_copy_number --in_files $k --working_dir ./${k:0:9}-pyclone-analysis --init_method connected --max_clusters 80 --num_iters 1000 1>${k:0:9}py.log 2>&1;
done

mv ./${k:0:9}-pyclone-analysis/tables/loci.tsv ./${k:0:9}-pyclone-analysis/tables/${k:0:9}_loci.tsv;
mv ./${k:0:9}-pyclone-analysis/tables/cluster.tsv ./${k:0:9}-pyclone-analysis/tables/${k:0:9}_cluster.tsv;

