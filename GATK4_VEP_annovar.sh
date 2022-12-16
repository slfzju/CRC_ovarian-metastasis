#Trimmomatic-GATK4-mutect2-vep/annovar
#Vi Trimmomatic.sh
#!/bin/bash
cd /fastq_rawdata/
ls|grep _R1.fastq.gz>gz1;
ls|grep _R2.fastq.gz>gz2;
paste gz1 gz2>config;
cat config |while read id;
do
      arr=(${id});
      fq1=${arr[0]};
      fq2=${arr[1]};
echo $fq1 $fq2
/anaconda3/envs/bin/java -jar ~/bin/trimmomatic.jar PE -threads 16 -phred33 $fq1 $fq2 out.$fq1 unpair.$fq1 out.$fq2 unpair.$fq2 ILLUMINACLIP:/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:50;done

#!/bin/bash
cd /fastq_rawdata/
ls|grep out_*fastq.gz>gz1;
ls|grep out_*R2.fastq.gz>gz2;
paste gz1 gz2>config2;
cat config2 |while read id;
do
      arr=(${id});
      fq1=${arr[0]};
      fq2=${arr[1]};
echo $fq1 $fq2
bwa mem -t 20 -M -Y -R '@RG\tID:$fq1WES\tSM:$fq1\tLB:WES\tPL:Illumina' /hg38ref/hg38.fa $fq1 $fq2 -o $fq1.sam.gz;done

##sam-to-bam
#!/bin/bash
ls|grep *sam.gz>config3;
cat config3 |while read id;
do
   i=${id}
   samtools view -Sb -t 10 ${i} >/media/sunlf/data/${id}%.bam.gz;
done

#!/bin/bash
ls *.bam.gz>config4;
cat config4 | while read id;
     do
    i=${id};
     echo $i
        samtools sort -@ 10 -l 9 -o ${i}.sorted.bam.gz ${i};done
        
#!/bin/bash/sorted
ls *.sorted.bam.gz>config5;
cat config5 | while read id;
     do
    i=${id};
        echo ${i}
java -jar /anaconda3/envs/share/gatk4-4.0.1.1-0/gatk4.jar ValidateSamFile -I $i --MODE SUMMARY | samtools quickcheck -qv $i > ${i:0:9}_bams.fofn \&& echo 'all ok'\ | samtools index $i;
done

#!/bin/bash
ls *.sorted.bam.gz>config5 ;
cat config5 | while read id;
     do
    i=${id};
         echo ${i}
java -jar /anaconda3/envs/share/gatk4-4.0.1.1-0/gatk4.jar MarkDuplicates --ASSUME_SORT_ORDER coordinate -I $i -O $i.marked.bam.gz -M $i.metrics;
done

#!/bin/bash
ls *.marked.bam.gz>config6 ;
cat config6 | while read id;
     do
    i=${id};
         echo ${i}
java -jar -Xmx10g /anaconda3/envs/share/gatk4-4.0.1.1-0/gatk4.jar FixMateInformation --SORT_ORDER coordinate -I $i -O ${i%%marked*}fixed.bam.gz;
done

#!/bin/bash
ls *.fixed.bam.gz>config7 ;
cat config7 | while read id;
     do
    j=${id};
         echo ${j}
samtools index -@ 16 $j;
done

#!/bin/bash
gatk=/anaconda3/envs/share/gatk4-4.0.1.1-0/gatk4.jar;
ref=/hg38ref/hg38.fa;
dbsnp=/hg38ref/dbsnp_146.hg38.vcf;
gold=/hg38ref/Mills_and_1000G_gold_standard.indels.hg38.vcf;
ls *.fixed.bam.gz>config7;
cat config7 | while read id;
do
        j=${id};
        echo ${j}
java -jar -Xmx20g $gatk BaseRecalibrator -R $ref -I $j --known-sites $dbsnp --known-sites $gold -O $j-recal.table -OBI true;done

ref=/hg38ref/hg38.fa;
ls *.fixed.bam.gz>gz7;
ls *-recal.table >gz8
paste gz7 gz8>config8
cat config8 | while read id;
do
  arr=(${id});
      fq1=${arr[0]};
      fq2=${arr[1]};
echo $fq1 $fq2
java -jar -Xmx20g $gatk ApplyBQSR -R $ref -I $fq1 -bqsr $fq2 -O ${fq1%%.fixed*}.bqsr.bam.gz --read-index coordinate -OBI true;done

ls *.bqsr.bam.gz>config9;
cat config9 | while read id;
     do
     j=${id};
         echo ${j}
samtools index -@ 16 $j | $gatk ValidateSamFile $j;
done

dbsnp=/hg38ref/dbsnp_146.hg38.vcf;
ls *.bqsr.bam.gz>config9;
cat config9 | while read id;
do
        j=${id};
        echo ${j}
java -jar -Xmx20g $gatk HaplotypeCaller -R $ref -I $j -dbsnp $dbsnp -O ${j%%.bqsr*}.raw.vcf;done

#!/bin/bash
gatk=/anaconda3/envs/share/gatk4-4.0.1.1-0/gatk4.jar;
ref=/hg38ref/hg38.fa;
dbsnp=/hg38ref/dbsnp_146.hg38.vcf;
omini=/hg38ref1000G_omni2.5.hg38.vcf;
snps=/hg38ref/1000G_phase1.snps.high_confidence.hg38.vcf;
gold=/hg38ref/Mills_and_1000G_gold_standard.indels.hg38.vcf;
hapmap=/hg38ref/2020hg38/hapmap_3.3.hg38.vcf;
ls *.raw.vcf>config10;
cat config10 | while read id;
do
        j=${id};
        echo ${j}
java -jar -Xmx16g $gatk VariantRecalibrator -R $ref -V $j --resource omni,known=false,training=true,truth=false,prior=12.0:$omini --resource 1000G,known=false,training=true,truth=false,prior=10.0:$snps --resource hapmap,known=false,training=true,truth=true,prior=15.0:$hapmap --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$dbsnp -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 --tranches-file ${j%%.raw*}.snps.tranches --rscript-file ${j%%.raw*}.snps.plots.R -O ${j%%.raw*}.snps.recal --max-gaussians 4;done

#!/bin/bash/snpvqsr
ls *.raw.vcf>gz9;
ls *.snps.tranches >gz10;
ls *snps.recal>gz11;
paste gz9 gz10 gz11>config11;
cat config11 | while read id;
do
      arr=(${id});
      fq1=${arr[0]};
      fq2=${arr[1]};
      fq3=${arr[2]};
echo $fq1 $fq2 $fq3
java -jar -Xmx20g $gatk ApplyVQSR -R $ref -V $fq1 --truth-sensitivity-filter-level 90.0 --tranches-file $fq2 --recal-file $fq3 -mode SNP -O ${fq1%%.raw*}.snps.VQSR.vcf;done

ls *.snps.VQSR.vcf >config12;
cat config12 | while read id;
do
        j=${id};
        echo ${j}
java -jar -Xmx20g $gatk VariantRecalibrator -R $ref -V $j --resource mills,known=true,training=true,truth=true,prior=12.0:$gold -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL --max-gaussians 6 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 --rscript-file ${j%%.snps*}.HC.snps.indels.plots.R --tranches-file ${j%%.snps*}.HC.snps.indels.tranches -O ${j%%snps*}.HC.snps.indels.recal;done

#!indelvqsr
ls *.snps.VQSR.vcf >gz12;
ls *.HC.snps.indels.tranches >gz13;
ls *.HC.snps.indels.recal >gz14;
paste gz12 gz13 gz14>config13;
cat config13 | while read id;
do
      arr=(${id});
      fq1=${arr[0]};
      fq2=${arr[1]};
      fq3=${arr[2]};
echo $fq1 $fq2 $fq3
java -jar -Xmx20g $gatk ApplyVQSR -R $ref -V $fq1 --truth-sensitivity-filter-level 90.0 --tranches-file $fq2 --recal-file $fq3 -mode INDEL -O ${fq1%%.snps*}.HQ.vcf;done

#GATK_mutect2
cd /fastq_rawdata/;
ls *HQ.vcf >config14;
cat config14 | while read id;
do
        i=${id};
         echo $i
gatk Mutect2 \
   -R $ref \
   -I $i \
   -tumor ${i%%.HQ*}
   -I normal.bam \
   -normal normal_sample_name \
   --germline-resource af-only-gnomad.vcf.gz \
   --af-of-alleles-not-in-resource 0.00003125 \
   --panel-of-normals pon.vcf.gz \
   -O ${i%%.HQ*}.somatic.vcf.gz;done
   
#VCF-to-VEP
#!/bin/bash
cd /ensembl-vep/;
ls /fastq_rawdata/*.somatic.vcf.gz >config15;
cat config15 | while read id;
do
        i=${id};
         echo $i
./vep --cache --offline --format vcf --vcf --force_overwrite --dir_cache /ensembl-vep/vep_data/cache/ --input_file $i --output_file /ensembl-vep/vep_data/output/${i%%.somatic*}-vep.vcf --dir_plugins /ensembl-vep/vep_data/plugins --plugin dbNSFP, /ensembl-vep/vep_data/plugins/dbNSFP4.0b2a.txt.gz,ALL;done

#vep2maf
#!/bin/bash/vep2maf
ls /ensembl-vep/vep_data/output/*vep.vcf >config1;
cat config16 | while read id;
   do
   i=${id};
   echo $i
perl vcf2maf.pl --input-vcf ${i} --output-maf ${i%%.vcf}.maf --filter-vcf 0 --ref-fasta $ref --inhibit-vep --tumor-id ${i%%-vep.vcf} --vcf-tumor-id ${i%%-vep.vcf} --ncbi-build GRCh38 --species homo_sapiens;
done

#Annovar
#!/bin/bash
dir=/annovar/converted-anno
find ./fastq_rawdata/*.somatic.vcf.gz >config17;
cat config17|while read id;
do
i=$id;
echo $id
perl convert2annovar.pl -format vcf4 -allsample $id -outfile anno;
mkdir ./converted-anno;
mv *.avinput ./converted-anno;
done

datadir=/hg38/hg38exome/;
find ./converted-anno -name "*.avinput" >config18;
cat config18|while read id;
do
i=$id;
echo $id
perl table_annovar.pl $id $datadir --buildver hg38 --outfile ${id%%.avinput} --remove --protocol refGene,esp6500siv2_all,Exac03,clinvar_20190305 --operation g,f,f,f --nastring . --polish --thread 50;
mv *.csv ./anno-to-maf;
done
