mkdir 1_RawData 2_CleanData 3_Mapping 4_SNP_InDel 5_Annotation 6_Disease

/home/genome/perl/Assembly_stat_evaluation.pl ref.fna ref.summary
###ref preparation
mkdir 0_Reference
 #构建 Index，bwa 的版本号为0.7.12-r1039
bwa index Reference_Sequence.fa
#创建 dict 文件
java -jar /home/genome/soft/picard-tools-1.124/picard.jar CreateSequenceDictionary R=ref.fasta O=ref.fasta.dict
#创建 fai 文件
samtools faidx Reference_Sequence.fa

mv *fastq.gz 1_RawData
cd 1_RawData
gzip -d *fastq.gz
mv *R1*fastq ${name}.R1.fastq
mv *R2*fastq ${name}.R2.fastq
mkdir fastQC
fastqc ${name}.R1.fastq -o fastQC -t 16 &
fastqc ${name}.R2.fastq -o fastQC -t 16 &

perl /home/genome/perl/Q20_Q30_GC.pl 33 {name}.summary ${name}.R1.fastq ${name}.R2.fastq &
###AdapterRemoval and quality trimming(2_RawData)
AdapterRemoval  --threads  24 --file1 $raw[0] --file2 $raw[1] --output1 ../2_CleanData/$name.R1.fq --output2 ../2_CleanData/$name.R2.fq
cd ../2_CleanData
mv ${name}.R1.fq  ${name}_AD_R1.fq
mv ${name}.R2.fq  ${name}_AD_R2.fq
perl /home/genome/perl/QualityTrim.pl ${name}_AD_R1.fq ${name}_AD_R2.fq ${name}_HQ
perl /storage/RawData/Jiama/blii/bin/Q20_Q30_GC.pl 33 ${name}.summary *_HQ_clean_R1.fq *_HQ_clean_R2.fq &
###Align files(3_Mapping)
cd ../3_Mapping
bwa aln -t 18 ../0_Reference/ref.fasta ../2_CleanData/SL.HQ_clean_R1.fq -f SL_HQ_R1.fq.sai &
bwa aln -t 18 ../0_Reference/ref.fasta ../2_CleanData/SL.HQ_clean_R2.fq -f SL_HQ_R2.fq.sai &

bwa sampe -a 400 -r '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' ../0_Reference/ref.fna SL_HQ_R1.fq.sai SL_HQ_R2.fq.sai ../2_CleanData/SL.HQ_clean_R1.fq ../2_CleanData/SL.HQ_clean_R2.fq > SL.sam &

##sort sam file
mkdir -p ../temp/picardSortSam
java -jar /home/genome/soft/picard-tools-1.124/picard.jar SortSam INPUT=${name}.sam OUTPUT=${name}.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=../temp/picardSortSam

##mark duplicates
mkdir ../temp/picardMarkDuplicates
java -jar /home/genome/soft/picard-tools-1.124/picard.jar MarkDuplicates INPUT=SL.sorted.bam OUTPUT=SL.dedup_reads.bam METRICS_FILE=SL.dedup.metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=../temp/picardMarkDuplicates

##get possible indels
$gatk -T RealignerTargetCreator -R ../0_Reference/ref.fasta -I SL.dedup_reads.bam -o SL.intervals &
##realign
nohup $gatk -T IndelRealigner -targetIntervals SL.intervals -R ../0_Reference/ref.fasta -I SL.dedup_reads.bam -o SL.realigned.bam &

##DepthOfCoverage
$gatk -T DepthOfCoverage -R ../0_reference/ref.fasta -I SL.realigned.bam -o SL.DepthOfCoverage -nt 4 --omitIntervalStatistics &
	
##统计maping情况
perl /home/genome/Shuyun/bin/stat_bwa_aln_v2.pl ../3_Mapping/SL.sorted.bam SL.sorted.bam.meminfo &
perl /home/genome/Shuyun/bin/depth_v2.pl -l 119667750 SL.sorted.bam SL.plot > SL.depthinfo &
samtools flagstat SL.realigned.bam > SL.realigned.bam.meminfo
###SNP InDel calling(4_SNP_InDel)
cd ../4_SNP_InDel	
##get raw variants
nohup $gatk -T HaplotypeCaller -R ../0_Reference/ref.fasta -I ../3_Mapping/SL.realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o SL.raw_variants.vcf &
	
##get raw snps
nohup $gatk -T SelectVariants -R ../0_Reference/ref.fasta -V SL.raw_variants.vcf -selectType SNP -o SL.raw_snps.vcf &
##get raw indels
nohup $gatk -T SelectVariants -R ../0_Reference/ref.fasta -V SL.raw_variants.vcf -selectType INDEL -o SL.raw_indels.vcf &
	
##filter snp
nohup $gatk -R ../0_Reference/ref.fasta -T VariantFiltration --filterExpression 'DP < 4 || DP > 2000 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName LowQualFilter --clusterWindowSize 10 -cluster 3 --missingValuesInExpressionsShouldEvaluateAsFailing --variant ./SL.raw_snp.vcf --logging_level ERROR -o SL.snp.step1.vcf

##filter indel	
nohup $gatk -R ../0_Reference/ref.fasta -T VariantFiltration -V ./SL.raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName Filter -o SL.indel.step1.vcf

##get final VCF file
system "grep -v 'LowQual' SL.snp.step1.vcf > SL.snp.step2.vcf";
system "grep -v 'LowQual' SL.indel.step1.vcf > SL.indel.step2.vcf";
system "grep -v '^#' SL.snp.step2.vcf|awk '{print \$10}'|cut -c 1-3|sort|uniq -c > snp.stat";
system "grep -v '^#' SL.snp.step2.vcf|wc -l >> snp.stat";
SNP_count=`grep -v '^#' ${name}.snp.final.vcf|wc -l`
awk 'BEGIN{printf "%.6f\n",'${SNP_count}'/'${Genome_size}'}' >> snp.count.stat
grep -v '#' SL.snp.step2.vcf|awk '{print $10}'|awk -F ":" '{print $1}'|sort|uniq -c > homo_heter_SNP_stat

system "grep -v '^#' SL.indel.step2.vcf|awk '{print \$10}'|cut -c 1-3|sort|uniq -c>indel.stat";
system "grep -v '^#' SL.indel.step2.vcf|wc -l >> indel.stat";
echo "Insertion:" >> indel.stat
awk 'length($4)<length($5)' SL.indel.step2.vcf|grep -v "#"|wc -l >> indel.stat
echo "Deletion:" >> indel.stat
awk 'length($4)>length($5)' SL.indel.step2.vcf|grep -v "#"|wc -l >> indel.stat

InDel_count=`grep -v '^#' SL.indel.step2.vcf|wc -l`
awk 'BEGIN{printf "%.6f\n",'${InDel_count}'/'${Genome_size}'}' >> indel.count.stat



/home/genome/soft/annovar/convert2annovar.pl -format vcf4 SL.snp.step2.vcf -outfile SL.snp.avinput -include -withzyg &

mkdir annoDB
cd annoDB
cp ../../0_Reference/Reference_Sequence.fa ./
cp ../../0_Reference/Reference_Sequence.gtf ./
/home/genome/LanHsiao/bin/gtfToGenePred -genePredExt Reference_Sequence.gtf Anno_refGene.txt
或者/home/genome/LanHsiao/bin/gff3ToGenePred ref.gff Anno_refGene.txt
perl /home/genome/soft/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile ref.fasta Anno_refGene.txt --out Anno_refGeneMrna.fa
cd ..
/home/genome/soft/annovar/annotate_variation.pl -geneanno -buildver Anno SL.snp.avinput ../annoDB
echo "All SNPs:" > snp.anno.stat
awk -F "\t" '{print $1}' SL.snp.avinput.variant_function|sort|uniq -c >> snp.anno.stat
awk -F "\t" '{print $1}' SL.snp.avinput.variant_function|wc -l >> snp.anno.stat
echo "exonic SNPs:" >> snp.anno.stat
awk -F "\t" '{print $2}' SL.snp.avinput.exonic_variant_function|sort|uniq -c >> snp.anno.stat
awk -F "\t" '{print $2}' SL.snp.avinput.exonic_variant_function|wc -l >> snp.anno.stat

##INDELs annotation
/home/genome/soft/annovar/convert2annovar.pl  -format vcf4 SL.indel.step2.vcf -outfile SL.indel.avinput -include -withzyg 
/home/genome/soft/annovar/annotate_variation_bacteria.pl -geneanno -buildver Anno SL.indel.avinput annoDB/

echo "All InDels:" > indel.anno.stat
awk -F "\t" '{print $1}' SL.indel.avinput.variant_function|sort|uniq -c >> indel.anno.stat
awk -F "\t" '{print $1}' SL.indel.avinput.variant_function|wc -l >> indel.anno.stat
echo "exonic InDels:" >> indel.anno.stat
awk -F "\t" '{print $2}' SL.indel.avinput.exonic_variant_function|sort|uniq -c >> indel.anno.stat
awk -F "\t" '{print $2}' SL.indel.avinput.exonic_variant_function|wc -l >> indel.anno.stat

##Indel Length distribution
mkdir stat_result
awk -F "\t" '$6=="-"' SL.indel.avinput.variant_function|awk -F "\t" '{print $7}' > ./stat_result/InDels_length_stat_tmp1
awk '{print $1"\t"length($1)}' ./stat_result/InDels_length_stat_tmp1 > ./stat_result/InDels_length_stat_1

awk -F "\t" '$7=="-"' SL.indel.avinput.variant_function|awk -F "\t" '{print $6}' > ./stat_result/InDels_length_stat_tmp2
awk '{print $1"\t"(-1)*length($1)}' ./stat_result/InDels_length_stat_tmp2 > ./stat_result/InDels_length_stat_2

cat InDels_length_stat_1 InDels_length_stat_2> InDels_length_stat_for_plot
Rscript /home/genome/Shuyun/bin/InDel_Length.R InDels_length_stat_for_plot

mkdir 6_SV
perl /storage/RawData/Jiama/blii/bin/breakdancer-1.1.2/perl/bam2cfg.pl -q 20 -c 4 -g -h $bam >$sample.bamcfg
/storage/RawData/Jiama/blii/bin/break/breakdancer-1.1_2011_02_21/cpp/breakdancer_max -q 20 -d ./ ${name}.bamcfg >${name}.ctx
sed '/^#/d' SL.ctx |cut -f 7|sort|uniq -c > SL.stat

##一起call snp
nohup /opt/bin/jdk1.7.0_80/bin/java -jar /opt/Software/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../AKM/0_reference/ref.fa -I ../AKM/3_mapping/AKM.fina.realigned.bam -I ../INH/3_mapping/INH.realigned.bam -I ../LOF/3_mapping/LOF.fina.realigned.bam -I ../CPM/3_mapping/CPM.fina.realigned.bam -I ../SM/3_mapping/SM.fina.realigned.bam -I ../EMB/3_mapping/EMB.fina.realigned.bam -stand_call_conf 50 -stand_emit_conf 30 -glm SNP --downsample_to_coverage 5000 -o all.SNP.vcf -nt 36 &

##SNP before_after position preference
#temp file
grep -v "#" SL.snp.step2.vcf|awk -F "\t" 'length($5)==1' > ./stat_result/SNP_positon_prefer_tmp
#mutation position in chromosome
awk '{print $1"\t"$2}' ./stat_result/SNP_positon_prefer_tmp > ./stat_result/SNP_ID
#get sequence of 5 bases upstream and downstream in ref
python /home/genome/Shuyun/bin/Sequence_Choose.py ../../0_Reference/ref.fasta ./stat_result/SNP_ID ./stat_result/SNP.fasta
#remove sequence id
grep -v ">" ./stat_result/SNP.fasta > ./stat_result/SNP_Sequence
#make columns
python /home/genome/Shuyun/bin/Sequence_col.py ./stat_result/SNP_Sequence  ./stat_result/SNP_Sequence_Table
#mutation base
awk '{print $5}' ./stat_result/SNP_positon_prefer_tmp > ./stat_result/SNP_Alt
#final columns
paste ./stat_result/SNP_Sequence_Table ./stat_result/SNP_Alt  --delimiters=""|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > SNP_Sequence_final_for_plot
#unify format
sed -i 's/a/A/g' SNP_Sequence_final_for_plot
sed -i 's/g/G/g' SNP_Sequence_final_for_plot
sed -i 's/c/C/g' SNP_Sequence_final_for_plot
sed -i 's/t/T/g' SNP_Sequence_final_for_plot
Rscript /home/genome/Shuyun/bin/SNP_position_pref.R SNP_Sequence_final_for_plot
	
	
##Circle
awk -F "\t" '{print $2}' ${name}.snp.avinput.exonic_variant_function|sort|uniq -c > snp.exon.anno.stat
grep -vP "\tsynonymous" ${name}.snp.avinput.exonic_variant_function|awk -F "\t" '{print $5"\t"$6}' > SNP_ID_for_tree
awk -F "\t" '{print $5"\t"$6}' ${name}.indel.avinput.exonic_variant_function > InDel_ID_for_tree
perl  /home/genome/Shuyun/bin/SNP_anno.pl SNP_ID_for_tree SNP_tree green 
perl  /home/genome/Shuyun/bin/SNP_anno.pl InDel_ID_for_tree InDel_tree red
cat SNP_tree InDel_tree > ${name}_SNP_InDel_tree
rm SNP_ID_for_tree InDel_ID_for_tree SNP_tree InDel_tree

	cd ../5_CNV
	


# mkdir chr
# cp ../0_Reference/Reference_Sequence.fa ./chr
# cd chr
# python /storage/RawData/Jiama/blii/bin/seq_split.py Reference_Sequence.fa
# rm Reference_Sequence.fa
# cd ..
# ~/tool/CNVnator/CNVnator_v0.3/src/cnvnator -root out_${name}.root -tree ../3_Mapping/${name}.realigned.bam
# ~/tool/CNVnator/CNVnator_v0.3/src/cnvnator -root out_${name}.root -his 100 -d chr/
# ~/tool/CNVnator/CNVnator_v0.3/src/cnvnator -root out_${name}.root -stat 100
# ~/tool/CNVnator/CNVnator_v0.3/src/cnvnator -root out_${name}.root -partition 100
# ~/tool/CNVnator/CNVnator_v0.3/src/cnvnator -root out_${name}.root -call 100 > ${name}.CNV
# awk 'BEGIN{FS="\t"}{if($4<=0.05||$4>=1.8){print}}' *.CNV |wc -l >CNV.stat
# awk 'BEGIN{FS="\t"}{if($4<=0.05||$4>=1.8){print}}' *.CNV>final.CNV


# cd ../6_SV
# perl /home/genome/soft/breakdancer-1.1.2/perl/bam2cfg.pl -q 20 -c 4 -g -h  ../3_Mapping/${name}.realigned.bam > ${name}.bamcfg
# /home/genome/soft/breakdancer-1.1.2/cpp/breakdancer_max -q 20 -d ./ ${name}.bamcfg >${name}.ctx
# sed '/^#/d' ${name}.ctx |cut -f 7|sort|uniq -c > ${name}.stat


# mkdir ../7_Tree
# cd ../7_Tree

# perl /home/yewx/tool/cgview/cgview_xml_builder/cgview_xml_builder.pl -sequence Reference_Sequence.gbk -output cog.xml
# cp ../4_SNP_InDel

# java -jar /home/yewx/tool/cgview/cgview.jar -i cog.xml -o SNP_InDel.png
	
awk -F ";" '{print $1"\t"$5"\t"$6}' final|cut -f 1,4,5,9-11|less -S 
awk -F ";" '{print $1"\t"$5"\t"$6}' final|less -S 
awk -F ";" '{print $1"\t"$5"\t"$6}' final|cut -f 1,5,6|less -S
less select_gene.pl 
less -S final 
 1127  awk -F ";" '{print $1"\t"$3"}' final|less -S


