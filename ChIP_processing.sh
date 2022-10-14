cd /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastqc/

fastqc /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R1_001.fastq.gz
fastqc /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R2_001.fastq.gz

cd /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/ 

trimmomatic PE -threads $SLURM_NTASKS -phred33 /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R1_001.fastq.gz /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R2_001.fastq.gz ${INFILE}_R1_trimmed.fastq.gz ${INFILE}_R1_failed.fastq.gz ${INFILE}_R2_trimmed.fastq.gz ${INFILE}_R2_failed.fastq.gz ILLUMINACLIP:/rds/user/cjr78/hpc-work/scripts/TruSeq3-SE.fa:2:30:10:1:true LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:30

cd /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastqc/
fastqc /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R1_trimmed.fastq.gz
fastqc /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/merged/${INFILE}_R2_trimmed.fastq.gz

#cd /rds/user/cjr78/hpc-work/ChIP/Foxa2
#mkdir bowtie
#cd bowtie
#mkdir ${INFILE}
#cd ${INFILE}
 
bowtie2 -p $SLURM_NTASKS --dovetail -X 2000 -t -x /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/${INFILE}_R1_trimmed.fastq.gz -2 /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/${INFILE}_R2_trimmed.fastq.gz | samtools view -b -F 4 -f 2 -q30 - > ${INFILE}.bam
#convert to bam
#-b ensures a bam output
#-F 4 ensures only mapped reads are carried over (-f 2 in other instances indicates to keep correctly paired reads)
#q30 only take good quality reads
samtools view -b -F 4 -f 2 -q30 ${INFILE}.sam > ${INFILE}.bam
python /rds/user/cjr78/hpc-work/scripts/fix_isize.py ${INFILE}.bam ${INFILE}fixed.bam

#remove blacklist
intersectBed -v -abam ${INFILE}fixed.bam -b /rds/user/cjr78/hpc-work/scripts/mm10.blacklist.bed > ${INFILE}_blacklist.bam

#sort bam
samtools sort ${INFILE}_blacklist.bam > ${INFILE}_blacklist_sorted.bam

#remove duplicates

picard MarkDuplicates INPUT=${INFILE}_blacklist_sorted.bam OUTPUT=${INFILE}_blacklist_sorted_markeddup.bam METRICS_FILE=${INFILE}_blacklist_sorted_dup_stats.txt 

#sort
samtools sort ${INFILE}_blacklist_sorted_markeddup.bam > ${INFILE}_final.bam

samtools index ${INFILE}_final.bam

rm *.sam
rm *_blacklist*.bam
rm *_blacklist*.bam.bai

#Obtain reads from Drosophila genome
bowtie2 -p $SLURM_NTASKS -t -x /rds/user/cjr78/hpc-work/References/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome -1 /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/${INFILE}_R1_trimmed.fastq.gz -2 /rds/user/cjr78/hpc-work/ChIP/Foxa2/fastq/${INFILE}_R2_trimmed.fastq.gz | samtools view -bF 4 -f 2 -q30 - > ${INFILE}_dm3.bam


##scale BAM files using spike-in ratios.
cd /rds/user/cjr78/hpc-work/ChIP/Foxa2/bowtie/*_S${SLURM_ARRAY_TASK_ID}

INFILE=`awk "NR==${SLURM_ARRAY_TASK_ID}" /rds/user/cjr78/hpc-work/ChIP/Foxa2/bowtie/scale_factors.txt`

samtools view -b -s $INFILE *_final.bam > ${SLURM_ARRAY_TASK_ID}_final_scaled.bam

#call peaks

macs2 callpeak -t FLFL_Foxa2_1_S1_final_scaled.bam FLFL_Foxa2_2_S2_final_scaled.bam --keep-dup all --outdir FLFL_merged --SPMR -g mm -B -n FLFL_merged 
macs2 callpeak -t CL1_Foxa2_1_S5_final_scaled.bam CL1_Foxa2_2_S6_final_scaled.bam --keep-dup all --outdir CL1_merged --SPMR -g mm -B -n CL1_merged 


