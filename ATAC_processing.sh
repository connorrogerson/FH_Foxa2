## trimmed fastq files were obtained from sequencing facility (Babaraham Institute)

mkdir bowtie
cd bowtie

mkdir ${SLURM_ARRAY_TASK_ID}
cd ${SLURM_ARRAY_TASK_ID}

bowtie2 -p $SLURM_NTASKS -X 2000 --dovetail -x ~/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 ~/ATAC/fastq/${SLURM_ARRAY_TASK_ID}_r1.fq.gz -2 ~/ATAC/fastq/${SLURM_ARRAY_TASK_ID}_r2.fq.gz -S ${SLURM_ARRAY_TASK_ID}.sam

samtools view -SbhF 4 -f2 -q30 $SLURM_ARRAY_TASK_ID.sam > $SLURM_ARRAY_TASK_ID.bam

intersectBed -v -abam $SLURM_ARRAY_TASK_ID.bam -b ~/scripts/mm10.blacklist.bed > ${SLURM_ARRAY_TASK_ID}_blacklist.bam

samtools sort ${SLURM_ARRAY_TASK_ID}_blacklist.bam > ${SLURM_ARRAY_TASK_ID}_blacklist_sorted.bam

java -Xmx2g -jar ~/.conda/envs/seq_analysis/share/picard-2.21.8-0/picard.jar MarkDuplicates INPUT=${SLURM_ARRAY_TASK_ID}_blacklist_sorted.bam OUTPUT=${SLURM_ARRAY_TASK_ID}_blacklist_sorted_markeddup.bam METRICS_FILE=${SLURM_ARRAY_TASK_ID}_blacklist_sorted_dup_stats

python ~/fix_isize.py ${SLURM_ARRAY_TASK_ID}_blacklist_sorted_markeddup.bam ${SLURM_ARRAY_TASK_ID}_blacklist_sorted_markeddup_fixed.bam

samtools sort ${SLURM_ARRAY_TASK_ID}_blacklist_sorted_markeddup_fixed.bam > ${SLURM_ARRAY_TASK_ID}_final_fixed.bam

samtools index ${SLURM_ARRAY_TASK_ID}_final_fixed.bam

rm *.sam
rm *blacklist.bam
rm *sorted.bam

cd /rds/user/cjr78/hpc-work/ATAC/

mkdir macs2
cd macs2
mkdir ${SLURM_ARRAY_TASK_ID}
cd ${SLURM_ARRAY_TASK_ID}

macs2 callpeak -t ~/ATAC/bowtie/${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_ID}_final_fixed.bam -n ${SLURM_ARRAY_TASK_ID}_final_fixed -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --call-summit --keep-dup all

macs2 callpeak -t ~/ATAC/bowtie/*/*_final_fixed.bam -n merged_final_fixed -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all

macs2 callpeak -t ~/ATAC/bowtie/1/1_final_fixed.bam ~/ATAC/bowtie/2/2_final_fixed.bam ~/ATAC/bowtie/3/3_final_fixed.bam -n merged_1_2_3 -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all

macs2 callpeak -t ~/ATAC/bowtie/4/4_final_fixed.bam ~/ATAC/bowtie/5/5_final_fixed.bam ~/ATAC/bowtie/6/6_final_fixed.bam -n merged_4_5_6 -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all

macs2 callpeak -t ~/ATAC/bowtie/7/7_final_fixed.bam ~/ATAC/bowtie/8/8_final_fixed.bam ~/ATAC/bowtie/9/9_final_fixed.bam -n merged_7_8_9 -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all

macs2 callpeak -t ~/ATAC/bowtie/10/10_final_fixed.bam ~/ATAC/bowtie/11/11_final_fixed.bam ~/ATAC/bowtie/12/12_final_fixed.bam -n merged_10_11_12 -q 0.01 -g hs -f BAM --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all
