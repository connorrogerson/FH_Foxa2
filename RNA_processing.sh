## Fh1 KO RNA-seq (Sciacovelli et al)
##remove adapters, trim low quality bases

INFILE=`awk "NR==${SLURM_ARRAY_TASK_ID}" ~/RNA-seq/Sciacovelli_et_al_2017/SLX-9386/samplesheet.txt`

cutadapt -q 15 -m 25 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $INFILE.trimmed.fq.gz $INFILE.fq.gz

#INFILE=`awk "NR==${SLURM_ARRAY_TASK_ID}" ~/RNA-seq/Sciacovelli_et_al_2017/SLX-9386/mm10_bam/in_samplesheet.txt`
#OUTFILE=`awk "NR==${SLURM_ARRAY_TASK_ID}" ~/RNA-seq/Sciacovelli_et_al_2017/SLX-9386/mm10_bam/out_samplesheet.txt`

STAR --runThreadN $SLURM_NTASKS --genomeDir ~/References/Mus_musculus/UCSC/mm10/Sequence/STARIndex/ --readFilesIn $INFILE --outFileNamePrefix $OUTFILE --outFilterMultimapNmax 100 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $SLURM_NTASKS

#removes unmapped reads and keeps only high quality alignments 
#samtools view -bF 4 -q1 ${OUTFILE}Aligned.sortedByCoord.out.bam > $OUTFILE.mapped.q1.bam
samtools index $OUTFILE.mapped.q1.bam
samtools view $OUTFILE.mapped.q1.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX -o $OUTFILE.final.bam
samtools index $OUTFILE.final.bam

featureCounts -T $NSLOTS -t exon -g gene_id  -Q 30 -a ~/.local/share/genomes/mm10/mm10.annotation.gtf -o counts_mm10.txt Aligned.sortedByCoord.out.bam

##siFoxa2 RNA-seq

INFILE=`awk "NR==${SLURM_ARRAY_TASK_ID}" ~/RNA-seq/siFoxa2/samplesheet.txt`

cutadapt -q 15 -m 25 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o fastq/$INFILE.trimmed.fastq.gz fastq/$INFILE.fastq.gz

fastqc -o ~/RNA-seq/siFoxa2/fastqc/ fastq/$INFILE.trimmed.fastq.gz

mkdir STAR
mkdir STAR/${SLURM_ARRAY_TASK_ID}
cd STAR/${SLURM_ARRAY_TASK_ID}

STAR --runThreadN $SLURM_NTASKS --genomeDir ~/References/Mus_musculus/UCSC/mm10/Sequence/STARIndex/ --readFilesIn ~/RNA-seq/siFoxa2/fastq/$INFILE.trimmed.fastq.gz --outFileNamePrefix ${SLURM_ARRAY_TASK_ID}. --outFilterMultimapNmax 100 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $SLURM_NTASKS

samtools view -bF 4 -q1 ${SLURM_ARRAY_TASK_ID}.Aligned.sortedByCoord.out.bam > $INFILE.mapped.q1.bam
samtools index $INFILE.mapped.q1.bam
samtools view $INFILE.mapped.q1.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX -o $INFILE.final.bam
samtools index $INFILE.final.bam

featureCounts -t exon -g gene_id -s 2 -a ~/References/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o siFoxa2_counts_mm10.txt \
~/RNA-seq/siFoxa2/STAR/1/FLFL_siNT_1_S1_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/FLFL_siNT_2_S2_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/FLFL_siNT_3_S3_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/1/FLFL_siFoxa2_1_S4_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/FLFL_siFoxa2_2_S5_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/FLFL_siFoxa2_3_S6_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/1/CL1_siNT_1_S7_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/CL1_siNT_2_S8_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/CL1_siNT_3_S9_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/1/CL1_siFoxa2_1_S10_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/CL1_siFoxa2_2_S11_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/CL1_siFoxa2_3_S12_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/1/CL19_siNT_1_S13_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/CL19_siNT_2_S14_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/CL19_siNT_3_S15_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/1/CL19_siFoxa2_1_S16_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/2/CL19_siFoxa2_2_S17_R1_001.final.bam \
~/RNA-seq/siFoxa2/STAR/3/CL19_siFoxa2_3_S18_R1_001.final.bam \
