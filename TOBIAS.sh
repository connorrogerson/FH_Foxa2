#!/bin/bash
#!
#! Example SLURM job script for Peta4-Skylake (Skylake CPUs, OPA)
#! Last updated: Mon 13 Nov 12:25:17 GMT 2017
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J TOBIAS
#! Which project should be charged:
#SBATCH -A FREZZA-SL3-CPU 
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=6
#! How much wallclock time will be required?
#SBATCH --time=6:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
##SBATCH --array=12

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by core number*walltime.
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*32 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.
#! Each task is allocated 1 core by default, and each core is allocated 5980MB (skylake)
#! and 12030MB (skylake-himem). If this is insufficient, also specify
#! --cpus-per-task and/or --mem (the latter specifies MB per node).

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:

#! Full path to application executable:
application=""

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
### ATAC-seq peaks

#cd ATACorrect

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_1_2_3.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --cores $SLURM_NTASKS 

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_4_5_6.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_7_8_9.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_10_11_12.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --cores $SLURM_NTASKS

#cd ../FootprintScores

#TOBIAS FootprintScores --signal ../ATACorrect/merged_1_2_3_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --output merged_1_2_3_corrected.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_4_5_6_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --output merged_4_5_6_corrected.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_7_8_9_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --output merged_7_8_9_corrected.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_10_11_12_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --output merged_10_11_12_corrected.bw --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak --outdir BINDetect/BINDetect_Fh1_fl_cl1_+pFh_timecourse --time-series --cores $SLURM_NTASKS

### H3K27ac clusters nfr


#cd ATACorrect

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_1_2_3.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --cores $SLURM_NTASKS 

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_4_5_6.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_7_8_9.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_10_11_12.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --cores $SLURM_NTASKS

#cd ../FootprintScores

#TOBIAS FootprintScores --signal ../ATACorrect/merged_1_2_3_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --output merged_1_2_3_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_4_5_6_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --output merged_4_5_6_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_7_8_9_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --output merged_7_8_9_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_10_11_12_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --output merged_10_11_12_corrected_footprints.bw --cores $SLURM_NTASKS


#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/merged_clusters_nfr.bed --outdir BINDetect/BINDetect_Fh1_fl_cl1_+pFh_timecourse_clusters_nfr --time-series --cores $SLURM_NTASKS

### all nfrs

#cd ATACorrect

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_1_2_3.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --cores $SLURM_NTASKS 

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_4_5_6.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_7_8_9.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_10_11_12.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --cores $SLURM_NTASKS

#cd ../FootprintScores

#TOBIAS FootprintScores --signal ../ATACorrect/merged_1_2_3_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --output merged_1_2_3_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_4_5_6_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --output merged_4_5_6_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_7_8_9_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --output merged_7_8_9_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_10_11_12_corrected.bw --regions /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --output merged_10_11_12_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/all_peaks.nfr.bed --outdir BINDetect/BINDetect_Fh1_fl_cl1_+pFh_timecourse_all_nfr --time-series --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals /rds/user/cjr78/hpc-work/ATAC/tobias/H3K27ac_nfr_all/FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/Cluster_1.nfr.bed --outdir BINDetect_cluster_1_output --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals /rds/user/cjr78/hpc-work/ATAC/tobias/H3K27ac_nfr_all/FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/Cluster_2.nfr.bed --outdir BINDetect_cluster_2_output --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals /rds/user/cjr78/hpc-work/ATAC/tobias/H3K27ac_nfr_all/FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/Cluster_3.nfr.bed --outdir BINDetect_cluster_3_output --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals /rds/user/cjr78/hpc-work/ATAC/tobias/H3K27ac_nfr_all/FootprintScores/merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ChIP/H3K27ac/macs2/broad_peak/Mfuzz/Cluster_4.nfr.bed --outdir BINDetect_cluster_4_output --cores $SLURM_NTASKS

## Differential ATAC-seq regions

#cd ATACorrect

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_1_2_3.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --cores $SLURM_NTASKS 

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_4_5_6.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_7_8_9.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --cores $SLURM_NTASKS

#TOBIAS ATACorrect --bam /rds/user/cjr78/hpc-work/ATAC/bowtie/merged/merged_10_11_12.bam --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --cores $SLURM_NTASKS

#cd ../FootprintScores

#TOBIAS FootprintScores --signal ../ATACorrect/merged_1_2_3_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --output merged_1_2_3_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_4_5_6_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --output merged_4_5_6_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_7_8_9_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --output merged_7_8_9_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS FootprintScores --signal ../ATACorrect/merged_10_11_12_corrected.bw --regions /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --output merged_10_11_12_corrected_footprints.bw --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/Jaspar_non-redundant_mm.txt --signals merged_{1_2_3,7_8_9,4_5_6}_corrected_footprints.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_updown.txt --outdir BINDetect/BINDetect_Fh1_fl_cl1_+pFh_timecourse_all_nfr --time-series --cores $SLURM_NTASKS

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/a_Fhfl_cl1_2x_up_q0.05_+diffH3K27ac.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/a_Fhfl_cl1_2x_up_q0.05_+diffH3K27ac_gimme/ --cores $SLURM_NTASKS 

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/b_Fhfl_cl1_2x_up_q0.05_-diffH3K27ac.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/b_Fhfl_cl1_2x_up_q0.05_-diffH3K27ac_gimme/ --cores $SLURM_NTASKS

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/c_Fhfl_cl1_2x_down_q0.05_+diffH3K27ac.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/c_Fhfl_cl1_2x_down_q0.05_+diffH3K27ac_gimme/ --cores $SLURM_NTASKS

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/d_Fhfl_cl1_2x_down_q0.05_-diffH3K27ac.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/d_Fhfl_cl1_2x_down_q0.05_-diffH3K27ac_gimme/ --cores $SLURM_NTASKS 

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/e_Fhcl_cl1_nonsig_+diffH3K27ac_up.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/e_Fhcl_cl1_nonsig_+diffH3K27ac_up_gimme/ --cores $SLURM_NTASKS 

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/f_Fhcl_cl1_nonsig_+diffH3K27ac_down.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/f_Fhcl_cl1_nonsig_+diffH3K27ac_down_gimme/ --cores $SLURM_NTASKS 

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/g_Fhcl_cl1_nonsig_-diffH3K27ac.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/g_Fhcl_cl1_nonsig_-diffH3K27ac_gimme/ --cores $SLURM_NTASKS 

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Open_activated_alluvial_CL1/ --cores $SLURM_NTASKS

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Open_alluvial_CL1/ --cores $SLURM_NTASKS

#TOBIAS BINDetect --skip-excel --motifs /home/cjr78/miniconda3/pkgs/gimmemotifs-0.14.4-py37h516909a_0/lib/python3.7/site-packages/data/motif_databases/gimme.vertebrate.v5.0.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Activated_alluvial_CL1/ --cores $SLURM_NTASKS

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Open_activated/ 

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Open/ 

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Activated/ 

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_closed_deactivated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Closed_deactivated/ 

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_deactivated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Deactivated/

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_closed.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Closed/

#TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant.pfm --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/Genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_stable.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial/Stable/ 

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Open_activated/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_open.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Open/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_activated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Activated/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_closed_deactivated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Closed_deactivated/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_deactivated.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Deactivated/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_closed.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Closed/ --cores $SLURM_NTASKS

TOBIAS BINDetect --motifs /rds/user/cjr78/hpc-work/Motifs/gimme.vertebrate.v5.0.tobias.meme --signals /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_7_8_9_corrected.bw /rds/user/cjr78/hpc-work/ATAC/tobias/ATAC_peaks/FootprintScores/merged_1_2_3_corrected.bw  --genome /rds/user/cjr78/hpc-work/References/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --peaks /rds/user/cjr78/hpc-work/ATAC/DEseq2/Alluvial_stable.bed --outdir /rds/user/cjr78/hpc-work/ATAC/tobias/Alluvial_gimmemotifs/Stable/ --cores $SLURM_NTASKS  




