#gimme maelstrom /rds/user/cjr78/hpc-work/ChIP/H3K27ac/gimmeMotifs/mfuzz_clusters/mfuzz_clusters.txt mm10 /rds/user/cjr78/hpc-work/ChIP/H3K27ac/gimmeMotifs/mfuzz_clusters/mfuzz_clusters_gimmemotifs.out

#gimme motifs -g mm10 -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_up_q0.05.txt gimmemotifs_ATAC_up2x_wg_bg 

#gimme motifs -g mm10 -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/DEseq2/Fhfl_cl1_2x_down_q0.05.txt gimmemotifs_ATAC_down2x_wg_bg 

#gimme motifs -g mm10 -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak ../DEseq2/Fhfl_cl1_2x_up_q0.05.txt gimmemotifs_ATAC_up2x_atac_bg 

#gimme motifs -g mm10 -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak ../DEseq2/Fhfl_cl1_2x_down_q0.05.txt gimmemotifs_ATAC_down2x_atac_bg 

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/DEseq2/mfuzz_clustering/list_of_peaks_in_cluster_all.txt mm10 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_atac_fuzz_clusters

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_alluvial/Alluvial_maelstrom_input.txt mm10 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maesltrom_alluvial/

#homerTools extract /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.narrowPeak /home/cjr78/.conda/envs/seq_analysis/share/homer-4.10-0/.//data/genomes/mm10/ -fa > /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.fa

#gimme motifs -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.fa -g mm10 Alluvial_open_activated_GM.5.0.Forkhead.0008.bed Alluvial_open_activated_GM.5.0.Forkhead.0008_gimmemotifs

#gimme motifs -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.fa -g mm10 Alluvial_open_GM.5.0.Forkhead.0008.bed Alluvial_open_GM.5.0.Forkhead.0008_gimmemotifs

#gimme motifs -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.fa --tools Homer -g mm10 Alluvial_open_activated_GM.5.0.Forkhead.0008.bed Alluvial_open_activated_GM.5.0.Forkhead.0008_gimmemotifs_homer

#gimme motifs -N $SLURM_NTASKS -b /rds/user/cjr78/hpc-work/ATAC/macs2/merged_all/merged_final_fixed_peaks.fa --tools Homer -g mm10 Alluvial_open_GM.5.0.Forkhead.0008.bed Alluvial_open_GM.5.0.Forkhead.0008_gimmemotifs_homer

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_forkhead/maelstrom_forkhead_input.txt mm10 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_test/ 

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_alluvial/input.table.txt mm10 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maesltrom_alluvial_jaspar2020_14.4/ -p /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant_validated.pfm

#gimme maelstrom -N $SLURM_NTASKS --filter_cutoff 0.6 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maelstrom_alluvial/input.table.txt mm10 /rds/user/cjr78/hpc-work/ATAC/gimmemotifs/maesltrom_alluvial_updated

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ChIP/H3.3/gimmemotifs/H3.3_maelstrom_input_vst.txt mm10 /rds/user/cjr78/hpc-work/ChIP/H3.3/gimmemotifs/maesltrom_H3.3_VST/ 

#gimme maelstrom -N $SLURM_NTASKS /rds/user/cjr78/hpc-work/ChIP/H3.3/gimmemotifs/H3.3_maelstrom_input_vst.txt mm10 /rds/user/cjr78/hpc-work/ChIP/H3.3/gimmemotifs/maesltrom_H3.3_VST_jaspar/ -p /rds/user/cjr78/hpc-work/Motifs/JASPAR2020_vertebrates_non-redundant_validated.pfm

#gimme maelstrom gimmemotifs_input.txt mm10 Foxa2_maelstrom

gimme motifs /rds/user/cjr78/hpc-work/ChIP/Foxa2/macs2/merged/CL1_Foxa2_merged_SE/CL1_Foxa2_SE_summits.bed -g mm10 CL1_Foxa2_summits -N $SLURM_NTASKS
