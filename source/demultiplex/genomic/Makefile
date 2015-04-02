BIN_PATH = ../../scripts/bin/

clean:
	rm script.out
	rm errorlog.csv

run_pre_align:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_prealignment; \

run_clean_pre_align:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_cleanprealignment; \	

run_align_fwd:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_alignment_fwd; \

run_align_rev:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_alignment_rev; \	

run_post_align:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_postalignment; \	

run_precompile_align:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_precompilealignment; \

run_compile_align:
	csh $(BIN_PATH)analysis_align_ops_v200.csh run_samples_compilealignment; \		

################################################################################
################################################################################
run_sam_to_bam:
	csh ../../scripts/bin/analysis_align_ops_v200.csh run_samples_sam_2_bam; \

run_sort_bam:
	csh ../../scripts/bin/analysis_align_ops_v200.csh run_samples_sort_bam; \

run_index_bam:
	csh ../../scripts/bin/analysis_align_ops_v200.csh run_samples_index_sortedbam; \

run_bam_to_tdf:
	csh ../../scripts/bin/analysis_align_ops_v200.csh run_samples_sortedbam_2_tdf; \
