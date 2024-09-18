# StablePE
This repository includes scripts used to analyze endogenous and StopPR screen sequencing data for Cirincione and Simpson et al.

For endogenous prime editing experiments (analyzing CRISPRessoBatch output), example commands used to call analyze_endogenous_editing.py are below.

For analyzing HEK3 endogenous site (related to Fig. 1 and Extended Data Fig. 1):
analyze_endogenous_editing.main(filepath_to_sample_CRISPRessoBatch_output_folder, filepath_to_desired_output_file,
                                sample_name, "", 20, 'A', 'GGCCCAGACTGAGCACGTGA', True, [9, 28], ['A', 'G'], 0)

For analyzing HSPA9 endogenous site (related to Fig. 3 and Extended Data Fig. 4):
analyze_endogenous_editing.main(filepath_to_sample_CRISPRessoBatch_output_folder, filepath_to_desired_output_file,
                                sample_name, "", 16, 'A', 'AATCCCTTCAAATTGAGCAC', False, 0, "", 0)



For StopPR prime editing dropout screen, example command used to call count_StopPR_epegRNAs.py is below.

count_StopPR_epegRNAs.main("StopPR_protospacer_sequences.txt", "StopPR_extension_sequences.txt", "StopPR_ps_ext_pairs.txt", 
                           filepath_to_read1_fastq, [1, -7], filepath_to_read3_fastq, [53, -8], 
                           filepath_to_desired_full_output_csv, filepath_to_desired_counts_output_file)
