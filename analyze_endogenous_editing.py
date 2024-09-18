# Author: Ann Cirincione

import pandas as pd
import argparse
import re
import numpy as np
from os.path import exists


# Returns:
# 0 = unedited, no indels
# 1 = edit, no indels
# 2 = indels
def classify_edit(align_seq, ref_seq, e_pos, e, s_present, s_pos, s, ignore_bp):
    has_edit = False
    has_indel = False

    # Has correct edit, will replace with N to signify this
    # If edit includes multiple bp substitutions:
    if type(e_pos) == list:
        has_edit = True
        for i in range(len(e_pos)):
            if align_seq[e_pos[i]] == e[i]:
                align_seq = align_seq[:e_pos[i]] + 'N' + align_seq[e_pos[i] + 1:]
            else:
                has_edit = False
    
    # If edit has only 1 bp substitution:
    else:
        if align_seq[e_pos] == e:
            has_edit = True
            align_seq = align_seq[:e_pos] + 'N' + align_seq[e_pos+1:]
    
    # If there is a SNP, will replace with N so it is not mistaken as a substitution
    if s_present:
        # If there is more than 1 SNP, go through each in the list:
        if type(s_pos) == list:
            for i in range(len(s_pos)):
                if align_seq[s_pos[i]] == s[i]:
                    align_seq = align_seq[:s_pos[i]] + 'N' + align_seq[s_pos[i] + 1:]
        # If only 1 SNP:
        else:
            if align_seq[s_pos] == s:
                align_seq = align_seq[:s_pos] + 'N' + align_seq[s_pos+1:]
    
    # If flagged to ignore a certain number of bp away from the cut/nick site (to account for high PCR error), change all bp outside that range to N
    # Note: not used in any current analyses for stable PE experiments
    if ignore_bp != 0:
        num_ignored = int((50-ignore_bp*2)/2)
        align_seq = 'N'*num_ignored + align_seq[num_ignored:len(align_seq)-num_ignored] + 'N'*num_ignored

    # Removing any indices that have at least one N to collapse algnments and determine how many have indels
    N_indices = [m.start(0) for m in re.finditer('N', align_seq)]
    align_seq_no_N = "".join([char for idx, char in enumerate(align_seq) if idx not in N_indices])
    ref_seq_no_N = "".join([char for idx, char in enumerate(ref_seq) if idx not in N_indices])

    # Assigning classification for this read
    if align_seq_no_N != ref_seq_no_N:
        has_indel = True

    if has_indel:
        classification = 2
    elif has_edit:
        classification = 1
    else:
        classification = 0

    return classification


# Returns: compiled prime editing efficiencies for full input file; calls above classify_edit function per read
def analyze_alignment_file(align_df, e_pos, e, s_present, s_pos, s, ignore_bp):
    align_df['edit_type'] = pd.Series(
        classify_edit(
            row.Aligned_Sequence, row.Reference_Sequence, e_pos, e, s_present, s_pos, s, ignore_bp)
        for row in align_df.itertuples()
    )

    edit_efficiencies = align_df.groupby('edit_type', as_index=False)['pct_Reads'].agg({'read_sum': sum})

    return edit_efficiencies


def main(input_folder, output_path, sample_name, extra_sample_id, edit_pos, edit, sgRNA_seq, snp_present, snp_pos, snp, ignore_outside_n_bp):
    if extra_sample_id != "":
        extra_sample_id = "_" + extra_sample_id

    # Read in input files (output from CRISPResso)
    rep1_alignments = pd.read_csv(input_folder + "_1" + extra_sample_id + "/Alleles_frequency_table_around_sgRNA_" + sgRNA_seq + ".txt",
                                  sep='\t').rename(columns={'#Reads': 'num_Reads', '%Reads': 'pct_Reads'})
    rep2_alignments = pd.read_csv(input_folder + "_2" + extra_sample_id + "/Alleles_frequency_table_around_sgRNA_" + sgRNA_seq + ".txt",
                                  sep='\t').rename(columns={'#Reads': 'num_Reads', '%Reads': 'pct_Reads'})

    # Analyzing prime editing efficiencies
    rep1_efficiencies = analyze_alignment_file(rep1_alignments, edit_pos, edit, snp_present, snp_pos, snp, ignore_outside_n_bp)
    rep2_efficiencies = analyze_alignment_file(rep2_alignments, edit_pos, edit, snp_present, snp_pos, snp, ignore_outside_n_bp)

    # If there are triplicate files for each sample
    if exists(input_folder + "_3" + extra_sample_id + "/Alleles_frequency_table_around_sgRNA_" + sgRNA_seq + ".txt"):
        rep3_alignments = pd.read_csv(input_folder + "_3" + extra_sample_id + "/Alleles_frequency_table_around_sgRNA_" + sgRNA_seq + ".txt",
                                      sep='\t').rename(columns={'#Reads': 'num_Reads', '%Reads': 'pct_Reads'})

        rep3_efficiencies = analyze_alignment_file(rep3_alignments, edit_pos, edit, snp_present, snp_pos, snp, ignore_outside_n_bp)

        merged_efficiencies = pd.merge(rep1_efficiencies, rep2_efficiencies,
                                       on='edit_type').merge(rep3_efficiencies, on='edit_type').rename( columns={'read_sum_x': 'rep_1_reads',
                                                                                                                 'read_sum_y': 'rep_2_reads',
                                                                                                                 'read_sum': 'rep_3_reads'})
    else:
        merged_efficiencies = pd.merge(rep1_efficiencies, rep2_efficiencies,
                                       on='edit_type').rename(columns={'read_sum_x': 'rep_1_reads',
                                                                       'read_sum_y': 'rep_2_reads'})

    merged_efficiencies['mean_pct'] = merged_efficiencies.drop('edit_type', axis=1).mean(axis=1)
    merged_efficiencies['std_pct'] = merged_efficiencies.drop(['edit_type', 'mean_pct'], axis=1).std(axis=1)
    merged_efficiencies['sample'] = sample_name + extra_sample_id

    merged_efficiencies.to_csv(output_path, mode='a+', sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze CRISPResso batch output to determine prime editing efficiency.")
    # Filepath args
    parser.add_argument("--input_folder", type=str,
                        help="Full path to input directory -- use name for one sample, replicates will be averaged (i.e. path to CRISPResso_on_R_11)")
    parser.add_argument("--output_path", type=str,
                        help="Full path to output file for editing efficiencies")
    parser.add_argument("--sample_name", type=str,
                        help="i.e. HEK3")
    parser.add_argument("--extra_sample_id", type=str, default="",
                        help="i.e. d7")
    parser.add_argument("--edit_pos", type=int,
                        help="Position of expected edit (currently only substitution edits incorporated)")
    parser.add_argument("--edit", type=str,
                        help="Expected edit (currently only substitution edits incorporated)")
    parser.add_argument("--sgRNA_seq", type=str,
                        help="Guide sequence/protospacer used in CRISPResso alignment")
    parser.add_argument("--snp_present", type=bool, default=False,
                        help="True if there is an expected SNP in the wild-type genome")
    parser.add_argument("--snp_pos", type=int, default=0,
                        help="Position of SNP if present")
    parser.add_argument("--snp", type=str, default="",
                        help="Expected SNP if present")
    parser.add_argument("--ignore_outside_n_bp", type=int, default=0,
                        help="To account for high PCR errors. Ex: inputted value of 10 -> ignore all mismatches 10bp away from the cut/nick site but converting them to N first")


    args = parser.parse_args()
    main(input_folder=args.input_folder, output_path=args.output_path, sample_name=args.sample_name, extra_sample_id=args.extra_sample_id,
         edit_pos=args.edit_pos, edit=args.edit, sgRNA_seq=args.sgRNA_seq, snp_present=args.snp_present, snp_pos=args.snp_pos, snp=args.snp,
         ignore_outside_n_bp=args.ignore_outside_n_bp)
    print("Done.")