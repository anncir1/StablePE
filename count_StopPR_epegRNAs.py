#!/usr/bin/env python
''' This script counts reads (FASTQ format) relative to a known list of sequences.
'''

#__author__ = "aelin@princeton.edu"
#__version__ = "0.1"

import argparse
import csv
import pandas as pd

def read_fastq(fastq):
    with open(fastq, 'r') as fastq_handle:
        record = []
        n = 0
        for line in fastq_handle:
            n += 1
            record.append(line.rstrip())
            if n == 4:
                yield record
                n = 0
                record = []

def align_reads(r1_fastq, start_end, dict_ps_seqToID):
    # Count all reads
    total_reads = 0
    # Count reads that don't match
    unmatched = 0

    dict_matches = {}
    for read in read_fastq(r1_fastq):
        total_reads += 1
        readname = read[0][1:-17] # strip first character and right most info, down to read name
        if start_end is not None:
            readseq = read[1][start_end[0]:start_end[1]] # trims bp(s) if needed (i.e., to result in 20bp protospacer)
        else:
            readseq = read[1]
        #readqs = read[3] # if quality scores are needed?
        if readseq in dict_ps_seqToID:
            dict_matches[readname] = dict_ps_seqToID[readseq]
        else:
            unmatched += 1

    print(str(total_reads-unmatched) + " reads mapped out of " + str(total_reads))

    return dict_matches, total_reads, total_reads-unmatched

def count_pegs(r1_matches, r2_matches, dict_peg_counts):
    # Count all reads
    mapped = 0
    # Count reads that don't match
    recombined = 0

    for read in r1_matches:
        if read in r2_matches: # if the read name was mapped to R2 as well
            mapped += 1
            if r1_matches[read]+"-"+r2_matches[read] in dict_peg_counts: # if there was no recombination
                dict_peg_counts[r1_matches[read]+"-"+r2_matches[read]] += 1
            else:
                recombined += 1

    print(str(mapped) + " reads matched, " + str(recombined) + " of those recombined")

    return dict_peg_counts, mapped, recombined


def main(ps_seqs, ext_seqs, ps_ext_matches, r1_fastq, r1_start_stop, r2_fastq, r2_start_stop, csv_counts, summary):
    print("Reading in dictionaries...")
    dict_ps_seqToID = pd.read_csv(ps_seqs, sep="\t").set_index('Protospacer').to_dict()['Protospacer_ID']
    dict_ext_seqToID = pd.read_csv(ext_seqs, sep="\t").set_index('Extension').to_dict()['Extension_ID']
    dict_peg_counts = dict.fromkeys(pd.read_csv(ps_ext_matches, sep="\t")['pegRNA_ID'], 0)

    print("Aligning R1...")
    r1_matches, r1_total_reads, r1_mapped_reads = align_reads(r1_fastq, r1_start_stop, dict_ps_seqToID)
    print("Aligning R2...")
    r2_matches, r2_total_reads, r2_mapped_reads = align_reads(r2_fastq, r2_start_stop, dict_ext_seqToID)

    print("Comparing read matches...")
    dict_peg_counts, combined_mapped_reads, recombined_reads = count_pegs(r1_matches, r2_matches, dict_peg_counts)

    print("Writing output files...")
    with open(csv_counts, 'w') as f:
        for key in dict_peg_counts.keys():
            f.write("%s,%s\n"%(key,dict_peg_counts[key]))
    f.close()

    with open(summary, 'w') as f:
        f.write("R1 mapped: " + str(r1_mapped_reads) + " out of " + str(r1_total_reads) + " reads.\n")
        f.write("R2 mapped: " + str(r2_mapped_reads) + " out of " + str(r2_total_reads) + " reads.\n")
        f.write("Combined R1 and R2, mapped: " + str(combined_mapped_reads) + " reads, " + str(recombined_reads) + " of which were recombined.\n")
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tool for counting reads per expected sequence (e.g., guide spacer, DNA barcode).")

    # Filepath args
    parser.add_argument('--ps_seqs', required=True, type=str,
                        help="Path to text file containing expected protospacer read sequences.")
    parser.add_argument('--ext_seqs', required=True, type=str,
                        help="Path to text file containing expected extension read sequences.")
    parser.add_argument('--ps_ext_matches', required=True, type=str,
                        help="Path to text file containing protospacer and extension combinations.")
    parser.add_argument('--r1_fastq', required=True, type=str,
                        help="Path to read 1 (protospacer) fastq file containing observed sequences.")
    parser.add_argument('--r1_start_stop', type=list, default=None,
                        help="List of start and stop indices to trim R1 sequences (optional).")
    parser.add_argument('--r2_fastq', required=True, type=str,
                        help="Path to read 2 (extension) fastq file containing observed sequences.")
    parser.add_argument('--r2_start_stop', type=list, default=None,
                        help="List of start and stop indices to trim R2 sequences (optional).")
    parser.add_argument('--csv_counts', required=True, type=str, help="Path to output csv file.")
    parser.add_argument('--summary', required=True, type=str, help="Path to summary text file.")

    args = parser.parse_args()
    main(ps_seqs=args.ps_seqs, ext_seqs=args.ext_seqs, ps_ext_matches=args.ps_ext_matches,
         r1_fastq=args.r1_fastq, r1_start_stop=args.r1_start_stop, r2_fastq=args.r2_fastq,
         r2_start_stop=args.r2_start_stop, csv_counts=args.csv_counts, summary=args.summary)
    print("Done.")