import os
import sys
import argparse
import numpy as np
import pandas as pd

def parse_args():
    ''' Parse and check command line arguments.'''

    # argument parsing
    parser = argparse.ArgumentParser(description = "Filter UK Biobank variants using ukb_snp_qc.txt file.")
    parser.add_argument("snp_qc_file",
            help = "Input SNP QC file, supplied by UK Biobank.")
    parser.add_argument("-o", "--output", default = "snp_qc_excl.txt",
            help = "Output file of SNPs to exclude. Default: snp_qc_excl.txt.")

    args = parser.parse_args()
    '''
    args = parser.parse_args(["/n/holylfs/LABS/xlin_lab_genetics/UKB/gwas_data/ukb_snp_qc.txt"])
    '''

    # input checking
    if not os.path.isfile(args.snp_qc_file):
        sys.exit("SNP QC file %s is not a file or does not exist."
                %args.snp_qc_file)

    return(args)

def main():
    args = parse_args()
    snp_qc_df = pd.read_csv(args.snp_qc_file, sep = " ")

    # define column names for QC columns
    batch_qc_cols = ["Batch_b" + str(batch_id).zfill(3) + "_qc" for
            batch_id in range(1, 96)]
    ukbileavax_qc_cols = ["UKBiLEVEAX_b" + str(batch_id) + "_qc" for
                batch_id in range(1, 12)]

    qc_cols = batch_qc_cols + ukbileavax_qc_cols
    snp_qc_cols = np.array(snp_qc_df[qc_cols])
    snp_idx = np.logical_not(np.all(snp_qc_cols, axis = 1))
    assert(snp_idx.shape[0] == snp_qc_df.shape[0])
    snp_excl = snp_qc_df.loc[snp_idx, :].rs_id

    num_removed = snp_excl.shape[0]
    num_orig = snp_qc_df.shape[0]
    print("Removing %d / %d (%0.1f%%) SNPs in SNP QC file for failing batch QC."
            %(num_removed, num_orig, 100 * num_removed / num_orig))
    print("Writing excluded SNPs to file %s..." %args.output)
    snp_excl.to_csv(args.output, index = False, header = False)

if __name__ == "__main__":
    main()
