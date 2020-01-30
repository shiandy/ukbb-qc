import os
import sys
import csv
import argparse
import numpy as np
import pandas as pd

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
        level = logging.INFO)

class set_node:
    ''' Object holding a node in the disjoint set data structure. '''
    def __init__(self, rank, dat):
        '''
        rank: Used for union by rank
        parent: The parent of this node. In the beginning, every node is
        its own parent.
        dat: Data for the node. Can be any python immutable.
        children: The children of this node, or nodes whose parent is
        this node. Also includes the current node. This is a python set.
        '''
        self.rank = rank
        self.parent = self
        self.data = dat
        self.children = set()

class DisjointSets:
    '''
    Disjoint set data structure implementing union find.
    Used to find clusters of related individuals from the pairwise
    relatedness information.
    '''
    def __init__(self):
        '''
        members is indexed by the set_node's data field, and maps
        a data point to a node. Assumes there are no overlaps in data
        points.
        '''
        self.members = dict()

    def make_set(self, dat):
        ''' Make a set from a data point. '''
        cur_node = set_node(0, dat)
        self.members[dat] = cur_node
        cur_node.children.add(cur_node)
        return cur_node

    def get(self, dat):
        ''' Return the node corresponding to a data value. '''
        if dat in self.members:
            return self.members[dat]
        else:
            return None

    def find(self, node):
        '''
        Given a node, find its root node.
        A root node is a node whose parent is itself.
        '''
        if node.parent != node:
            node.parent = self.find(node.parent)
        return node.parent

    def link(self, node1, node2):
        '''
        Changes parent pointer of the node with smaller rank to point to
        node with bigger rank.
        node1, node2: two set_node objects.
        '''
        if node1.rank == node2.rank:
            node2.rank += 1
            node1.parent = node2
            node2.children.update(node1.children)
        elif node1.rank < node2.rank:
            node1.parent = node2
            node2.children.update(node1.children)
        else:
            assert(node2.rank < node1.rank)
            node2.parent = node1
            node1.children.update(node2.children)

    def union(self, node1, node2):
        '''
        Creates a union between two nodes, so they exist under one root
        node.
        node1, node2: two set_node objects.
        '''
        self.link(self.find(node1), self.find(node2))

def load_data(pheno_file, sample_file):
    ''' Load data from phenotype and sample file, and merge them.'''
    # mapping for columns to use
    pheno_cols_touse = {
            "eid": "eid",
            "31-0.0": "phenotyped_sex",
            "22001-0.0": "genetic_sex",
            "22021-0.0": "genetic_kinship",
            "22006-0.0": "genetic_ethnicity",
            "22019-0.0": "sex_chr_aneuploidy",
            "22027-0.0": "het_missing_outliers",
            "22005-0.0": "missingness"
    }
    logging.info("Reading phenotype file...")
    ukb_phenotypes = pd.read_csv(pheno_file,
            usecols = pheno_cols_touse.keys(), engine = "c")
    ukb_phenotypes.rename(columns = pheno_cols_touse, inplace = True)

    logging.info("Reading sample file...")

    # need to skip the 2nd row of the file because it contains
    # information about the data types of the columns, which will mess
    # up inferring the data types.
    ukb_sample = pd.read_csv(sample_file, sep = " ", header = 0,
            skiprows = [1])
    assert(all(ukb_sample['ID_1'] == ukb_sample['ID_2']))

    # merge the phenotype information and sample information
    ukb_sqc = pd.merge(ukb_phenotypes, ukb_sample, left_on = "eid",
            right_on = "ID_1")
    logging.info("%d samples remain after merge" % ukb_sqc.shape[0])
    return ukb_sqc

def basic_filter(dat, miss_thresh):
    ''' Filter based on sex check, sex chromosome aneuploidy,
    heterozygosity & missingness outliers, and missingness.'''

    assert(miss_thresh >= 0 and miss_thresh <= 1)

    logging.info("%d samples had phenotyped sex != genetic sex." %
            np.sum(dat['phenotyped_sex'] != dat['genetic_sex']))
    logging.info("%d samples had >=1 relative." %
            np.sum(dat['genetic_kinship'] == 0))
    logging.info("%d samples are of white British ancestry." %
            np.nansum(dat['genetic_ethnicity'] == 1))
    logging.info("%d samples have sex chromosome aneuploidy." %
            np.nansum(dat['sex_chr_aneuploidy'] == 1))
    logging.info("%d samples are heterozygosity/missingness outliers." %
            np.nansum(dat['het_missing_outliers'] == 1))
    logging.info("%d samples have missingness >= %f" %
            (np.nansum(dat['missingness'] >= miss_thresh), miss_thresh))

    keep_arr = np.array((dat['phenotyped_sex'] == dat['genetic_sex'],
        pd.isna(dat['sex_chr_aneuploidy']),
        pd.isna(dat['het_missing_outliers']),
        dat['missingness'] < miss_thresh))

    keep_inds = np.all(keep_arr, axis = 0)
    samp_touse = dat.loc[keep_inds, :]

    logging.info("%d samples remain after filtering by sex check, "
    "sex chromosome aneuploidy, het.missing.outliers, and missingness" %
            samp_touse.shape[0])

    return(samp_touse)

def related_clusters(rel_file, sep = " ", col1 = "ID1", col2 = "ID2"):
    '''
    Find related clusters using the union-find algorithm.
    rel_file: A file with the columns col1 and col2. This file contains
    at least two columns, and each row of this file is a related pair.
    sep: Delimiter between columns in rel_file.
    col1: ID column of one of the related pair. Assumed to be an
    integer.
    col2: ID column of the other of the related pair. Assumed to be an
    integer.
    '''

    disjoint_set = DisjointSets()
    with open(rel_file, "r") as csvfile:
        csv_read = csv.DictReader(csvfile, delimiter = sep)
        for row in csv_read:
            # assume IDs are integers
            id1 = int(row[col1])
            id2 = int(row[col2])
            node1 = disjoint_set.get(id1)
            node2 = disjoint_set.get(id2)

            # Make new sets if nodes do not exist. This prevents us from
            # adding duplicates.
            if node1 is None:
                node1 = disjoint_set.make_set(id1)
            if node2 is None:
                node2 = disjoint_set.make_set(id2)
            # Since each row describes a related pair, union them.
            disjoint_set.union(node1, node2)

    # list of related clusters. Each element is a set of the node's
    # data.
    clusters = list()
    for (id_key, node) in disjoint_set.members.items():
        if node.parent == node:
            child_data = map(lambda child_node: child_node.data,
                    node.children)
            clusters.append(set(child_data))
    return(clusters)

def relatedness_filter(dat, relatedness):
    ''' Filters out related samples.'''
    logging.info("Calculating relatedness clusters...")
    rel_clusters = related_clusters(relatedness)

    # map each cluster to a cluster id
    clust_mapping = dict()
    for (clust_id, clust) in enumerate(rel_clusters):
        for eid in clust:
            clust_mapping[eid] = clust_id

    logging.info("Filtering out related samples...")
    temp = dat[dat['genetic_kinship'] <= 0]
    temp2 = dat[dat['genetic_kinship'] >= 1]

    # group by cluster ID, and then get minimum missingness per cluster
    # source:
    # https://datascience.stackexchange.com/questions/26308/after-grouping-to-minimum-value-in-pandas-how-to-display-the-matching-row-resul
    # see answer by Bon Ryu
    groups = temp2.set_index("eid").groupby(clust_mapping)
    idx_to_keep = groups['missing'].idxmin()
    related_touse = temp2.loc[temp2['eid'].isin(idx_to_keep)]

    # take people with no relatedness and add on filtered people with
    # relatives
    samp_touse = temp.append(related_touse)
    logging.info("%d samples remain after filtering out related "
        "individuals" % samp_touse.shape[0])
    return samp_touse

def parse_args():
    ''' Parse and check command line arguments.'''

    # argument parsing
    parser = argparse.ArgumentParser(description = "Filter UK Biobank Samples")
    parser.add_argument("phenotype",
            help = "Input UK Biobank phenotype csv file")
    parser.add_argument("sample",
            help = "Input UK Biobank sample file (.sample)")
    parser.add_argument("-r", "--relatedness",
            help = "Input UK Biobank relatedness file (ukb[#]_rel_[#].dat.")
    parser.add_argument("-o", "--output", default = "filtered_samples.txt",
            help = "Output file for sample IDs to keep")
    parser.add_argument("--white-british", action = "store_true",
            default = False,
            help = "Include only individuals with white British ancestry.")
    parser.add_argument("-m", "--missingness", type = float,
            default = 0.05,
            help = "Exclude samples with missingness above this number.  Default: 0.05")

    args = parser.parse_args()

    '''
    args = parser.parse_args([
        "/n/holylfs/LABS/xlin_lab_genetics/UKB/phenotypes/ukb37696.csv",
        "/n/holylfs/LABS/xlin_lab_genetics/UKB/gwas_data/imputed/ukb52008_imp_chr22_v3_s487314.sample",
        "-r",
        "/n/holylfs/LABS/xlin_lab_genetics/UKB/gwas_data/ukb52008_rel_s488282.dat",
        "--white-british"
    ])
    '''

    # error checking
    if (args.missingness < 0 or args.missingness > 1):
        sys.exit("Missingness threshold must be between 0 and 1.")

    # check if files exist
    if not os.path.isfile(args.phenotype):
        sys.exit("Phenotype file %s is not a file or does not exist." %
                args.phenotype)
    if not os.path.isfile(args.sample):
        sys.exit("Sample file %s is not a file or does not exist." %
                args.sample)
    if (not (args.relatedness is None)) and (not os.path.isfile(args.relatedness)):
        sys.exit("Relatedness file %s is not a file or does not exist."
                % args.relatedness)

    return(args)

def main():
    args = parse_args()

    print("Sample filtering arguments used:")
    for (arg_name, arg_value) in vars(args).items():
        print("%s: %s" %(arg_name, arg_value))

    # load data
    ukb_sqc = load_data(args.phenotype, args.sample)

    # sample filtering
    ukb_samp_touse = basic_filter(ukb_sqc, args.missingness)

    if args.white_british:
        # filter out non white british
        ukb_samp_touse = ukb_samp_touse.loc[ukb_samp_touse['genetic_ethnicity'] == 1]
        logging.info("%d samples remain after filtering out "
        "individuals without white British ancestry" %
        ukb_samp_touse.shape[0])

    if not (args.relatedness is None):
        ukb_samp_touse = relatedness_filter(ukb_samp_touse,
                args.relatedness)

    out_dir = os.path.dirname(args.output)
    if not os.path.isdir(out_dir):
        logging.info("Output directory %s doesn't exist, creating it."
                %out_dir)
        os.mkdir(out_dir)

    ukb_samp_touse['eid'].to_csv(args.output, header = False,
            index = False)
    logging.info("Samples to include written to file %s" %args.output)

if __name__ == "__main__":
    main()
