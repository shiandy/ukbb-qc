# =====================================================================
# UK Biobank QC Pipeline
# =====================================================================
#
#
import os

# =====================================================================
# Constants
# =====================================================================
# list of chromosomes
CHROM = list(range(1, 23)) + ["X", "XY"]

# =====================================================================
# CONFIGURATION
# =====================================================================
# directory for output files
# if you want to use $SCRATCH/directory, replace $SCRATCH with
# os.environ['SCRATCH']
#out_root = os.path.join(os.environ['SCRATCH'], "xlin", "ashi",
#    "uk-biobank-qc")
# PLEASE EDIT THIS!
out_root = "/path/to/your/output/directory"

# sample missingness: filter out samples who are missing >= this
# proportion of variants. Must be a number between 0 and 1.
samp_missing = 0.05

# variant missingness: filter out variants that are missing above this
# threshold in samples. Must be a number between 0 and 1.
var_missing = 0.05

# variant INFO filter: filter out variants with INFO score less than
# this cutoff.
var_info = 0.5

# variant HWE filter: filter out variants with HWE p-values less than
# this cutoff.
var_hwe = 1e-10

# variant MAF filter: filter out variants with MAF less than this
# threshold.
var_maf = 1e-4

# For MAF filtering: set this to True to use the MAF computed by plink
# on the unrelated individuals of white British ancestry. Otherwise,
# if this is False, uses the MAF computed in the MFI files.
use_plink_maf = True

# Set to True if you wish to only include individuals with white British
# ancestry in your final QC-ed dataset.
use_white_british = True

# Set to True if you wish to only include unrelated individuals in your
# final QC-ed dataset.
use_unrelated = True

# =====================================================================
# Paths to the data
# =====================================================================
ukbb_root = "/n/holylfs/LABS/xlin_lab_genetics/UKB/"
gwas_root = os.path.join(ukbb_root, "gwas_data")
imputed_root = os.path.join(gwas_root, "imputed")
calls_root = os.path.join(gwas_root, "calls")
calls_base = os.path.join(calls_root, "ukb_cal_allchr_v2")

sample_file = os.path.join(imputed_root, "ukb52008_imp_chr22_v3_s487314.sample")

localrules: filter_samples, imputed_stats, sample_symlinks,
    filter_variants, snp_qc_filter, all, localall

rule all:
    input:
        expand(os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.gds"), chrom = CHROM)
# run snakemake localall to run all the local (non-cluster) jobs.
rule localall:
    input:
        incl_samples = os.path.join(out_root, "incl_samples.txt"),
        incl_samples_hwe = os.path.join(out_root, "incl_samples_hwe.txt"),
        autosomes_samples = expand(os.path.join(imputed_root,"ukb52008_imp_chr{chrom}_v3.sample"), chrom = range(1, 23)),
        x_samples = os.path.join(imputed_root, "ukb52008_imp_chrX_v3.sample"),
        xy_samples = os.path.join(imputed_root, "ukb52008_imp_chrXY_v3.sample"),
        excl_lists = expand(os.path.join(out_root, "excl_imp_chr{chrom}.txt"),
            chrom = CHROM),
        snp_qc = os.path.join(out_root, "excl_snp_qc.txt")


# =====================================================================
# Sample filtering
# =====================================================================

# Use these samples for your final QC-ed dataset.

# configure the relatedness flag
def make_rel_flag(wildcards, input):
    if use_unrelated:
        rel_flag = "--relatedness " + input.relatedness
    else:
        rel_flag = ""
    return(rel_flag)

rule filter_samples:
    input:
        phenotype = os.path.join(ukbb_root, "phenotypes", "ukb37696.csv"),
        sample = sample_file,
        relatedness = os.path.join(gwas_root, "ukb52008_rel_s488282.dat")
    output: os.path.join(out_root, "incl_samples.txt")
    params:
        samp_missing = samp_missing,
        white_british = "--white-british" if use_white_british else "",
        relatedness = make_rel_flag
    log:
        "logs/filter_samples/filter_samples.log"
    shell:
        '''
        python filter_samples.py {input.phenotype} {input.sample} \
            {params.relatedness} \
            --missingness {params.samp_missing} \
            {params.white_british} \
            --output {output} > {log} 2>&1
        '''

# Use these samples for plink calculations (HWE, MAF, missingness). By
# default, we get rid of samples with the wrong sex, with excessive
# missingness, and with sex chromosome aneuploidy. We also only retain
# unrelated white British individuals.
rule filter_samples_hwe:
    input:
        phenotype = os.path.join(ukbb_root, "phenotypes", "ukb37696.csv"),
        sample = sample_file,
        relatedness = os.path.join(gwas_root, "ukb52008_rel_s488282.dat")
    output: os.path.join(out_root, "incl_samples_hwe.txt")
    params:
        samp_missing = samp_missing
    log:
        "logs/filter_samples/filter_samples_hwe.log"
    shell:
        '''
        python filter_samples.py {input.phenotype} {input.sample} \
            --relatedness {input.relatedness} \
            --missingness {params.samp_missing} \
            --white-british \
            --output {output} > {log} 2>&1
        '''

# =====================================================================
# Variant filtering
# =====================================================================
# run statistics on genotype calls
rule calls_stats:
    input:
        incl_samps = os.path.join(out_root, "incl_samples_hwe.txt"),
        bed = calls_base + ".bed",
        bim = calls_base + ".bim",
        fam = calls_base + ".fam"
    output:
        hwe = calls_base + "_stats.hwe",
        maf = calls_base + "_stats.frq",
        imiss = calls_base + "_stats.imiss",
        lmiss = calls_base + "_stats.lmiss"
    params:
        input_prefix = calls_base,
        output_prefix = calls_base + "_stats",
        memory_mb = 8000
    threads: 4
    shell:
        '''
        plink -bfile {params.input_prefix} \
            --keep-fam {input.incl_samps} \
            --hardy --freq --missing \
            --memory {params.memory_mb} \
            --threads {threads} \
            --out {params.output_prefix}
        '''

# run all the HWE, missingness, and frequency calculations for the
# imputation data
rule imputed_stats:
    input:
        hwe = expand(os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.hardy"), chrom = CHROM),
        freq = expand(os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.afreq"), chrom = CHROM),
        smiss = expand(os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.smiss"), chrom = CHROM),
        vmiss = expand(os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.vmiss"), chrom = CHROM)

# symlink sample files. Makes it easier to run batch jobs on all of
# them. Caution: autosomes and x/xy sample files have different names.
rule sample_symlinks:
    input:
        sample_auto = os.path.join(imputed_root, "ukb52008_imp_chr22_v3_s487314.sample"),
        sample_x = os.path.join(imputed_root, "ukb52008_imp_chrX_v3_s486663.sample"),
        sample_xy = os.path.join(imputed_root, "ukb52008_imp_chrXY_v3_s486349.sample")
    output:
        autosomes = expand(os.path.join(imputed_root,"ukb52008_imp_chr{chrom}_v3.sample"),
        chrom = range(1, 23)),
        x = os.path.join(imputed_root, "ukb52008_imp_chrX_v3.sample"),
        xy = os.path.join(imputed_root, "ukb52008_imp_chrXY_v3.sample")
    run:
        for dest in output.autosomes:
            os.symlink(input.sample_auto, dest)
        os.symlink(input.sample_x, output.x)
        os.symlink(input.sample_xy, output.xy)

# convert bgen files to pgen files for plink
rule convert_pgen:
    input:
        bgen = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3.bgen"),
        sample = os.path.join(imputed_root, "ukb52008_imp_chr{chrom}_v3.sample")
    output:
        # output gets sent to a file on the scratch disk by default.
        plink_pgen = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pgen"),
        plink_pvar = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pvar.zst"),
        plink_psam = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.psam")
    threads: 4
    params:
        memory_mb = 8000,
        output_prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
    shell:
        '''
        plink2 --bgen {input.bgen} ref-first \
            --sample {input.sample} \
            --memory {memory_mb} --threads {threads} \
            --make-pgen vzs --out {params.output_prefix}
        '''

# run statistics on imputation data using plink
# need to separate out X chromosome since it generates hardy.x file.
rule imputation_chrX_stats:
    input:
        plink_pgen = os.path.join(out_root, "ukb_imp_chrX_v3.pgen"),
        plink_pvar = os.path.join(out_root, "ukb_imp_chrX_v3.pvar.zst"),
        plink_psam = os.path.join(out_root, "ukb_imp_chrX_v3.psam"),
        incl_samps = os.path.join(gwas_root, "incl_samples_hwe.txt")
    output:
        hwe = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats.hardy"),
        freq = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats.afreq"),
        smiss = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats.smiss"),
        vmiss = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats.vmiss")
    threads: 4
    params:
        input_prefix = os.path.join(out_root, "ukb_imp_chrX_v3"),
        output_prefix = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats"),
        hwe_x = os.path.join(imputed_root, "ukb_imp_chrX_v3_stats.hardy.x"),
        memory_mb = 8000
    shell:
        '''
        plink2 --pfile {params.input_prefix} vzs \
            --memory {params.memory_mb} --threads {threads} \
            --keep-fam {input.incl_samps} \
            --hardy --freq --missing \
            --out {params.output_prefix}
        ln -s {params.hwe_x} {output.hwe}
        '''

rule imputation_chr_stats:
    input:
        plink_pgen = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pgen"),
        plink_pvar = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pvar.zst"),
        plink_psam = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.psam"),
        incl_samps = os.path.join(gwas_root, "incl_samples_hwe.txt")
    output:
        hwe = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.hardy"),
        freq = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.afreq"),
        smiss = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.smiss"),
        vmiss = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.vmiss")
    threads: 4
    params:
        input_prefix = lambda wildcards, input: os.path.splitext(input[0])[0],
        output_prefix = lambda wildcards, output: os.path.splitext(output[0])[0],
        memory_mb = 8000
    shell:
        '''
        plink2 --pfile {params.input_prefix} vzs \
            --memory {params.memory_mb} --threads {threads} \
            --keep-fam {input.incl_samps} \
            --hardy --freq --missing \
            --out {params.output_prefix}
        '''

# generate list of variants to exclude
# do this using statistics from the imputed data

# configure the freq flag
def make_freq_flag(wildcards, input):
    if use_plink_maf:
        freq_flag = "--freq " + input.freq
    else:
        freq_flag = ""
    return(freq_flag)

rule filter_variants:
    input:
        mfi = os.path.join(imputed_root, "ukb_mfi_chr{chrom}_v3.txt"),
        hwe = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.hardy"),
        freq = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.afreq"),
        vmiss = os.path.join(imputed_root, "ukb_imp_chr{chrom}_v3_stats.vmiss")
    output: os.path.join(out_root, "excl_imp_chr{chrom}.txt")
    params:
        info = var_info,
        hwe = var_hwe,
        maf = var_maf,
        missing = var_missing,
        freq = make_freq_flag
    log:
        "logs/filter_variants/filter_variants_chr{chrom}.log"
    shell:
        '''
        python filter_variants.py {input.mfi} {input.hwe} {input.vmiss} \
            {params.freq} \
            --output {output} \
            --info {params.info} \
            --hwe {params.hwe} \
            --maf {params.maf} \
            --missingness {params.missing} > {log}
        '''

# and also with the genotyped SNPs QC filter
rule snp_qc_filter:
    input: os.path.join(gwas_root, "ukb_snp_qc.txt")
    output: os.path.join(out_root, "excl_snp_qc.txt")
    log:
        "logs/filter_variants/snp_qc_filter.log"
    shell:
        '''
        python filter_snp_qc.py {input} --output {output} > {log}
        '''

# =====================================================================
# Final filtering
# =====================================================================

rule final_filter:
    input:
        plink_pgen = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pgen"),
        plink_pvar = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.pvar.zst"),
        plink_psam = os.path.join(out_root, "ukb_imp_chr{chrom}_v3.psam"),
        incl_samples = os.path.join(out_root, "incl_samples.txt"),
        excl_variants = os.path.join(out_root, "excl_imp_chr{chrom}.txt"),
        excl_snpqc = os.path.join(out_root, "excl_snp_qc.txt")
    output:
        bgen = os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.bgen"),
        sample = os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.sample")
    log:
        "logs/final_filter/final_filter_chr{chrom}.log"
    threads: 4
    params:
        # get plink input / output prefix
        input_prefix = lambda wildcards, input: os.path.splitext(input[0])[0],
        output_prefix = lambda wildcards, output: os.path.splitext(output[0])[0],
        memory_mb = 8000
    shell:
        '''
        plink2 --pfile {params.input_prefix} vzs \
            --memory {params.memory_mb} --threads {threads} \
            --keep-fam {input.incl_samples} \
            --exclude {input.excl_variants} {input.excl_snpqc} \
            --export bgen-1.3 \
            --out {params.output_prefix}
        '''

# =====================================================================
# Convert bgen files to GDS
# =====================================================================

rule bgen_to_gds:
    input:
        bgen = os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.bgen"),
        sample = os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.sample")
    output: os.path.join(out_root, "ukb_imp_chr{chrom}_v3_qc.gds")
    log:
        "logs/bgen_to_gds/bgen_to_gds_chr{chrom}.log"
    threads: 16
    shell:
        '''
        Rscript --vanilla bgen2gds.R --bgen {input.bgen} \
            --output {output} \
            --threads {threads} > {log} 2>&1
        '''
