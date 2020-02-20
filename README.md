# UK Biobank Imputed Genotype QC

This is the UK Biobank QC procedure for the imputed genotypes.

## QC Procedure

This is the QC procedure as of 2020-01-31.

### Sample QC

1. Start: 487314 samples
2. 367 removed for failing sex check (genotyped sex != phenotyped sex)
3. 651 removed for sex chromosome aneuploidy
4. After retaining samples with only white british ancestry, 408186
   samples remain.
5. To filter out related samples, we used the UK Biobank provided
   relatedness file to identify clusters of relatives. For each cluster,
   we selected 1 individual who had the lowest missingness rate. After
   this step, 332430 samples remained.

### Variant QC

We filter based on MAF, HWE p-value, and missingness rate for each SNP.
For all the SNPs in the imputed dataset, we re-compute the MAF, HWE
p-value, and missingness rate on the sample of unrelated white british
ancestry individuals. We then filter out SNPs with:

1. missingness > 5%
2. INFO < 0.5
3. MAF < 1e-4
4. HWE p-value < 1e-10

The number of variants for each chromosome is reported in the
`variant_counts.txt` file.


## Usage

The QC is implemented using
[Snakemake](https://snakemake.readthedocs.io/en/stable/). Snakemake is a
workflow management system that helps create reproducible and scalable
data analyses. Snakemake can automatically determine what steps to run.
It won't repeat steps that are already run, which saves computation
time, and also keeps track of what needs to be run, saving programmer
time.

### Installation

This QC pipeline depends on:

1. Python 3
2. Snakemake
3. Numpy and Pandas
4. R
5. gds2bgen R package
6. Plink
7. Plink2

### Installation of Python Components

On the FAS cluster, you need to first load the Anaconda module with:

    module load Anaconda3/5.1.0

On your local computer, make sure
[Miniconda3](https://docs.conda.io/en/latest/miniconda.html) is
installed.

You should install the Python components into an isolated software
environment using Anaconda:

    conda create -c bioconda -c conda-forge -n ukbb-qc snakemake-minimal numpy pandas

Then, this conda environment needs to be activated:

    conda activate ukbb-qc
    snakemake --help

### Installation of gds2bgen

Install [gds2bgen](https://github.com/zhengxwen/gds2bgen) by following
the instructions here: https://github.com/zhengxwen/gds2bgen

On the FAS Cluster, you will need to load R

    module load R/3.6.1-fasrc02

### Plink and Plink2

We also need Plink and Plink2. Plink is used to calculate statistics on
the genotype calls, while Plink2 is needed to calculate statistics on
the imputed data and to do the bulk of the filtering. Plink is available
on the FAS Cluster using

    module load plink/1.90-fasrc01

For Plink2, you can download the executable from the [Plink2
website](https://www.cog-genomics.org/plink/2.0/). On the FAS Cluster, I
used the Linux AVX2 Intel version.

    $ plink2 --version
    PLINK v2.00a2LM AVX2 Intel (28 Nov 2019)

### Configuration

The commands that get run for the QC live in `Snakefile`. This is a
special file, similar to a Makefile, that gets read by Snakemake and
automatically executed. The output directory for the QC-ed files are in
this directory and needs to be tweaked before you can run the QC.
There are also some QC parameters and input directories in this file
that you might want to change.

#### Specifying the output

The `out_root` variable under the `CONFIGURATION` section controls the
output of the QC-ed data. Please change this to the directory where you
want to put the data.

#### Specifying the input

There are 6 paths that you may need to change. These are under the
section `Paths to the data` in `Snakefile`.

1. `ukbb_root`: Root directory to UK Biobank data. The phenotype file
    should live in the phenotypes subdirectory.
2. `gwas_root`: Root directory for UK Biobank GWAS data. The sample
   relatedness file should live here.
3. `imputed_root`: The imputed genotypes should live here.
4. `calls_root`: The genotype chip data should live here.
5. `calls_base`: Basename for the chip data file (i.e. without the .bed,
   .bim, .fam extension). We merged all the chip data for each chromosme
   together into one file, and this is the basename for the merged file.

TODO: This could be cleaned up a bit.

#### QC Configuration

The QC parameters are located near the top of the `Snakefile`, under the
section `CONFIGURATION`. Here, you can control the different QC
thresholds. Underneath that, there is a section where you can control the
input file folder.

#### Cluster Configuration

Snakemake uses the cluster to run some long commands and to run many
commands at once. The cluster we use is the FAS Cluster, which uses the
Slurm scheduler. The cluster configuration is found in `cluster.json`.
It is a JSON file. You should change the partition to your own lab's
partition instead of the Lin Lab partitions. For more information on
cluster configurations, see the [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#snakefiles-cluster-configuration)
or this [tutorial](https://hpc-carpentry.github.io/hpc-python/17-cluster/).


### Running the QC

First, you need to make sure that Anaconda, plink 1.90, and R are
activated on the FAS cluster. Make sure that you can call the plink2
executable. Then, activate the `conda` environment.

The QC process uses python scripts to generate lists of variants and
samples to keep / exclude. Then, it uses plink2 to perform calculations
and filtering. First, the UK Biobank bgen files are converted to the
plink2 pgen format. Then, the MAF and HWE p-value are calculated. Next,
we use plink2 to output a finalized bgen file that contains the QC-ed
variants and samples. Finally, this bgen file is converted to a gds
file.

The master file used by Snakemake to execute the commands is the
`Snakefile`.

To get a list of commands to run, use

    snakemake -pn

This does a dry run and prints out all of the commands that need to be
run. Snakemake will automatically look for a file called `Snakefile` to
use. If you renamed your `Snakefile` to something else, you can specify
the file to use using the `-s` flag, e.g.

    snakemake -pn -s my_snakefile

This QC was run on the FAS cluster. We have "local" commands which can
be run on an interactive node, and then we have "cluster" commands which
get sent off to the Slurm scheduler by Snakemake to be run.

First, run the local commands (on an interactive node) with

    snakemake localall

or

    snakemake localall -s my_snakefile

To run the "cluster" commands, use

    nohup snakemake -j $N --cluster-config cluster.json \
        --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -n {cluster.n} -N {cluster.N} -p {cluster.partition}" 2>&1 &

Change `$N` to the number of jobs to submit simultaneously. The jobs are
parallelized by chromosome, so 24 is a good choice. Snakemake will
automatically submit the jobs to the cluster. The cluster configuration
is found in `cluster.json`. See the Configuration section below for
configuration options.

#### Nohup

It will take some time to go through all the cluster commands. Above, we
ran Snakemake with nohup, so you can safely log out of the cluster and
your commands will keep running. The output of Snakemake is stored in
the file `nohup.out` by default. Each command that Snakemake runs also
writes its own log file, which are stored in the `logs` directory. You
can monitor those logs to get a sense of the progress. Additionally,
cluster jobs will generate a `slurm-#######.out` file, where `######` is
the job ID. This file can also be checked to monitor progress.

#### Screen / Tmux

Alternatively, instead of using nohup, you can use screen or tmux to
keep your snakemake command going. Once you detach your screen or tmux
session, you can log out of the cluster. Then, you can log back in and
attach your session. Take note of which login node you used. If you just
ssh to login.rc.fas.harvard.edu, you will be automatically redirected to
one of the multiple login nodes. Instead, use

    uname -a

to find out which login node you're currently on before you log out.
When you log back in, ssh specifically to that login node, e.g.

    ssh username@boslogin01.rc.fas.harvard.edu
