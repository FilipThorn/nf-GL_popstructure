# nf-GL_popstructure

[Nextflow](https://www.nextflow.io/) pipeline that calculates genotype
likelihoods in [ANGSD](https://github.com/ANGSD/angsd) from a list of bam files,
admixture through [NGSadmix](https://github.com/aalbrechtsen/NGSadmix) and
PCAs through [PCAngsd](https://github.com/Rosemeis/pcangsd).

v2.0 includes pseudo-linkage pruning (by thinning) and a speed-optimization update.

(**TODO**: More description here of the actual workflow (repeated k's, repeated runs, etc).

## Quick start

(**TODO**: Update install instructions in file [INSTALL](INSTALL))

1. Install prerequisites (see file [INSTALL](INSTALL)).

2. Download or git clone this repository:

        $ git clone https://github.com/FilipThorn/nf-GL_popstructure

3. Prepare input files ([see below](#input-files))

4. Run the nextflow pipeline (**Note**: needs a nextflow version that can run
   dsl1.).

        $ nextflow run main.nf \
            --bamlist_tsv sample.tsv \
            --outdir outfolder \
            --chr_ref chr.list

## Input files

1. `sample.tsv`: Tab separated file with three columns: `name`: named set of
   bam files to analyze, `subset`: file with absolute paths to bam files in the
   named set, `ancestral`: a number `n` (integer) with starting number of
   ancestral populations (`k`). Analyses will be repeated with
   `k={n-2,n-1,n,n+1,n+2}`.

        name    subset  ancestral
        SUBSET1 /path/to/SUBSET1/bam.list n
        SUBSET2 /path/to/SUBSET2/bam.list n

2.  `bam.list`: File with absolute paths to bam files in the set(s) defined in `sample.tsv`

        /path/to/IndvXXXX/IndvXXXX_sorted.bam
        /path/to/IndvXXXX/IndvXXXX_sorted.bam
        /path/to/IndvXXXXIndvXXXX_sorted.bam
        /path/to/IndvXXXX/IndvXXXX_sorted.bam

3. `chr.list`: Subset of chromosomes or scaffolds present in your bam files. Example:

        chr1
        chr2
        chr3
        chr4
        chr5

## Output files

Example:

        output
        ├── 01.GL
        │   ├── cat
        │   │   ├── mt_k1to5
        │   │   │   ├── mt_k1to5_all.beagle.gz
        │   │   │   └── mt_k1to5_prune.beagle.gz
        │   │   └── mt_k6to10
        │   │       ├── mt_k6to10_all.beagle.gz
        │   │       └── mt_k6to10_prune.beagle.gz
        │   └── split
        │       ├── mt_k1to5
        │       │   └── mt_k1to5_OZ187420.1.beagle.gz
        │       └── mt_k6to10
        │           └── mt_k6to10_OZ187420.1.beagle.gz
        ├── 02.NGSadmix
        │   ├── mt_k1to5_prune
        │   │   ├── mt_k1to5_prune_k1_permutate1.log
        │   │   ├── mt_k1to5_prune_k1_permutate1.qopt
        │   │   ├── mt_k1to5_prune_k1_permutate2.log
        │   │   ├── mt_k1to5_prune_k1_permutate2.qopt
        ...
        │   │   ├── mt_k1to5_prune_k9_permutate10.log
        │   │   └── mt_k1to5_prune_k9_permutate10.qopt
        ...
        │   └── mt_k6to10_prune
        │       ├── mt_k6to10_prune_k10_permutate1.log
        │       ├── mt_k6to10_prune_k10_permutate1.qopt
        │       ├── mt_k6to10_prune_k10_permutate2.log
        │       ├── mt_k6to10_prune_k10_permutate2.qopt
        ...
        │       ├── mt_k6to10_prune_k9_permutate10.log
        │       └── mt_k6to10_prune_k9_permutate10.qopt
        └── 03.PCAngsd
            ├── mt_k1to5_prune
            │   └── mt_k1to5_prune.cov
            └── mt_k6to10_prune
                └── mt_k6to10_prune.cov

(**TODO**: describe the output files)

## Plotting the output

For plotting the output, we highly recommend pong
(<https://github.com/ramachandran-lab/pong>).

(**TODO**: add a worked example on using pong)

## HPC environment

Use of a HPC is highly recommended. Create a nextflow config profile that
matches your cluster set-up
[`profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
and add it in the folder `profile/` and to the
[`nextflow.config`](nextflow.config).

### Example running on [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/dardel-hpc-system)

(**TODO**: finish these instructions)

Current version (20 Mar 2025) uses hard coded paths in `main.nf` to the conda
environment on dardel.  Hence, the conda environment needs to be build before
execution, and paths adjusted.

The example below assumes mandatory arguments are provided in the
[`nextflow.config`](nextflow.config) file.

    $ screen -S glpop
    $ ml PDC/23.12 bioinfo-tools Nextflow/22.10.1
    $ export NXF_OPTS='-Xms1g -Xmx4g'
    $ export NXF_HOME=/cfs/klemming/projects/supr/nrmdnalab_storage/src/NFX_HOME
    $ nextflow run \
        -w $SNIC_TMP/nf-GL_popstructure/work \
        main.nf \
        -name GL_popstructure \
        -with-report GL_popstructure.html \
        -profile dardel \
        --project 'naiss2024-22-1518'

