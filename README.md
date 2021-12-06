# nf-GL_popstructure

Nextflow pipeline that calculates genotype likelihoods in angsd from a list of
bamfiles and plots admixture through NGSadmix and PCAs through PCAngsd.

**NOTE: please see the documentation for branch `jn` for further instructions**

## Quick start

1) Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04)
2) Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10) 
3) Download (git clone) this repository:
   ```bash
   git clone https://github.com/FilipThorn/nf-GL_popstructure
   ```
4) Download and install [`PCAngsd`](https://github.com/Rosemeis/pcangsd)
5) Run nextflowpipeline:
   ```bash
   nextflow run GL_popstr.nf --bams /PATH/TO/BAMFILELIST/'*.list' --outdir /PATH/TO/RESULTS/ --chr_ref /PATH/TO/CHROMOSOMELIST
   ```

## Input files

1)  bam file list example: 

        /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam
        /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam
        /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam
        /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam

**Lables in plots are based on the subdirectory name**

/results/**Indv0001**/Indv0001_sorted.bam

if you have a different file structure you can run the pipeline with the flag **--skip_plots true** and create your plots on your own

2) chrosome reference file example:

        chr1
        chr2
        chr3
        chr4
        chr5

**Subset of scaffolds present in your bamfiles** 

## HPC enviroment

Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)

