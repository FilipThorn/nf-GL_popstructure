# nf-GL_popstructure
Nextflow pipeline that calculates genotype likelihoods in angsd from a list of bamfiles and plots admixture trhough NGSadmix and PCAs through PCAngsd.

## Quick start
1) Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04)
2) Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10) 
3) Download git clone of this repository:
   ```bash
   git clone https://github.com/FilipThorn/nf-GL_popstructure
   ```
4) Download and install [`PCAngsd`](https://github.com/Rosemeis/pcangsd)
5) Run nextflowpipeline:
   ```bash
   nextflow run GL_popstr.nf --bams /PATH/TO/BAMFILELIST/'*.list' --outdir /PATH/TO/RESULTS/ --chr_ref /PATH/TO/CHRSOMELIST
   ```
## Input files

1)  bam file list example: 
    
    /results/IndvXXXX/IndvXXXX_sorted.bam<br> 
    /results/IndvXXXX/IndvXXXX_sorted.bam<br> 
    /results/IndvXXXXIndvXXXX_sorted.bam<br> 
    /results/IndvXXXX/IndvXXXX_sorted.bam<br> 

    **Lables in plots are depend on the subdirectory name**
    
    /results/**Indv0001**/Indv0001_sorted.bam<br> 
    
    if you have a different file structure you can run the pipeline with the flag **--skip_plots true** and create your plots on your own
    
2) chrosome reference file exampel:
  
   chr1<br> 
   chr2<br> 
   chr3<br> 
   chr4<br> 
   chr5<br> 
 
   **Subset of scaffolds present in your bamfiles** 
 
 
 ## HPC enviroment
 Use of a HPC is recomended. create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
 
 
##Output examples 
[Alt text](/example_plots/AsteAmay_NGSadmix_k2.png?raw=true)

  
