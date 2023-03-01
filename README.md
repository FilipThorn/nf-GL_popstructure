#V2.0 update to include psuedo linkage pruning and speed optimisation update

# nf-GL_popstructure
Nextflow pipeline that calculates genotype likelihoods in angsd from a list of bamfiles and admixture through NGSadmix and PCAs through PCAngsd.

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

1)  sample.tsv example:

    name    subset  ancestral
    SUBSET1     /Absolute/PATH/SUBSET1/bam.list   n<br>
    SUBSET2    /Absolute/PATH/SUBSET2/bam.list  n<br>
    
    
***name = subsets name***<br>
***subset = subsets bamlist location***<br>
***ancestral = clustering to fit. n +1 -1 will be calculated***<br>
    

2)  bam.list example: 
    
    /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam<br> 
    /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam<br> 
    /Absolute/PATH/IndvXXXXIndvXXXX_sorted.bam<br> 
    /Absolute/PATH/IndvXXXX/IndvXXXX_sorted.bam<br> 
    
    
3) chrosome reference file exampel:
  
   chr1<br> 
   chr2<br> 
   chr3<br> 
   chr4<br> 
   chr5<br> 
 
   **Subset of scaffolds present in your bamfiles** 
 
 
 ## HPC enviroment
 Use of a HPC is recomended. Create a nextflow config profile that matches your cluster set-up [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles)
 

  
