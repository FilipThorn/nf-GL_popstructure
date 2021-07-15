#!/usr/bin/env nextflow

if (params.help) {
    log.info """\
            __________________
            |                |
            | |```````| |`````
            | |____   | |
            |     |   | |
            | |````   | |
            | |ilip   | |hörn     
            –––––––––––––––––––––––––––––––––––––––
            GenotypeLikelihood Population Structure 
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            nextflow run GL_popstr.nf --bams /PATH/'*.list' --outdir /PATH/TO/RESULTS/ --chr_ref /PATH/TO/CHRSOMELIST 
         
            'Mandatory arguments:'
            --bams         FILE      Path to file containing a list of bam paths. Extention .list  
            --outdir       PATH      Path to output directory where results will be should be stored 
            --chr_ref      FILE      Path to file containing a subset of chromosomse present in bam files 
            
            'OPTIONS'
            --help                   Outputs this help log      
            --k            INTEGER   Number of Ancestral populations to test. Defaults to k=2   
            --skip_plots   BOOLEAN   Run pipeline without plots (true/false)       
            -resume                  Nextflow cmd to resume modified workflow
            
            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 

            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            Consult http://www.popgen.dk/software/index.php/ANGSD for ANGSD
            Consult http://www.popgen.dk/software/index.php/NgsAdmix for NGSadmix
            Consult http://www.popgen.dk/software/index.php/PCAngsd for PCangsd
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         GenotypeLikelihood Population Structure 
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         bams         : ${params.bams}
         outdir       : ${params.outdir}
         k            : ${params.k}
         chr_ref      : ${params.chr_ref}
         """
         .stripIndent()


// Channel
   Channel.fromPath( params.bams)
         .ifEmpty { error "Cannot find any path matching: ${params.bams}" }
         .map { it -> [it.name -  ~/\.list/, it] }
         .set { input }
        
         input.into{ bams_ch; bams_list1_ch; bams_list2_ch }

process Genotypelikelihoods {

   label 'RAM_high'
    
   publishDir "${params.outdir}/01.GL/$subset", mode:'copy'

   input:
   tuple val(subset), file(bam) from bams_ch

   output:
   tuple val(subset), file("${subset}.beagle.gz") into GL_ch

   script:
   """
   angsd -nThreads ${task.cpus} -bam $bam -rf ${params.chr_ref} \
   -uniqueOnly 1 -minMapQ 20 -minQ 20 -GL 2 -doGlf 2 \
   -doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -minMaf 0.05 \
   -SNP_pval 1e-6 -doCounts 1 -setMinDepth 20 -out $subset
   """
}

// Split channel into multiple
   GL_ch.into { GL_ch1; GL_ch2 }    

process NGSadmix {
   
    label 'RAM_high'
   
   publishDir "${params.outdir}/02.NGSadmix/$subset", mode:'copy'

   List k_list = 2..params.k
   
   input:
   tuple val(subset), file(GL) from GL_ch1
   each anc from k_list

   output:
   tuple val(subset), file("${subset}_k${anc}.qopt") into admixture_ch

   script:
   """
   NGSadmix -likes $GL -K $anc -P ${task.cpus} -o ${subset}_k${anc}
   """
} 
 
process PCAngsd {

  label 'RAM_high'

  publishDir "${params.outdir}/03.PCAngsd/$subset", mode:'copy'

  input:
  set val(subset), file(GL) from GL_ch2

  output:
  set val(subset), file("${subset}.cov") into covariance_ch

  script:
  """
  ${params.PCAngsd} -beagle $GL -o ${subset} -threads ${task.cpus}
  """
}

if (!params.skip_plots){
  
    admixture_ch.combine( bams_list2_ch, by: 0 )
    .set { admix_comb }

    process NGSadmix_plot {

    label 'FAST'
    
    publishDir "${params.outdir}/04.plots/$subset", mode:'copy'
    
    List k_list = 2..params.k

    input:
    each anc from k_list
    tuple val(subset), file("${subset}_k${anc}.qopt"), file(name) from admix_comb

    output:
    file("*.pdf")

    script:
    """
    #!/usr/bin/env Rscript 

    library("ggplot2") 

    bam_list<-read.table("$name", header = FALSE)

    pop<-data.frame(indiv=character(0))
      for( i in 1:length(bam_list[["V1"]])){
        line<-strsplit(bam_list[i,], split = '/')
        name<-rev(unlist(line))[2]
        pop[i,]<-name
    } 

    colpanel <- c("antiquewhite3", "azure3", "cadetblue", "chartreuse3", "cornflowerblue", "darkgoldenrod3", "darkolivegreen3", "darkorchid3", "deeppink1", "deepskyeblue3")

    admix <- t(as.matrix(read.table("${subset}_k${anc}.qopt")))
    K <- nrow(admix)
    ord<-order(pop[,1])

    plot_nam<-paste0("${subset}_NGSadmix_k",K ,".pdf")
    pdf(plot_nam)
      bar<-barplot(admix[,ord],col=colpanel[1:K],space=0,border=NA,ylab="Admixture proportion")
      axis(1, labels = pop[ord,1], at = bar,     las = 2, cex.axis = 0.6)

    dev.off()

    """
  }

covariance_ch.combine( bams_list1_ch, by: 0 )
    .set { covariance_comb }   
 
  process PCAngsd_plot {
    
    label 'FAST'
    
    publishDir "${params.outdir}/04.plots/$subset", mode:'copy'

    input:
    tuple val(subset), file("${subset}.cov"), file(name) from covariance_comb

    output:
    file("${subset}_PCA.pdf")


    script:
    """
    #!/usr/bin/env Rscript  
    
    library("ggplot2")
    library("ggrepel")

    bam_list<-read.table("${name}", header = FALSE )
    
    pop<-data.frame(indiv=character(0))
    for( i in 1:length(bam_list[["V1"]])){
      line<-strsplit(bam_list[i,], split = '/')
      name<-rev(unlist(line))[2]
      pop[i,]<-name
    }
    pop<-unlist(list(pop[["indiv"]]))
    
    C <- as.matrix(read.table("${subset}.cov"))
    e <-as.data.frame(eigen(C)[["vectors"]])
    
    p = ggplot(aes_(x=e[["V1"]], y=e[["V2"]]), data=e)+geom_point()  + 
      theme_classic() +
      ggtitle("${subset}") +  xlab("PC1") + ylab("PC2") +
      geom_label_repel(aes_(x=e[["V1"]], y=e[["V2"]], label=pop),
                        point.padding = 0.1, label.size=0.1, box.padding=0.35,
                        label.padding = 0.1, show.legend = FALSE, size=4,
                        min.segment.length = 0.1, max.overlaps = 100)

    plot(p)   
   
    ggsave(
      "${subset}_PCA.pdf",
      plot = last_plot(),
      device = "pdf",
      scale=1,
      width = 15,
      height = 15,
      unit= "in")
    """
  }
  
}

