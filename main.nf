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
            nextflow run main.nf --bamlist_tsv --outdir /PATH/TO/RESULTS/ --chr /PATH/TO/chr.list
         
            'Mandatory arguments:'
            --bams         FILE      Path to file containing a list of bam paths. Extention .list  
            --outdir       PATH      Path to output directory where results will be should be stored 
            --chr          FILE      Path to file containing a subset of chromosomse present in bam files 
            
            'OPTIONS'
            --help                   Outputs this help log      
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
         outdir       	: ${params.outdir}
         chr      	: ${params.chr}
         """
         .stripIndent()


// Channel
	channel.fromPath(params.bamlist_tsv)
		.splitCsv(header:true, sep:'\t')
		.view()
		.map { row -> tuple(row.name, file(row.subset), row.ancestral ) }
		.set { subset_ch }

// make chromosome list
chromo = file(params.chr).readLines()

//make K interval list
interval = ['-1', '0', '1']

process GenerateGL {

   label 'HIGH_RAM'
   
   tag "$name"   
    
   publishDir "${params.outdir}/01.GL/split/$name", mode:'copy'

   input:
   tuple val(name), file(subset), val(ancestral) from subset_ch
   each chr from chromo

   output:
   tuple val(name), file("${name}_${chr}.beagle.gz"), val(ancestral)  into GL_split_ch

   script:
   """
   indLen=\$( wc -l $subset | awk '{print \$1}')	
    
   
    angsd \
	-nThreads ${task.cpus} \
	-bam $subset \
	-out ${name}_${chr} \
	-r $chr \
	-uniqueOnly 1 \
	-minMapQ 20 \
	-minQ 20 \
	-GL 2 \
	-doGlf 2 \
	-doMajorMinor 1 \
	-skipTriallelic 1 \
	-doMaf 1 \
	-minMaf 0.05 \
	-SNP_pval 1e-6 \
	-doCounts 1 \
	-minInd \$indLen \
	-setMinDepthInd 2 \
	-setMinDepth 20
   """
}


process MergeGL {

   tag "$name"
   
   publishDir "${params.outdir}/01.GL/cat/$name", mode:'copy'

   input:
   tuple val(name), file(beagle), val(ancestral) from GL_split_ch.groupTuple(by:[0,2], size:chromo.size()).map { prefix, beagle, surfix -> tuple( prefix, beagle.sort{it.name}, surfix) }

   output:
   tuple val(name), file("${name}_all.beagle.gz"), val(ancestral) into GL_merge_ch
   tuple val("${name}_prune"), file("${name}_prune.beagle.gz"), val(ancestral) into GL_prune_ch

   script:
   """
   zcat $beagle | head -n 1 > header.txt

   zcat $beagle | grep -v marker > ${name}_noheader.beagle 
   cat header.txt ${name}_noheader.beagle | gzip -c > ${name}_all.beagle.gz
   
   #Psuedo Pruning
   awk '(NR%${params.pruneDist}==1)' ${name}_noheader.beagle > ${name}_noheader_prune.beagle
   cat header.txt ${name}_noheader_prune.beagle | gzip -c > ${name}_prune.beagle.gz
   """
}

// split channel

GL_merge_ch.mix(GL_prune_ch).view().into { GL_pca_ch; GL_admix_ch }


process NGSadmix {
   
   tag "$name"

   label 'EXTRA'
   
   publishDir "${params.outdir}/02.NGSadmix/$name", mode:'copy'
   
   input:
   tuple val(name), file(GL), val(ancestral) from GL_admix_ch
   each number from interval

   output:
   file("${name}_k*.qopt")

   script:
   """
   echo $number
   
   k=\$(($number + $ancestral))

   NGSadmix -likes $GL -K \$k -P ${task.cpus} -o ${name}_k\${k}
   """
} 

process PCANGSD {

   tag "$name"
   
   label 'EXTRA'
   
   publishDir "${params.outdir}/03.PCAngsd/$name", mode:'copy'
   
   input:
   tuple val(name), file(GL), val(ancestral) from GL_pca_ch

   output:
   file("${name}.cov") 
   
   script:
   """
   pcangsd.py -beagle $GL -o ${name} -threads ${task.cpus}
   """
}


