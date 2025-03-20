#!/usr/bin/env nextflow

nextflow.enable.dsl=1

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
             Dardel version (JN)
             –––––––––––––––––––––––––––––––––––––––
             'USAGE'
             1. Prepare input data
             2. Make sure params are set to your liking in nextflow.config
             3. Run (with a nextflow version capable of running dsl1)
                NXF_VER=22.10.1 nextflow run main.nf --bamlist_tsv bams.tsv --chr chr.list --outdir outfolder

             'Mandatory arguments:'
             --bamlist_tsv FILE   Path to tsv FILE with three columns:
                                  name  n  path_to_file_with_list_of_bam_files
             --chr         FILE   Path to FILE containing a list of chromosomes to use
             --outdir      FOLDER Path to output FOLDER where results will be stored

             'OPTIONS'
             --help               Outputs this help log
             -resume              Nextflow cmd to resume modified workflow

             'HPC'
             -profile      FILE   If the intention is to run the workflow on HPC, please
                                  provide a suitable profile in the nextflow.config file
                                  and in the profile folder

             'SUPPORT'
             Email Filip.Thorn@NRM.se for questions on script
             Consult http://www.popgen.dk/software/index.php/ANGSD for ANGSD
             Consult http://www.popgen.dk/software/index.php/NgsAdmix for NGSadmix
             Consult http://www.popgen.dk/software/index.php/PCAngsd for PCangsd
             """
    .stripIndent()
    exit 1
}

log.info """\
         –––––––––––––––––––––––––––––––––––––––
         GenotypeLikelihood Population Structure
         NEXTFLOW   P I P E L I N E
         Dardel version (JN)
         –––––––––––––––––––––––––––––––––––––––
         bamlist_tsv : ${params.bamlist_tsv}
         chr         : ${params.chr}
         outdir      : ${params.outdir}
         """
         .stripIndent()

// Channel
channel.fromPath(params.bamlist_tsv)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple(row.name, file(row.subset), row.ancestral) }
    .set { subset_ch }

// Make chromosome list
chromo = file(params.chr).readLines()

// Make K interval list
interval = ['-2', '-1', '0', '1', '2']

// Run angsd
process GenerateGL {

    tag "$name"
    label 'RAM_high'
    //conda 'environment.yml'
    conda '/cfs/klemming/projects/supr/nrmdnalab_storage/src/miniforge3/envs/nf-GL_popstructure'

    publishDir "${params.outdir}/01.GL/split/$name", mode: 'copy'

    input:
    tuple val(name), file(subset), val(ancestral) from subset_ch
    each chr from chromo

    output:
    tuple val(name), file("${name}_${chr}.beagle.gz"), val(ancestral) into GL_split_ch

    script:
    """
    indLen=\$(wc -l < $subset)
    angsd \
        -nThreads ${task.cpus} \
        -out ${name}_${chr} \
        -GL 1 \
            -doGlf 2 \
            -doMajorMinor 1 \
            -skipTriallelic 1 \
        -doMaf 1 \
            -minMaf $params.minMaf \
            -SNP_pval 1e-6 \
        -doCounts 1 \
            -setMinDepthInd $params.setMinDepthInd \
            -setMinDepth $params.setMinDepth \
            -minInd \$indLen \
            -minQ $params.minQ \
        -bam $subset \
            -uniqueOnly 1 \
            -minMapQ $params.minMapQ \
            -r $chr \
            -only_proper_pairs 1 \
            -remove_bads 1
    """
}

process MergeGL {

    tag "$name"
    //conda 'environment.yml'
    conda '/cfs/klemming/projects/supr/nrmdnalab_storage/src/miniforge3/envs/nf-GL_popstructure'

    publishDir "${params.outdir}/01.GL/cat/$name", mode:'copy'

    input:
    tuple val(name), file(beagle), val(ancestral) from GL_split_ch.groupTuple(by:[0,2], size:chromo.size()).map { prefix, beagle, surfix -> tuple(prefix, beagle.sort{it.name}, surfix) }

    output:
    tuple val(name), file("${name}_all.beagle.gz"), val(ancestral) into GL_merge_ch
    tuple val("${name}_prune"), file("${name}_prune.beagle.gz"), val(ancestral) into GL_prune_ch

    script:
    """
    zcat $beagle | head -n 1 > header.txt
    zcat $beagle | grep -v marker > ${name}_noheader.beagle
    cat header.txt ${name}_noheader.beagle | gzip -c > ${name}_all.beagle.gz
    awk '(NR%${params.pruneDist}==1)' ${name}_noheader.beagle > ${name}_noheader_prune.beagle
    cat header.txt ${name}_noheader_prune.beagle | gzip -c > ${name}_prune.beagle.gz
    """
}

// Split channel
if (params.prune == false) {
    GL_merge_ch.view().into { GL_pca_ch; GL_admix_ch }
}

if (params.prune == true ) {
    GL_prune_ch.view().into { GL_pca_ch; GL_admix_ch }
}

process NGSadmix {

    tag "$name"
    label 'EXTRA'
    //conda 'environment.yml'
    conda '/cfs/klemming/projects/supr/nrmdnalab_storage/src/miniforge3/envs/nf-GL_popstructure'

    publishDir "${params.outdir}/02.NGSadmix/$name", mode:'copy'

    List permutate = 1..10

    input:
    tuple val(name), file(GL), val(ancestral) from GL_admix_ch
    each number from interval
    each boot from permutate

    output:
    file("${name}_k*.qopt")
    file("*.log")

    script:
    """
    echo $number
    k=\$(($number + $ancestral))
    NGSadmix -likes $GL -K \$k -P ${task.cpus} -seed \$RANDOM -o ${name}_k\${k}_permutate${boot}
    """
}

process PCANGSD {

    tag "$name"
    label 'EXTRA'
    //conda 'environment.yml'
    conda '/cfs/klemming/projects/supr/nrmdnalab_storage/src/miniforge3/envs/nf-GL_popstructure'

    publishDir "${params.outdir}/03.PCAngsd/$name", mode:'copy'

    input:
    tuple val(name), file(GL), val(ancestral) from GL_pca_ch

    output:
    file("${name}.cov")

    script:
    """
    pcangsd -b $GL -o ${name} -t ${task.cpus}
    """
}

