#!/usr/bin/env NXF_HOME=.nextflow/home NXF_WORK=.nextflow/work nextflow -log data/results/Cristina3/nextflow.log run sRNAflow.nf -resume --name Cristina3 --input data/test/*.fa -c bin/sRNAflow.config -with-dag data/results/Cristina3/flowchart.html -with-timeline data/results/Cristina3/timeline.html -with-report data/results/Cristina3/report.html

workflow.runName = "${params.name}"

genomes  = Channel.from( "homo_sapiens","mus_musculus" )
datasets = Channel.fromPath("data/test/*.fa").map { file -> tuple(file.simpleName, file) }

results_path = "${workflow.projectDir}/data/results/${params.name}"

// Closure saveClosure = {
//     file -> relative_target_dir = file.tokenize(".")[0]
//     relative_target_path = relative_target_dir + "/" + file
// }

// fq_queue = Channel.fromPath("data/test/*.fq","data/test/*.fastq").map { file -> tuple(file.simpleName, file) }

// Step x. Builds the genome
process downloadGenome {
	storeDir 'db/genomes'
	input:
	val genome from Channel.from( "homo_sapiens","mus_musculus" )
    file "sRNAflow_downloadGenome.R" from file("bin/sRNAflow_downloadGenome.R")
	output:
	file "${genome}.fa" into genomes_fasta
	file "${genome}.gtf" into genomes_gtf
	"""
#!/usr/bin/env Rscript --vanilla 
specie <- "${genome}"
print(specie)
source("sRNAflow_downloadGenome.R")
	"""
}

// Step x. Builds the genome
process build_bowtie2_index {
	storeDir 'db/genomes'
	input:
	file genome from genomes_fasta
    file "sRNAflow_bowtie2_index.sh" from file("bin/sRNAflow_bowtie2_index.sh")
	tag "$genome.baseName"
	output:
	file "${genome}.*.bt2" into genome_index
	"""
	source ./sRNAflow_bowtie2_index.sh ${genome} ${task.cpus}
	"""
}

// // Step x. ShortStack
// process ShortStack {
// 	storeDir 'bin'
// 	input:
//     file "ShortStack.pl" from file("https://raw.githubusercontent.com/zajakin/ShortStack/master/ShortStack")
// 	output:
// 	file "ShortStack.pl" into ShortStack
// 	"""
// #!/usr/bin/perl
// #do ./ShortStack.pl --threads ${task.cpus}
// 	"""
// }

/*
process tutorial_shell {
    echo true
    // container 'rocker/r-apt:bionic'
    // containerOptions = '--volume /data/db:/db'
    input:
        file 'tutorial.sh' from file("bin/tutorial.sh")
        set datasetID, file(datasetFile) from datasets
    output:
        set datasetID, file("${datasetID}.size") into sizes_files
    publishDir "$results_path/${datasetID}"  //, mode: 'move', saveAs: saveClosure
//    script:
    shell:
    '''
    echo !{params.name}
    cat !{datasetFile} | wc -l > !{datasetID}.size
    source ./tutorial.sh
   '''
}

process tutorial_R {
    publishDir "${workflow.projectDir}/data/output/", mode: "move"
    echo true
    container "rocker/r-apt:bionic"
    containerOptions = "--volume `pwd`:/db"
    input:
        file "tutorial.R" from file("bin/tutorial.R")
        file aaa from sizes_files
    // output:
    //     file "foo"
    '''
#!/usr/bin/env Rscript --vanilla 
source("tutorial.R")
    '''
}
*/