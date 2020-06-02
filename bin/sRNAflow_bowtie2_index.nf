#!/usr/bin/env NXF_HOME=.nextflow/home NXF_WORK=.nextflow/work nextflow -log data/results/Cristina3/nextflow.log run sRNAflow.nf -resume --name Cristina3 --input data/test/*.fa -c bin/sRNAflow.config -with-dag data/results/Cristina3/flowchart.html -with-timeline data/results/Cristina3/timeline.html -with-report data/results/Cristina3/report.html

workflow.runName = "${params.name}"

genomes  = Channel.from( "homo_sapiens","mus_musculus" )
datasets = Channel.fromPath("data/test/*.fa").map { file -> tuple(file.simpleName, file) }

results_path = "${workflow.projectDir}/data/results/${params.name}"

//	output:
//	file "${genome}.fa" into genomes_fasta
// Step 4. Builds the genome
process build_bowtie2_index {
	storeDir 'db/genomes'
	input:
	file genome from genomes_fasta
    file "sRNAflow_bowtie2_index.R" from file("bin/sRNAflow_bowtie2_index.R")
	tag "$genome.baseName"
	output:
	file "${genome}.*.bt2*" into genome_index
	"""
#!/usr/bin/env Rscript --vanilla ${genome} ${task.cpus}
source("sRNAflow_bowtie2_index.R")
	"""
}
