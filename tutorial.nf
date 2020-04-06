vegetable_datasets = Channel
                .fromPath("data/test/*.fa")
                .map { file -> tuple(file.baseName, file) }
                
results_path = "$PWD/results"
process clustalw2_align {
    publishDir "$results_path/$datasetID"
    container 'fbartusch/clustalw2:latest'
    input:
    set datasetID, file(datasetFile) from vegetable_datasets

    output:
    set datasetID, file("${datasetID}.aln") into aligned_files

    script:
    """
    clustalw2 -INFILE=${datasetFile} -OUTFILE=${datasetID}.aln
    """
}

//  workflow.projectDir
// NXF_HOME=.nextflow/home NXF_WORK=.nextflow/work nextflow -log .nextflow/my.log run tutorial.nf -c bin/tutorial.config  -with-timeline .nextflow/timeline.html -with-report .nextflow/report.html

