
params {
	core=4
	origout=analysis_long.Very_sensitive
	WD="."
	DV="homo_sapiens"
	DB="$HOME/data/db/$DV"
	out="$WD/$origout/$DV"
	samtools="$HOME/conda/bin/samtools "
	shortstack="$HOME/bin/ShortStack "
	rRNAdb="$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-id95.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-db:$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s-id98.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s:$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-dna-db:$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-dna-db"
}
	
process perlTask {
	output:
	stdout randNums
	shell:
	'''
#!/usr/bin/env perl
	'''
}

process pyTask {
	echo true
	input:
	stdin randNums
	
	'''
#!/usr/bin/env python
	import sys
	// for line in sys.stdin:
	'''
}

/*
	* Defines some parameters in order to specify the refence genomes
* and read pairs by using the command line options
*/
	params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
	params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
	params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
	params.outdir = 'results'
	
	log.info """\
         R N A T O Y   P I P E L I N E    
         =============================
         genome: ${params.genome}
         annot : ${params.annot}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
	.stripIndent()
	
	/*
		* the reference genome file
	*/
		genome_file = file(params.genome)
	annotation_file = file(params.annot)
	
	/*
		* Create the `read_pairs` channel that emits tuples containing three elements:
		* the pair ID, the first read-pair file and the second read-pair file 
	*/
		Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 
	
	/*
		* Step 1. Builds the genome index required by the mapping process
	*/
		process buildIndex {
			tag "$genome_file.baseName"
			
			input:
				file genome from genome_file
			
			output:
				file 'genome.index*' into genome_index
			
			"""
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
		}
	
	/*
		* Step 2. Maps each read-pair by using Tophat2 mapper tool
	*/
		process mapping {
			tag "$pair_id"
			
			input:
				file genome from genome_file 
			file annot from annotation_file
			file index from genome_index
			set pair_id, file(reads) from read_pairs
			
			output:
				set pair_id, "accepted_hits.bam" into bam
			
			"""
    tophat2 -p ${task.cpus} --GTF $annot genome.index $reads
    mv tophat_out/accepted_hits.bam .
    """
		}
	
	/*
		* Step 3. Assembles the transcript by using the "cufflinks" tool
	*/
		process makeTranscript {
			tag "$pair_id"
			publishDir params.outdir, mode: 'copy'  
			
			input:
				file annot from annotation_file
			set pair_id, file(bam_file) from bam
			
			output:
				set pair_id, file('transcript_*.gtf') into transcripts
			
			"""
    cufflinks --no-update-check -q -p $task.cpus -G $annot $bam_file
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
		}
	
	workflow.onComplete { 
		println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
	}
	