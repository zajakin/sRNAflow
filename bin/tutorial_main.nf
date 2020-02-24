#!/usr/bin/env nextflow

params {
	core=4
	origout="analysis"
	WD="."
	DV="homo_sapiens"
	DB="$HOME/data/db/$DV"
	genomefile="$DB/$DV.fa"
	out="$WD/$origout/$DV"
	samtools="$HOME/conda/bin/samtools "
	shortstack="$HOME/bin/ShortStack "
	rRNAdb="$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-id95.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-db:$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s-id98.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s:$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-dna-db:$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-dna-db"
}
	

	
/*
 * Step x. Builds the genome
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
 * Step x. 
 */
process mycount {
	tag "$genome_file.baseName"
	input:
	file genome from genome_file
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	echo "number of mapped reads Htseq" #  >> $shdir/stat.txt
	htseq_opt="htseq-count -f bam -a 0 -s no" #--additional-attr=gene_name
	$htseq_opt -q $shdir/$f.bam $DB/$DV.gtf > $shdir/htseq-count_$f.txt  #  -s reverse
	printf '%s\t%i\n' $DV `cat $shdir/htseq-count_$f.txt | grep -v "^__" | awk '{s+=$2} END {print s}'`
	stLog="$shdir/htseq_${f}__no_feature.txt $shdir/htseq_${f}__ambiguous.txt $shdir/htseq_${f}__too_low_aQual.txt $shdir/htseq_${f}__not_aligned.txt $shdir/htseq_${f}__alignment_not_unique.txt" 
	rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
	touch $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
	type="miRBase_hairpin_mergedFeatures"
# for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lincRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
	for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","lincRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","vaultRNA_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
	do
		cat $DB/gtf_biotypes/$type.gtf >> $shdir/tmp.gtf
# (grep "^#" file.gtf; grep -v "^#" file.gtf | sort -k1,1 -k4,4n) > file_sorted.gtf
		$htseq_opt -q $shdir/$f.sam $shdir/tmp.gtf --samout=$shdir/priority.sam > $shdir/tmp.txt  #  -s reverse
		
		grep -v "XF:Z:__no_feature" $shdir/priority.sam | grep -v "XF:Z:__not_aligned" > $shdir/priority_$f.$type.sam
		grep -v "^@" $shdir/priority_$f.$type.sam | awk '{print $1}' > $shdir/tmpReads.txt
# egrep -v -f $shdir/testedID.txt $shdir/tmp.txt > $shdir/htseq-count_$f.$type.txt # | sed '$ s/\n//' | tr '\n' '|'
		grep -v -Ff $shdir/testedID.txt $shdir/tmp.txt > $shdir/htseq-count_$f.$type.txt
		grep -v -Ff $shdir/testedReads.txt $shdir/tmpReads.txt > $shdir/priority_$f.$type.txt
		echo "$type `grep __no_feature $shdir/tmp.txt  | awk '{print $2}'`" >> $shdir/htseq_${f}__no_feature.txt
		echo "$type `grep __ambiguous  $shdir/tmp.txt  | awk '{print $2}'`" >> $shdir/htseq_${f}__ambiguous.txt
		echo "$type `grep __too_low_aQual $shdir/tmp.txt  | awk '{print $2}'`" >> $shdir/htseq_${f}__too_low_aQual.txt
		echo "$type `grep __not_aligned $shdir/tmp.txt  | awk '{print $2}'`" >> $shdir/htseq_${f}__not_aligned.txt
		echo "$type `grep __alignment_not_unique $shdir/tmp.txt  | awk '{print $2}'`" >> $shdir/htseq_${f}__alignment_not_unique.txt
		printf '%s\t%i\n' "$type" `cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | awk '{s+=$2} END {print s}'`
	cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | awk '{print $1}' >> $shdir/testedID.txt
		cat $shdir/priority_$f.$type.txt >> $shdir/testedReads.txt
		cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" >> $shdir/htseq-count_$f.priority.txt
		done
	grep -v -Ff $shdir/testedID.txt $shdir/htseq-count_$f.txt > $shdir/htseq-count_$f.Other_types.txt
	cat $shdir/htseq-count_$f.Other_types.txt | awk '$2 > 0' > $shdir/htseq-count_$f.Other_types.morethen0.txt
	cat $shdir/htseq-count_$f.Other_types.morethen0.txt | awk '{print $1}' > $shdir/tmp.txt
	grep -Ff $shdir/tmp.txt $DB/$DV.gtf > $shdir/htseq-count_$f.Other_types.gtf
	rename $shdir/tmp.txt $shdir/htseq-count_$f.All.txt
	rm $shdir/tmp.gtf
# GTFmiRNA=`ls $GTFPATH*.gtf | grep "miRNA"`
	grep -vh -P "^__" $shdir/htseq-count_$f.piRNAdb_mergedFeatures.txt $shdir/htseq-count_$f.piRNAbank_mergedFeatures.txt $shdir/htseq-count_$f.piRBase_mergedFeatures.txt  > $shdir/htseq-count_$f.all_piRNA_mergedFeatures.txt
	grep -vh -P "^__" $shdir/htseq-count_$f.lincRNA_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_hc_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_mergedFeatures.txt  > $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt
	grep -vh -P "^__" $shdir/htseq-count_$f.miRBase_mature_mergedFeatures.txt $shdir/htseq-count_$f.miRBase_hairpin_mergedFeatures.txt $shdir/htseq-count_$f.miRNA_mergedFeatures.txt  > $shdir/htseq-count_$f.all_miRNA_mergedFeatures.txt
	grep -vh -P "^__" $shdir/htseq-count_$f.MT_mergedFeatures.txt $shdir/htseq-count_$f.Mt_rRNA_mergedFeatures.txt $shdir/htseq-count_$f.Mt_tRNA_mergedFeatures.txt $shdir/htseq-count_$f.other\(MT\)_mergedFeatures.txt  > $shdir/htseq-count_$f.all_MT_mergedFeatures.txt
	grep -vh -P "^__" $shdir/htseq-count_$f.YRNA\(misc_RNA\)_mergedFeatures.txt $shdir/htseq-count_$f.noYorPiwi\(misc_RNA\)_mergedFeatures.txt  > $shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt
	grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 1' > $shdir/htseq-count_$f.morethen1.txt
	grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 0' > $shdir/htseq-count_$f.morethen0.txt
	cat $shdir/priority_$f.lincRNA_mergedFeatures.txt $shdir/priority_$f.lncipedia_hc_mergedFeatures.txt $shdir/priority_$f.lncipedia_mergedFeatures.txt | sort | uniq > $shdir/priority_$f.all_lncRNA_mergedFeatures.txt
	cat $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt | gzip > $WD/htseq-count_$f.all_lncRNA_mergedFeatures.txt.gz
	cat $shdir/htseq-count_$f.protein_coding_mergedFeatures.txt | gzip > $WD/htseq-count_$f.protein_coding_mergedFeatures.txt.gz
# grep -A4 -Ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt $ff | gzip > $WD/$f.all_lncRNA.fq.gz
# grep -A4 -Ff $shdir/priority_$f.protein_coding_mergedFeatures.txt $ff | gzip > $WD/$f.protein_coding.fq.gz
	seqtk subseq $ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt | gzip > $WD/$f.all_lncRNA.fq.gz
	seqtk subseq $ff $shdir/priority_$f.protein_coding_mergedFeatures.txt | gzip > $WD/$f.protein_coding.fq.gz
	"""
}


/*
 * Step x. 
 */
process countOver {
	tag "$genome_file.baseName"
	input:
	file f from alignments_files
	file genome from genome_file
	file gtf from genome_file
	file DV from specie
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	if [ "$DV" == "homo_sapiens" ]
	then 
		type="miRBase_mature_mergedFeatures"
	htseq_opt="htseq-count -f sam -a 0 -s no"
	for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lincRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Ensembl_genes_mergedFeatures"}
	do
		ann=$DB/gtf_biotypes/$type.gtf
		if [ ! -e $ann ]; then ann="$DB/gtf_biotypes/$type.gff"; fi
# if [ "$type" == "miRBase" ]; then ann="$ann --type miRNA --idattr Name"; fi
# $htseq_opt --nonunique all -q $a1 $ann -o $a2/tmp.a.sam > $a2/htseq-nopriority_$f.$type.a.txt
			$htseq_opt -q $shdir/$f.sam $ann -o $shdir/tmp.sam > $shdir/htseq-nopriority_$f.$type.txt
			grep -v __no_feature $shdir/tmp.sam | grep -v __not_aligned  > $shdir/htseq_$type.sam
			rm $shdir/tmp.sam
			grep -v "^@" $shdir/htseq_$type.sam | awk '{print $1}' > $shdir/htseq_reads_list.$type.txt
			done
			else
				grep -v "^@" $shdir/$f.sam | awk '{print $1}' > $shdir/htseq_reads_list.$DV.txt
				fi
	"""
}


/*
 * Step x. 
 */
process send2blast { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $fasta $out/$f 1 4; filename=${fasta##*/} ; fa=$out/$f/${filename%.*}_random${tsize}.1.faTab
	tag "$genome_file.baseName"
input:
	file genome from genome_file
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	tsize="2000"
	filename=${1##*/}
	fa=faTab/${filename%.*}_random${tsize}.${3}.faTab
		mem=`expr 80000 + $4 \* 2000`
	mkdir -p "$2/faTab"
	fasta_formatter -i $1 -t | shuf | head -n $tsize > $2/$fa
	# cat $1 | tr '\n' '\t' | sed 's/\t>/\n>/g' | shuf | head -n $tsize | sort | tr '\t' '\n' > $2/$fa
		sleep 2
	bsub -M $mem -n $4 -R "rusage[mem=$mem,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c ; export BLASTDB=$HOME/data/db/blast ; R --vanilla < $HOME/bin/parse_blast_output.R $2 $fa $4 > $2/$fa.log"
	"""
}
# export BLASTDB=$HOME/data/db/blast ; R --vanilla < $HOME/bin/parse_blast_output.R /homes/pavelz/data/shared/Cristina_merged_43-91/analysis_long/homo_sapiens/38N/ /homes/pavelz/data/shared/Cristina_merged_43-91/analysis_long/homo_sapiens/38N/38N_random10kB.1.fasta 4

/*
 * Step x. 
 */
process send2kraken { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $ff $out/$f 1 4
tag "$genome_file.baseName"
input:
	file genome from genome_file
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	filename=${1##*/}
	fa=$2/${filename%.*}
	bsub -M 42000 -n $4 -R "rusage[mem=42000,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c; kraken2 --threads $4 --db $HOME/data/db/kraken2 --fasta-input $1 2> $fa.kraken2.log | cut -f2,3 > $fa.kraken2_4krona.txt ; ktImportTaxonomy $fa.kraken2_4krona.txt -o $fa.kraken2_report.html 2>&1 >> $fa.kraken2.log"
	"""
}

/*
 * Step x. 
 */
process send2sourmash { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $ff $out/$f 1 4
tag "$genome_file.baseName"
input:
	file genome from genome_file
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	filename=${ff##*/}
	fa=$out/${f##*/}/${filename%.*}
	# trim-low-abund.py -C 3 -Z 18 -V -M 2e9 -o ${fa}.abundtrim $ff
	bsub -M 42000 -n 4 -R "rusage[mem=42000,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c; sourmash compute --scaled 1000 -k 21 ${fa}.abundtrim --merge ${filename%.*} -o ${fa}-reads.sig; sourmash gather -k 21 ${fa}-reads.sig $HOME/data/db/sourmash/refseq-d2-k21.sbt.json -o ${fa}-reads.refseq21.out; sourmash gather -k 21 ${fa}-reads.sig $HOME/data/db/sourmash/genbank-d2-k21.sbt.json -o  ${fa}-reads.genbank21.out; sourmash compute --scaled 1000 -k 21 ${ff} --merge ${filename%.*} -o ${fa}.sig; sourmash gather -k 21 ${fa}.sig $HOME/data/db/sourmash/refseq-d2-k21.sbt.json -о ${fa}.refseq21.out > $fa.sourmash.log"
	"""
}

/*
 * Step x. 
 */
process runBowtie2 {
	tag "$genome_file.baseName"
	input:
	file genome from genome_file
	
	output:
	file 'genome.index*' into genome_index
	
	"""
	i=""
	ff="$i"
	ff="38preU.fq"
	ff="$1"
	fn="${ff##*/}"
	f=${fn%.*}
	# ff=data/$f.fq
	inFasta=""
	txtLog=$out/$f/stat.txt
	# if [ ! -e "$out/$f/Unmapped_$f.fa" ]; then 
		mkdir $out/$f
	#   echo > $txtLog
		if [ "${ff}" == "data/$f.fa" ]; then 
	#      ff=data/$f.fa
			inFasta="-f"
		else
	#   # ad3=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC # Illumina
	# ad3=ATCACCGACTGCCCATAGAGAG  # Ion Torrent
	# pr3=AGGCTGAGACTGCCAAGGCACACAGGGGATAGG
	# ad5=CCAAGGCG
	# cutadapt --quiet --quality-cutoff=10,10 -O 9  -e 0.1 -a $ad3             -o $out/$f/${f}_r1.fastq $ff
	# cutadapt --quiet -O 9  -e 0.1 -g $ad3                -o $out/$f/${f}_r2.fastq $out/$f/${f}_r1.fastq
	# cutadapt --quiet -O 13 -e 0   -g ^ATCACCGACTGCCCATAG -o $out/$f/${f}_r3.fastq $out/$f/${f}_r2.fastq
	# cutadapt --quiet -O 13 -e 0   -g ^ATCACCGACTGCC      -o $out/$f/${f}_r4.fastq $out/$f/${f}_r3.fastq
	# cutadapt --quiet -m 10 -M 42                         -o $out/$f/$f.fastq      $out/$f/${f}_r4.fastq
	# ff=$out/$f/$f.fastq
	# fastqc -o qc $ff
	# #     filename=${ff##*/}
	# #     fasta=$out/$f/${filename%.*}.fasta
			fasta=$out/$f/$f.fasta
			fastq_to_fasta -i $ff -o $fasta
	# send2kraken $fasta $out/$f 0 $core
	# # send2sourmash $fasta $out/$f 0 $core
	# # send2blast $fasta $out/$f 1 $core
	# # send2blast $fasta $out/$f 2 $core
			fi 
	# sortmerna --ref $rRNAdb --reads $ff --aligned $out/$f/$f\_with_rRNA2 --fastx --other $out/$f/$f\_without_rRNA2 -m 4096 --log -a $core
	# mv data/${i}.fq.fastq data/${i}.fq; fastqc -o qc data/${i}.fq
			
	# if [ ! -e "$out/$f/Unmapped_$f.fa" ]; then 
			shdir="$out/$f/ShortStack.Bowtie2_N1_Seed20_wGaps_trimm"
			shfile="$out/$f/$f.bam"
	# rm -rf "$shdir"
	# bowtie2opt="--time --end-to-end -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal"
	# bowtie2 $bowtie2opt -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam >> $txtLog  2>&1
	# # bowtie2opt="--time --local -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal --score-min G,1,10"
	# # bowtie2 $bowtie2opt -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam >> $txtLog  2>&1
	# $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile
	# $shortstack --readfile $shfile --genomefile $DB/$DV.fa --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 21 --ranmax 20 --keep_quals --inbam --nohp
	# # rm $shfile $shdir.sam
	# $samtools index $shdir/$f.bam
	# $samtools idxstats $shdir/$f.bam  > ${shdir}/mapped.txt
	# $samtools sort $shdir/$f.bam -o $shdir/$f.sam
	# echo "number of mapped reads" `awk '{s+=$3} END {print s}' ${shdir}/mapped.txt` >> $txtLog
	# echo "number of unmapped reads" `awk '{s+=$4} END {print s}' ${shdir}/mapped.txt` >> $txtLog
			mycount $shdir/$f.bam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
				awk '!/^\*\t/ {print $1 "\t" $3}' ${shdir}/mapped.txt > ${shdir}/htseq-count_$DV.txt
	# countOver $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
	#   fi
	# # fi
	# ff="$out/$f/Unmapped_$f.fa"
	# rm $ff
	# if [ "$inFasta" == "" ]; then
	#   fastq_to_fasta -i $out/$f/Unmapped_$f.fq -o $ff
	# else
	#   ln -rfs $out/$f/Unmapped_$f.fq $ff
	# fi
	# send2kraken $ff $out/$f 0 $core
	# send2blast $ff $out/$f 1 $core
	# send2blast $ff $out/$f 2 $core
	# send2sourmash $fasta $out/$f 0 $core
	# # # rm -rf "$out/$f.1" "$out/$f.2"  # --blast '1 cigar qcov qstrand'
	# # sortmerna --ref $rRNAdb --reads $out/$f/Unmapped_$f.fq --sam --SQ --blast 0 --aligned $out/$f/Unmapped_$f\_with_rRNA --fastx --other $out/$f/Unmapped_$f\_without_rRNA -m 4096 --log -a $core
	# # fastq_to_fasta -i $out/$f/Unmapped_$f\_with_rRNA.fq -o $out/$f/Unmapped_$f\_with_rRNA.fasta
	# # $samtools view -uhS -F4 $out/$f/Unmapped_$f\_with_rRNA.sam | $samtools sort -@ $core - -o $out/$f/Unmapped_$f\_with_rRNA.bam
	# # $samtools index $out/$f/Unmapped_$f\_with_rRNA.bam
	# # $samtools idxstats $out/$f/Unmapped_$f\_with_rRNA.bam | sort -k3 -n -r > $out/$f/mapped_rRNA.txt
			
	"""
}


/*
 * STEP 8 - MultiQC
 */
process multiqc {
	publishDir "${params.outdir}/MultiQC", mode: 'copy'
	
	when:
	!params.skipQC && !params.skipMultiqc
	
	input:
	file multiqc_config from ch_multiqc_config
	file ('fastqc/*') from fastqc_results.collect()
	file ('trim_galore/*') from trimgalore_results.collect()
	file ('mirtrace/*') from mirtrace_results.collect()
	file ('samtools/*') from ch_sort_bam_flagstat_mqc.collect()
	file ('samtools_genome/*') from ch_genome_bam_flagstat_mqc.collect().ifEmpty([])
	file ('software_versions/*') from software_versions_yaml.collect()
	file workflow_summary from create_workflow_summary(summary)
	
	output:
	file "*multiqc_report.html" into multiqc_report
	file "*_data"
	
	script:
	rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
	rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
	"""
	multiqc . -f $rtitle $rfilename --config $multiqc_config -m bowtie1 -m samtools -m cutadapt -m fastqc -m custom_content
	"""
}


/*
 * STEP 9 - Output Description HTML
 */
process output_documentation {
	publishDir "${params.outdir}/pipeline_info", mode: 'copy'
	
	input:
	file output_docs from ch_output_docs
	
	output:
	file "results_description.html"
	
	script:
	"""
	markdown_to_html.r $output_docs results_description.html
	"""
}




/*
 * Completion e-mail notification
 */
workflow.onComplete {
	def subject = "[sRNAflow] Successful: $workflow.runName"
	if(!workflow.success){
		subject = "[sRNAflow] FAILED: $workflow.runName"
	}
	sendMail{
		to: 'you@gmail.com'
        from: ''
        subject: subject
        body: 'Hi, how are you!'
        text: subject
		attach '/some/file.txt', fileName: 'manuscript.txt'
		attach '/other/image.png', disposition: 'inline'
		if (workflow.success) {
			mqc_report = multiqc_report.getVal()
			if (mqc_report.getClass() == ArrayList){
				log.warn "[sRNAflow] Found multiple reports from process 'multiqc', will use only one"
				mqc_report = mqc_report[0]
			}
		}
	}
}


