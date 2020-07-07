#!/bin/bash
function mycount {
      htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q" #--additional-attr=gene_name
      # htseq_opt="featureCounts -s no"
      $htseq_opt $shdir/$f.sam $DB/$DV.gtf > $shdir/htseq-count_$f.txt  #  -s reverse
      printf '%s\t%i\n' $DV `cat $shdir/htseq-count_$f.txt | grep -v "^__" | awk '{s+=$2} END {print s}'`
      stLog="$shdir/htseq_${f}__no_feature.txt $shdir/htseq_${f}__ambiguous.txt $shdir/htseq_${f}__too_low_aQual.txt $shdir/htseq_${f}__not_aligned.txt $shdir/htseq_${f}__alignment_not_unique.txt" 
      rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt > /dev/null 2>&1
      touch $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
      type="miRBase_hairpin_mergedFeatures"
      # for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
      for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","vaultRNA_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
          do
            cat $DB/gtf_biotypes/$type.gtf | sed "s/gene_id \"/gene_id \"${type}_/" >> $shdir/tmp.gtf
            # (grep "^#" file.gtf; grep -v "^#" file.gtf | sort -k1,1 -k4,4n) > file_sorted.gtf
            $htseq_opt $shdir/$f.sam $shdir/tmp.gtf --samout=$shdir/priority.sam > $shdir/tmp.txt  #  -s reverse
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
      grep -v -Ff $shdir/testedID.txt $shdir/htseq-count_$f.txt > $shdir/htseq-count_$f.Lost_reads.txt
      printf '%s\t%i\n' "Lost_reads" `cat $shdir/htseq-count_$f.Lost_reads.txt | grep -v "^__" | awk '{s+=$2} END {print s}'`
      cat $shdir/htseq-count_$f.Lost_reads.txt | awk '$2 > 0' > $shdir/htseq-count_$f.Lost_reads.morethen0.txt
      cat $shdir/htseq-count_$f.Lost_reads.morethen0.txt | awk '{print $1}' > $shdir/tmp.txt
      grep -Ff $shdir/tmp.txt $DB/$DV.gtf > $shdir/htseq-count_$f.Lost_reads.gtf
      mv $shdir/tmp.txt $shdir/htseq-count_$f.All.txt
      # rm $shdir/tmp.gtf
      # GTFmiRNA=`ls $GTFPATH*.gtf | grep "miRNA"`
      grep -vh -P "^__" $shdir/htseq-count_$f.piRNAdb_mergedFeatures.txt $shdir/htseq-count_$f.piRNAbank_mergedFeatures.txt $shdir/htseq-count_$f.piRBase_mergedFeatures.txt  > $shdir/htseq-count_$f.all_piRNA_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.lncRNA_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_hc_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_mergedFeatures.txt  > $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.miRBase_mature_mergedFeatures.txt $shdir/htseq-count_$f.miRBase_hairpin_mergedFeatures.txt $shdir/htseq-count_$f.miRNA_mergedFeatures.txt  > $shdir/htseq-count_$f.all_miRNA_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.MT_mergedFeatures.txt $shdir/htseq-count_$f.Mt_rRNA_mergedFeatures.txt $shdir/htseq-count_$f.Mt_tRNA_mergedFeatures.txt $shdir/htseq-count_$f.other\(MT\)_mergedFeatures.txt  > $shdir/htseq-count_$f.all_MT_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.YRNA\(misc_RNA\)_mergedFeatures.txt $shdir/htseq-count_$f.noYorPiwi\(misc_RNA\)_mergedFeatures.txt  > $shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt
      grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 1' > $shdir/htseq-count_$f.morethen1.txt
      grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 0' > $shdir/htseq-count_$f.morethen0.txt
      cat $shdir/priority_$f.lncRNA_mergedFeatures.txt $shdir/priority_$f.lncipedia_hc_mergedFeatures.txt $shdir/priority_$f.lncipedia_mergedFeatures.txt | sort | uniq > $shdir/priority_$f.all_lncRNA_mergedFeatures.txt
      cat $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt | gzip > $out/htseq-count_$f.all_lncRNA_mergedFeatures.txt.gz
      cat $shdir/htseq-count_$f.protein_coding_mergedFeatures.txt | gzip > $out/htseq-count_$f.protein_coding_mergedFeatures.txt.gz
      # grep -A4 -Ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt $ff | gzip > $out/$f.all_lncRNA.fq.gz
      # grep -A4 -Ff $shdir/priority_$f.protein_coding_mergedFeatures.txt $ff | gzip > $out/$f.protein_coding.fq.gz
      # seqtk subseq $ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt | gzip > $out/$f.all_lncRNA.fq.gz
      # seqtk subseq $ff $shdir/priority_$f.protein_coding_mergedFeatures.txt | gzip > $out/$f.protein_coding.fq.gz
}

function countOver {
  # if [ "$DV" == "homo_sapiens" ]
  # then 
    # type="miRBase_mature_mergedFeatures"
    htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q"
    for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
        do
          ann=$DB/gtf_biotypes/$type.gtf
          if [ ! -e $ann ]; then ann="$DB/gtf_biotypes/$type.gff"; fi
          # if [ "$type" == "miRBase" ]; then ann="$ann --type miRNA --idattr Name"; fi
          # $htseq_opt --nonunique all -q $a1 $ann -o $a2/tmp.a.sam > $a2/htseq-nopriority_$f.$type.a.txt
          $htseq_opt $shdir/$f.sam $ann -o $shdir/tmp.sam > $shdir/htseq-nopriority_$f.$type.txt
          grep -v __no_feature $shdir/tmp.sam | grep -v __not_aligned  > $shdir/htseq_$type.sam
          rm $shdir/tmp.sam
          grep -v "^@" $shdir/htseq_$type.sam | awk '{print $1}' > $shdir/htseq_reads_list.$type.txt
          awk -v type="$type" '/__ambiguous/{print "Lost_in_" type "\t" $2}' $shdir/htseq-nopriority_$f.$type.txt
    done
  # else
  #   grep -v "^@" $shdir/$f.sam | awk '{print $1}' > $shdir/htseq_reads_list.$DV.txt
  # fi
}

core=`nproc`
export PATH=$HOME/bin:$HOME/conda/bin:$HOME/.local/bin:$PATH
# origout=analysis_long.withRepeats
out="$2"
DV="genomes"
DB="$out/$DV"
samtools="samtools "
gencore="$HOME/conda/envs/gencore/bin/gencore "
shortstack="$(pwd)/bin/ShortStack "
rRNAdb="$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-id95.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-18s-id95:$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s-id98.fasta,$HOME/data/db/sortmerna/rRNA_databases/silva-euk-28s-id98:$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5s-database-id98-dna:$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-database-id98-dna.fasta,$HOME/data/db/sortmerna/rRNA_databases/rfam-5.8s-database-id98-dna"

cd $out
mkdir -p $out/qc $out/qc_raw

ff="$1"
fn="${ff##*/}"
f=${fn%.*}
# ff=data/$f.fq
inFasta=""
txtLog=$out/$f/stat.txt
# if [ ! -e "$out/$f/Unmapped_$f.fa" ]; then 
  # mkdir $out/$f
  echo -e "File\t$f" > $txtLog
  if [ "${ff}" == "data/$f.fa" ]; then 
#      ff=data/$f.fa
     inFasta="-f"
     fasta_formatter -i $ff -t > $out/$f/${f}.faTab
     echo -e "Raw\t`grep -c \"^>\" $ff`" >> $txtLog
     echo -e "QC_and_adapter3\t`grep -c \"^>\" $ff`" >> $txtLog
     gawk 'length($NF) < 43' $out/$f/${f}.faTab > $out/$f/${f}_r4.faTab
     echo -e "Filter_more42\t`wc -l $out/$f/${f}_r4.faTab`" >> $txtLog
     gawk 'length($NF) > 9 {print ">" $1 "\n" $NF}' $out/$f/${f}_r4.faTab > $out/$f/${f}.fasta
     echo -e "Filter_less10\t`grep -c \"^>\" $out/$f/${f}.fasta`" >> $txtLog
  else
    ad3=TGGAATTCTCGGGTGCCAAGG # Illumina TruSeq Small RNA
    ad5=GTTCAGAGTTCTACAGTCCGACGATC # Illumina TruSeq Small RNA
     # ad3=ATCACCGACTGCCCATAGAGAG  # Ion Torrent
     # pr3=AGGCTGAGACTGCCAAGGCACACAGGGGATAGG
     # ad5=CCAAGGCG
     fastqc -o $out/qc_raw $ff  > /dev/null
     echo -e "Raw\t`grep -c ^@ $ff`" >> $txtLog
     cutadapt --quiet --quality-cutoff=20,20 -a $ad3 -m 1 -o $out/$f/${f}_r3.fastq $ff
     echo -e "QC_and_adapter3\t`grep -c ^@ $out/$f/${f}_r3.fastq`" >> $txtLog
     # cutadapt --quiet --quality-cutoff=10,10 -O 9  -e 0.1 -a $ad3 -o $out/$f/${f}_r1.fastq $ff
     # cutadapt --quiet -O 9  -e 0.1 -g $ad3                -o $out/$f/${f}_r2.fastq $out/$f/${f}_r1.fastq
     # cutadapt --quiet -O 13 -e 0   -g XATCACCGACTGCCCATAG -o $out/$f/${f}_r3.fastq $out/$f/${f}_r2.fastq
     # cutadapt --quiet -O 8 -e 0   -g ^ATCACCGACTGCC      -o $out/$f/${f}_r4.fastq $out/$f/${f}_r3.fastq
     # cutadapt --quiet -m 43                         -o $out/$f/${f}_long.fastq      $out/$f/${f}_r3.fastq
     # cutadapt --quiet -M 42                         -o $out/$f/${f}_r4.fastq      $out/$f/${f}_r3.fastq
     # echo -e "Filter_more42\t`grep -c ^@ $out/$f/${f}_r4.fastq`" >> $txtLog
     cutadapt --quiet -m 10                         -o $out/$f/$f.fastq      $out/$f/${f}_r3.fastq
     echo -e "Filter_less10\t`grep -c ^@ $out/$f/${f}.fastq`" >> $txtLog
     ff=$out/$f/$f.fastq
     fastqc -o $out/qc $ff > /dev/null
    fasta=$out/$f/$f.fasta
    fastq_to_fasta -i $ff -o $fasta
    gawk '!/^>/{ i++; seq[i] = length($1); sum += length($1) } END { asort(seq,sorted,"@val_num_asc"); print("reads_after_trimm\t" i "\nmedian_reads_length\t" sorted[int(i/2)] "\nmean_reads_length\t" sum / i) }' $out/$f/$f.fasta >> $txtLog
  fi 

  # if [ ! -e "$out/$f/Unmapped_$f.fa" ]; then 
    shdir="$out/$f/ShortStack"
    shfile="$out/$f/$f.bam"
    rm -rf "$shdir"
    bowtie2opt="--time --end-to-end -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal"
    bowtie2 $bowtie2opt -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam > $out/$f/bowtie2.log 2>&1
    # # bowtie2opt="--time --local -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal --score-min G,1,10"
    # # bowtie2 $bowtie2opt -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam >> $txtLog  2>&1
    $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile
    $shortstack --readfile $shfile --genomefile $DB/$DV.fa --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp >> $out/$f/Shortstack.log 2>&1
    # # rm $shfile $shdir.sam
    dd=""
    if [ `echo $f | grep -c dd` == 1 ]; then
      dd="_dd"
      $gencore -i $shdir/$f.bam -h $shdir/${f}_dd.html -o $shdir/$f$dd.bam -r $DB/$DV.fa -s 1 >> $txtLog
    fi # if  UMI
    $samtools index $shdir/$f$dd.bam
    $samtools idxstats $shdir/$f$dd.bam  > ${shdir}/mapped.txt
    $samtools sort $shdir/$f$dd.bam -o $shdir/$f.sam


    gawk '/reads; of these/{ print "For_mapping\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/were unpaired; of these/{ print "Unpaired\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/were unpaired; of these/{ print "Unpaired_%\t" $2 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned 0 times/{ print "Unmapped\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned 0 times/{ print "Unmapped_%\t" $2 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned exactly 1 time/{ print "Uniq_mapped\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned exactly 1 time/{ print "Uniq_mapped_%\t" $2 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned >1 times/{ print "Multimapped\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/aligned >1 times/{ print "Multimapped_%\t" $2 }' $out/$f/bowtie2.log >> $txtLog
    gawk '/overall alignment rate/{ print "Overall_alignment_rate_%\t" $1 }' $out/$f/bowtie2.log >> $txtLog
    awk '{l+=$3; s+=$4 } END {print "mapped_reads\t" l+s "\nmapped_reads_less20\t" l "\nmapped_reads_over20\t" s}' ${shdir}/mapped.txt >> $txtLog
    mycount $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
    tail -n 1 ${shdir}/htseq_${f}__no_feature.txt | awk '{print "No_feature\t" $2}' >> $txtLog
    awk '!/^\*\t/ {print $1 "\t" $3}' ${shdir}/mapped.txt > ${shdir}/htseq-count_$DV.txt
    # countOver $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
#   fi
# # fi
    awk '!/^@/ {if($3=="*") $3=0; $3=gensub(/_.*/, "",1,$3); print $1 "\t" $3}' $shdir/$f.sam > $out/$f/forKrona/$f.forKrona.txt
    awk '{print $2}' $out/$f/forKrona/$f.forKrona.txt | sort -n | uniq -c > $out/$f/forKrona/$f.counts.txt
    $HOME/conda/bin/ktImportTaxonomy $out/$f/forKrona/$f.forKrona.txt -o $out/$f/output/$f.report.htm
    $HOME/conda/bin/ktImportTaxonomy -c $out/$f/forKrona/$f.forKrona2.txt -o $out/$f/output/$f.report3.htm


