#!/bin/bash
while getopts ":s:v:n:r:f:t:o:" opt; do
  case $opt in
    s) strategy="$OPTARG"
    ;;
    v) specie="$OPTARG"
    ;;
    n) f="$OPTARG"
    ;;
    r) rawfile="$OPTARG"
    ;;
    f) ff="$OPTARG"
    ;;
    t) ftype="$OPTARG"
    ;;
    o) out="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

function mycount {
  htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q" #--additional-attr=gene_name
  # htseq_opt="featureCounts -s no"
  $htseq_opt $shdir/$f.sam $DB/$DV.gtf > $shdir/htseq-count_$f.txt  #  -s reverse
  # printf '%s\t%i\n' $DV `cat $shdir/htseq-count_$f.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  stLog="$shdir/htseq_${f}__no_feature.txt $shdir/htseq_${f}__ambiguous.txt $shdir/htseq_${f}__too_low_aQual.txt $shdir/htseq_${f}__not_aligned.txt $shdir/htseq_${f}__alignment_not_unique.txt"
  rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt > /dev/null 2>&1
  touch $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
  type="miRBase_hairpin_mergedFeatures"
  for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","vaultRNA_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
  do
    pigz -cd $DB/gtf_biotypes/$type.gtf.gz | sed "s/gene_id \"/gene_id \"${type}_/" >> $shdir/tmp.gtf
    # (grep "^#" file.gtf; grep -v "^#" file.gtf | sort -k1,1 -k4,4n) > file_sorted.gtf
    $htseq_opt $shdir/$f.sam $shdir/tmp.gtf --samout=$shdir/priority.sam > $shdir/tmp.txt  #  -s reverse
    grep -v "XF:Z:__no_feature" $shdir/priority.sam | grep -v "XF:Z:__not_aligned" > $shdir/priority_$f.$type.sam
    grep -v "XF:Z:__no_feature" $shdir/priority.sam | grep -v "XF:Z:__not_aligned" | grep -v "^@" | gawk '{print $1}' > $shdir/tmpReads.txt
    # egrep -v -f $shdir/testedID.txt $shdir/tmp.txt > $shdir/htseq-count_$f.$type.txt # | sed '$ s/\n//' | tr '\n' '|'
    grep -v -Ff $shdir/testedID.txt $shdir/tmp.txt > $shdir/htseq-count_$f.$type.txt
    grep -v -Ff $shdir/testedReads.txt $shdir/tmpReads.txt > $shdir/priority_$f.$type.txt
    echo "$type `grep __no_feature $shdir/tmp.txt  | gawk '{print $2}'`" >> $shdir/htseq_${f}__no_feature.txt
    echo "$type `grep __ambiguous  $shdir/tmp.txt  | gawk '{print $2}'`" >> $shdir/htseq_${f}__ambiguous.txt
    echo "$type `grep __too_low_aQual $shdir/tmp.txt  | gawk '{print $2}'`" >> $shdir/htseq_${f}__too_low_aQual.txt
    echo "$type `grep __not_aligned $shdir/tmp.txt  | gawk '{print $2}'`" >> $shdir/htseq_${f}__not_aligned.txt
    echo "$type `grep __alignment_not_unique $shdir/tmp.txt  | gawk '{print $2}'`" >> $shdir/htseq_${f}__alignment_not_unique.txt
    # printf '%s\t%i\n' "$type" `cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
    cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | gawk '{print $1}' >> $shdir/testedID.txt
    cat $shdir/priority_$f.$type.txt >> $shdir/testedReads.txt
    cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" >> $shdir/htseq-count_$f.priority.txt
  done
  # grep -v -Ff $shdir/testedID.txt $shdir/htseq-count_$f.txt > $shdir/htseq-count_$f.Lost_reads.txt
  # cat $shdir/htseq-count_$f.Lost_reads.txt | gawk '$2 > 0' > $shdir/htseq-count_$f.Lost_reads.morethen0.txt
  # cat $shdir/htseq-count_$f.Lost_reads.morethen0.txt | gawk '{print $1}' > $shdir/tmp.txt
  # grep -Ff $shdir/tmp.txt $DB/$DV.gtf > $shdir/htseq-count_$f.Lost_reads.gtf
  # mv $shdir/tmp.txt $shdir/htseq-count_$f.All.txt
  rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/tmpReads.txt $shdir/tmp.txt $shdir/priority.sam $shdir/testedReads.txt > /dev/null 2>&1
  # GTFmiRNA=`ls $GTFPATH*.gtf | grep "miRNA"`
  grep -vh -P "^__" $shdir/htseq-count_$f.piRNAdb_mergedFeatures.txt $shdir/htseq-count_$f.piRNAbank_mergedFeatures.txt $shdir/htseq-count_$f.piRBase_mergedFeatures.txt  > $shdir/htseq-count_$f.all_piRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.lncRNA_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_hc_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_mergedFeatures.txt  > $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.miRBase_mature_mergedFeatures.txt $shdir/htseq-count_$f.miRBase_hairpin_mergedFeatures.txt $shdir/htseq-count_$f.miRNA_mergedFeatures.txt  > $shdir/htseq-count_$f.all_miRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.MT_mergedFeatures.txt $shdir/htseq-count_$f.Mt_rRNA_mergedFeatures.txt $shdir/htseq-count_$f.Mt_tRNA_mergedFeatures.txt $shdir/htseq-count_$f.other\(MT\)_mergedFeatures.txt  > $shdir/htseq-count_$f.all_MT_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.YRNA\(misc_RNA\)_mergedFeatures.txt $shdir/htseq-count_$f.noYorPiwi\(misc_RNA\)_mergedFeatures.txt  > $shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt
  # grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | gawk '$2 > 1' > $shdir/htseq-count_$f.morethen1.txt
  # grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | gawk '$2 > 0' > $shdir/htseq-count_$f.morethen0.txt
  # cat $shdir/priority_$f.lncRNA_mergedFeatures.txt $shdir/priority_$f.lncipedia_hc_mergedFeatures.txt $shdir/priority_$f.lncipedia_mergedFeatures.txt | sort | uniq > $shdir/priority_$f.all_lncRNA_mergedFeatures.txt
}

function mystat {
  printf 'all_Ensembl\t%i\n' `cat $shdir/htseq-count_$f.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  for type in {"all_miRNA_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","all_MT_mergedFeatures","all_piRNA_mergedFeatures","all_lncRNA_mergedFeatures","vaultRNA_mergedFeatures","misc_RNA_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
  do
    printf '%s\t%i\n' "$type" `cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  done
  # printf '%s\t%i\n' "Lost_reads" `cat $shdir/htseq-count_$f.Lost_reads.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  # cat $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt | pigz > $out/htseq-count_$f.all_lncRNA_mergedFeatures.txt.gz
  # cat $shdir/htseq-count_$f.protein_coding_mergedFeatures.txt | pigz > $out/htseq-count_$f.protein_coding_mergedFeatures.txt.gz
  # grep -A4 -Ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt $ff | pigz > $out/$f.all_lncRNA.fq.gz
  # grep -A4 -Ff $shdir/priority_$f.protein_coding_mergedFeatures.txt $ff | pigz > $out/$f.protein_coding.fq.gz
  # seqtk subseq $ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt | pigz > $out/$f.all_lncRNA.fq.gz
  # seqtk subseq $ff $shdir/priority_$f.protein_coding_mergedFeatures.txt | pigz > $out/$f.protein_coding.fq.gz
}

function countOver {
  htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q"
  for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures","miRBase_mature","miRBase_hairpin","miRNA","piRNAdb","piRNAbank","piRBase","YRNA(misc_RNA)","snRNA","snoRNA","tRF","tRNAhalves","GtRNAdb","rRNA","MT","Mt_rRNA","Mt_tRNA","other(MT)","vaultRNA","lncRNA","lncipedia_hc","lncipedia","noYorPiwi(misc_RNA)","protein_coding","processed_pseudogene","Other_types","RepeatMasker_tRNA","RepeatMasker_rRNA","RepeatMasker","Ensembl_genes"}
  do
    ann=$DB/gtf_biotypes/$type.gtf.gz
    if [ ! -e $ann ]; then ann="$DB/gtf_biotypes/$type.gtf"; fi
    # $htseq_opt --nonunique all -q $a1 $ann -o $a2/tmp.a.sam > $a2/htseq-nopriority_$f.$type.a.txt
    $htseq_opt $shdir/$f.sam $ann -o $shdir/tmp.sam > $shdir/htseq-nopriority_$f.$type.txt
    grep -v __no_feature $shdir/tmp.sam | grep -v __not_aligned  > $shdir/htseq_$type.sam
    rm $shdir/tmp.sam
    grep -v "^@" $shdir/htseq_$type.sam | gawk '{print $1}' > $shdir/htseq_reads_list.$type.txt
    gawk -v type="$type" '/__ambiguous/{print "Lost_in_" type "\t" $2}' $shdir/htseq-nopriority_$f.$type.txt
  done
}

core=`nproc`
export PATH=$HOME/bin:$HOME/conda/bin:$HOME/.local/bin:$PATH
DV="genomes"
DB="$out/$DV"
DFA="$(pwd)/www/db/genomes/${specie}.fa"
DBW="$(pwd)/www/db/genomes/bowtie/${specie}/${specie}"
samtools="samtools "
gencore="$(pwd)/bin/gencore "
shortstack="$(pwd)/ShortStack/ShortStack "
isomiRSEA="$(pwd)/bin/isomiR-SEA "
taxonomy="ktImportTaxonomy -tax $(pwd)/www/db/taxonomy "

pushd $out
mkdir -p $out/qc # $out/qc_raw

if [ "${ftype}" == "fa" ]; then
  inFasta="-f"
  # mark="^>"
  # adremoved=$out/$f/${f}_r3.fasta
  # more42=$out/$f/${f}_r4.fasta
  # fasta=$ff
else
  inFasta=""
  # mark="^@"
  # adremoved=$out/$f/${f}_r3.fastq
  # more42=$out/$f/${f}_r4.fastq
  # fastqc -o $out/qc $ff > /dev/null 2>&1
  # fasta=$out/$f/$f.fasta
  # sed -n '1~4s/^@/>/p;2~4p' $ff > $fasta
fi
txtLog=$out/$f/stat.txt
cp -f $out/$f/logs/trimm.txt $txtLog
# echo -e "File\t$f" > $txtLog
# echo -e "Raw\t`grep -c \"$mark\" $rawfile`" >> $txtLog
# echo -e "QC_and_adapter3\t`grep -c \"$mark\" $adremoved`" >> $txtLog
# echo -e "Filter_long\t`grep -c \"$mark\" $more42`" >> $txtLog
# echo -e "Filter_short\t`grep -c \"$mark\" $ff`" >> $txtLog
# gawk '!/^>/{ i++; seq[i] = length($1); sum += length($1) } END { asort(seq,sorted,"@val_num_asc"); print("reads_after_trimm\t" i "\nmedian_reads_length\t" sorted[int(i/2)] "\nmean_reads_length\t" sum / i) }' $fasta >> $txtLog

shdir="$out/$f/ShortStack"
shfile="$out/$f/$f.bam"
rm -rf "$shdir"
dd=""
bowtie2opt="--time --end-to-end -p $core --mm $inFasta --no-unal -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"
if [ "${strategy}" == "successively" ]; then
  bowtie2 $bowtie2opt -k 21 -x $DBW --un $out/$f/Unmapped_2main_$f.fq -U $ff -S $shdir.sam > $out/$f/logs/bowtie2main.log 2>&1
  $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1
  rm $shdir.sam
  $shortstack --readfile $shfile --genomefile $DFA --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp >> $out/$f/logs/Shortstack_2main.log 2>&1
  rm $shfile
  if [ `echo $f | grep -c dd` == 1 ]; then
    dd="_dd"
    $gencore -i $shdir/$f.bam -h $shdir/${f}_dd.html -o $shdir/$f$dd.bam -r $DFA -s 1 >> $txtLog
  fi # if  UMI
  $samtools index $shdir/$f$dd.bam
  $samtools idxstats $shdir/$f$dd.bam  > ${shdir}/mapped.txt
  $samtools sort $shdir/$f$dd.bam -o $shdir/$f.sam
  
  gawk '/reads; of these/{ print "For_mapping\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/were unpaired; of these/{ print "Unpaired\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/were unpaired; of these/{ print "Unpaired_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned 0 times/{ print "Unmapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned 0 times/{ print "Unmapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned exactly 1 time/{ print "Uniq_mapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned exactly 1 time/{ print "Uniq_mapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned >1 times/{ print "Multimapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/aligned >1 times/{ print "Multimapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '/overall alignment rate/{ print "Overall_alignment_rate_%\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
  gawk '{l+=$3; s+=$4 } END {print "mapped_reads\t" l+s "\nmapped_reads_less200\t" l "\nmapped_reads_over200\t" s}' ${shdir}/mapped.txt >> $txtLog

  bowtie2 $bowtie2opt -k 201 -x $DB/$DV --un $out/$f/Unmapped_$f.fq -U $out/$f/Unmapped_2main_$f.fq -S $shdir.sam > $out/$f/logs/bowtie2.log 2>&1
  $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1
  rm $shdir.sam
  $shortstack --readfile $shfile --genomefile $DFA --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp >> $out/$f/logs/Shortstack.log 2>&1
  $samtools view -o $out/$f/$f.sam.gz $shfile
else
  bowtie2 $bowtie2opt -k 201 -x $DB/$DV --un $out/$f/Unmapped_$f.fq -U $ff -S $shdir.sam > $out/$f/logs/bowtie2.log 2>&1
  $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1
  rm $shdir.sam
  $shortstack --readfile $shfile --genomefile $DB/$DV.fa --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp >> $out/$f/logs/Shortstack.log 2>&1
  $samtools view -o $out/$f/$f.sam.gz $shfile
fi
rm $shfile
if [ `echo $f | grep -c dd` == 1 ]; then
  dd="_dd"
  $gencore -i $shdir/$f.bam -h $shdir/${f}_dd.html -o $shdir/$f$dd.bam -r $DB/$DV.fa -s 1 >> $txtLog
fi # if  UMI
$samtools index $shdir/$f$dd.bam
$samtools idxstats $shdir/$f$dd.bam  > ${shdir}/mapped.txt
$samtools sort $shdir/$f$dd.bam -o $shdir/$f.sam

gawk '/reads; of these/{ print "For_mapping\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/were unpaired; of these/{ print "Unpaired\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/were unpaired; of these/{ print "Unpaired_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned 0 times/{ print "Unmapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned 0 times/{ print "Unmapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned exactly 1 time/{ print "Uniq_mapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned exactly 1 time/{ print "Uniq_mapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned >1 times/{ print "Multimapped\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/aligned >1 times/{ print "Multimapped_%\t" $2 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '/overall alignment rate/{ print "Overall_alignment_rate_%\t" $1 }' $out/$f/logs/bowtie2.log >> $txtLog
gawk '{l+=$3; s+=$4 } END {print "mapped_reads\t" l+s "\nmapped_reads_less200\t" l "\nmapped_reads_over200\t" s}' ${shdir}/mapped.txt >> $txtLog
mycount $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV
mystat $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
tail -n 1 ${shdir}/htseq_${f}__no_feature.txt | gawk '{print "No_feature\t" $2}' >> $txtLog
# gawk '!/^\*\t/ {print $1 "\t" $3}' ${shdir}/mapped.txt > ${shdir}/htseq-count_$DV.txt
# countOver $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
gawk '!/^@/ {if($3=="*") $3=0; $3=gensub(/_.*/, "",1,$3); print $1 "\t" $3}' $shdir/$f.sam > $out/$f/forKrona/$f.forKrona.txt
gawk '/^@/ {$1=gensub(/^@/,"",1,$1); print $1 "\t0"}' $out/$f/Unmapped_$f.fq >> $out/$f/forKrona/$f.forKrona.txt

gawk '{print $2}' $out/$f/forKrona/$f.forKrona.txt | sort -n | uniq -c > $out/$f/forKrona/$f.counts.txt
# rm $shdir/$f.sam
$taxonomy $out/$f/forKrona/$f.forKrona.txt -o $out/species_diagrams/$f.report.htm
popd

rm -rf $out/$f/isomiR-SEA
mkdir -p $out/$f/isomiR-SEA
cp www/db/genomes/mature.fa $out/$f/isomiR-SEA/mature.txt
pushd $out/$f/isomiR-SEA
(sed -n '1~4s/^@/>/p;2~4p' $out/$f/Unmapped_$f.fq | gawk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$1);N++;next;} {printf("%s",$0);} END {printf("\n");}' | gawk '{print $NF}'; \
  cat $out/$f/ShortStack/priority_$f.miR*.sam | gawk '{print $10}') \
  | sort | uniq -c | sort -nr | gawk '{print $2 "\t" $1}' > $out/$f/isomiR-SEA/$f.txt
$isomiRSEA -s hsa -l 16 -b 4 -i $out/$f/isomiR-SEA/ -p $out/$f/isomiR-SEA -ss 6 -h 11 -m mature -t $f > summary.txt
popd

# rm -rf $out/$f/miTAR
# mkdir -p $out/$f/miTAR
# (cat $out/$f/ShortStack/priority_$f.miR*.sam | \
#   gawk '{$2=gensub(".*_hsa-","hsa-",1,$NF); $2=gensub(/\[.*[-+.]\]/,"","g",$2); $2=gensub("]$","",1,$2); $2=gensub(/\+.*/,"",1,$2); $2=gensub(".*_mergedFeatures_","",1,$2); print ">"$2"\n"$10}'; \
#   cat $out/$f/isomiR-SEA/out_result_mature_tag_unique.txt | gawk -F'\t' '!/^tag_index/{$2=gensub("U","T","g",$2); print $6"\n"$2}') | \
#   gawk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$1);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sort | uniq | gawk '{print $1"\n"$2}' > $out/$f/miTAR/$f.fa
# # export PYTHONPATH=~/.local/lib/python3.8/site-packages:/usr/local/lib/python3.8/dist-packages:$PYTHONPATH
# # sudo pip install --upgrade numpy
# # sed -i 's/tf.set_random_seed(sdnum)/tf.random.set_seed(sdnum)/p' $(pwd)/miTAR/predict_multimiRmultimRNA.py
# # /usr/bin/python3 $(pwd)/miTAR/predict_multimiRmultimRNA.py -i1 $out/$f/miTAR/$f.fa -i2 $DB/$DV.fa -o $out/$f/miTAR/$f\_predictedTar.fa -s 22 -p 0.8 -ns 1
# /usr/bin/python3 $(pwd)/miTAR/predict_multimiRmultimRNA.py -i1 $out/$f/miTAR/$f.fa -i2 $(pwd)/www/db/genomes/homo_sapiens.fa -o $out/$f/miTAR/$f\_predictedTar.fa -s 22 -p 0.8 -ns 1
# $(pwd)/bin/miranda $out/$f/miTAR/$f.fa www/db/genomes/homo_sapiens.fa -out $out/$f/miTAR/$f\_predicted_miRanda.txt
# head $(pwd)/www/db/genomes/homo_sapiens.fa
rm $shdir/*.sam $shdir/mapped.txt
