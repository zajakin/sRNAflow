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

# strategy="successively"; specie="homo_sapiens"; f="07NS-v1"; rawfile="Cristina_202005/data/07NS-v1.fastq.gz"; ff="/home/pawel/Desktop/sRNAflow/www/results/Prostate/07NS-v1/07NS-v1.fastq.gz"; ftype="fq"; out="/home/pawel/Desktop/sRNAflow/www/results/Prostate"
function mycount {
  Rscript --vanilla $rsubread $shdir $shdir/$f.sam $DB/$DV$ext $shdir/htseq-count_$f.txt -
  stLog="$shdir/htseq_${f}__no_feature.txt $shdir/htseq_${f}__ambiguous.txt $shdir/htseq_${f}__too_low_aQual.txt $shdir/htseq_${f}__not_aligned.txt $shdir/htseq_${f}__alignment_not_unique.txt"
  rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt > /dev/null 2>&1
  touch $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
  type="miRBase_hairpin_mergedFeatures"
  for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","vault_RNA_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
  do
    pigz -cd $DB/gtf_biotypes/$type$ext | sed "s/gene_id \"/gene_id \"${type}_/" >> $shdir/tmp.gtf
    # (grep "^#" file.gtf; grep -v "^#" file.gtf | sort -k1,1 -k4,4n) > file_sorted.gtf
    Rscript --vanilla $rsubread $shdir $shdir/$f.sam $shdir/tmp.gtf $shdir/tmp.txt $shdir/priority.sam
    grep -v "XF:Z:__no_feature" $shdir/priority.sam | grep -v "XF:Z:__not_aligned" > $shdir/priority_$f.$type.sam
    echo "$type `grep -c 'XF:Z:__no_feature' $shdir/priority.sam`" >> $shdir/htseq_${f}__no_feature.txt
    echo "$type `grep -c 'XF:Z:__not_aligned' $shdir/priority.sam`" >> $shdir/htseq_${f}__not_aligned.txt
    echo "$type `grep -c 'XF:Z:__ambiguous' $shdir/priority.sam`" >> $shdir/htseq_${f}__ambiguous.txt
    echo "$type `grep -c 'XF:Z:__too_low_aQual' $shdir/priority.sam`" >> $shdir/htseq_${f}__too_low_aQual.txt
    echo "$type `grep -c 'XF:Z:__alignment_not_unique' $shdir/priority.sam`" >> $shdir/htseq_${f}__alignment_not_unique.txt
    grep -v "XF:Z:__no_feature" $shdir/priority.sam | grep -v "XF:Z:__not_aligned" | grep -v "^@" | gawk '{print $1}' > $shdir/tmpReads.txt
    grep -v -Ff $shdir/testedID.txt $shdir/tmp.txt > $shdir/htseq-count_$f.$type.txt
    grep -v -Ff $shdir/testedReads.txt $shdir/tmpReads.txt > $shdir/priority_$f.$type.txt
    cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | gawk '{print $1}' >> $shdir/testedID.txt
    cat $shdir/priority_$f.$type.txt >> $shdir/testedReads.txt
    cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" >> $shdir/htseq-count_$f.priority.txt
  done
  rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/tmpReads.txt $shdir/tmp.txt $shdir/priority.sam $shdir/testedReads.txt > /dev/null 2>&1
  grep -vh -P "^__" $shdir/htseq-count_$f.piRNAdb_mergedFeatures.txt $shdir/htseq-count_$f.piRNAbank_mergedFeatures.txt $shdir/htseq-count_$f.piRBase_mergedFeatures.txt  > $shdir/htseq-count_$f.all_piRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.lncRNA_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_hc_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_mergedFeatures.txt  > $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.miRBase_mature_mergedFeatures.txt $shdir/htseq-count_$f.miRBase_hairpin_mergedFeatures.txt $shdir/htseq-count_$f.miRNA_mergedFeatures.txt  > $shdir/htseq-count_$f.all_miRNA_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.MT_mergedFeatures.txt $shdir/htseq-count_$f.Mt_rRNA_mergedFeatures.txt $shdir/htseq-count_$f.Mt_tRNA_mergedFeatures.txt $shdir/htseq-count_$f.other\(MT\)_mergedFeatures.txt  > $shdir/htseq-count_$f.all_MT_mergedFeatures.txt
  grep -vh -P "^__" $shdir/htseq-count_$f.YRNA\(misc_RNA\)_mergedFeatures.txt $shdir/htseq-count_$f.noYorPiwi\(misc_RNA\)_mergedFeatures.txt  > $shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt
}

function mystat {
  printf 'all_Ensembl\t%i\n' `cat $shdir/htseq-count_$f.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  for type in {"all_miRNA_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","all_MT_mergedFeatures","all_piRNA_mergedFeatures","all_lncRNA_mergedFeatures","vault_RNA_mergedFeatures","misc_RNA_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
  do
    printf '%s\t%i\n' "$type" `cat $shdir/htseq-count_$f.$type.txt | grep -v "^__" | gawk '{s+=$2} END {print s}'`
  done
}

function countOver {
  for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vault_RNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures","miRBase_mature","miRBase_hairpin","miRNA","piRNAdb","piRNAbank","piRBase","YRNA(misc_RNA)","snRNA","snoRNA","tRF","tRNAhalves","GtRNAdb","rRNA","MT","Mt_rRNA","Mt_tRNA","other(MT)","vault_RNA","lncRNA","lncipedia_hc","lncipedia","noYorPiwi(misc_RNA)","protein_coding","processed_pseudogene","Other_types","RepeatMasker_tRNA","RepeatMasker_rRNA","RepeatMasker","Ensembl_genes"}
  do
    Rscript --vanilla $rsubread $shdir $shdir/$f.sam $DB/gtf_biotypes/$type$ext $shdir/htseq-nopriority_$f.$type.txt $shdir/tmp.sam
    grep -v __no_feature $shdir/tmp.sam | grep -v __not_aligned  > $shdir/htseq_$type.sam
    rm $shdir/tmp.sam
    grep -v "^@" $shdir/htseq_$type.sam | gawk '{print $1}' > $shdir/htseq_reads_list.$type.txt
    gawk -v type="$type" '/__ambiguous/{print "Lost_in_" type "\t" $2}' $shdir/htseq-nopriority_$f.$type.txt
  done
}

core=`nproc`
#export PATH=$HOME/bin:$HOME/conda/bin:$HOME/.local/bin:$PATH
DV="genomes"
DB="$out/$DV"
DFA="$(pwd)/www/db/genomes/${specie}.fa"
DBW="$(pwd)/www/db/genomes/bowtie2/${specie}/${specie}"
samtools="samtools "
gencore="$(pwd)/bin/gencore "
rsubread="$(pwd)/bin/sRNAflow_rsubread.R "
shortstack="$(pwd)/ShortStack/ShortStack "
isomiRSEA="$(pwd)/bin/isomiR-SEA "
taxonomy="ktImportTaxonomy -tax $(pwd)/www/db/taxonomy "

pushd $out
mkdir -p $out/qc

inFasta=""
[ "${ftype}" == "fa" ] && inFasta="-f"
tax=""
[ "${strategy}" == "successively" ] && tax="_tax_filtRC"
txtLog=$out/$f/stat${tax}.txt
cp -f $out/$f/logs/trimm.txt $txtLog

shdir="$out/$f/ShortStack"
shfile="$out/$f/$f.bam"
# rm -rf "$shdir"
dd=""
ext="_tax.gtf.gz"
# bowtie2opt="--time --end-to-end -p $core --mm $inFasta --no-unal --very-sensitive --no-1mm-upfront";						   #tax="_tax_no_1mm"
bowtie2opt="--time --end-to-end -p $core --mm $inFasta --no-unal --very-sensitive --no-1mm-upfront --score-min L,-1.15,-0.24"; #tax="_tax_filterEC"; tax="_tax_filtHard"; tax="_tax_selByBowtie_withQC"
if [ "${strategy}" == "successively" ]; then
  ext=".gtf.gz"
  if [ ! -d "$shdir" ]; then
	  bowtie2 $bowtie2opt -k 21 -x $DBW --un $out/$f/Unmapped_2main_$f.fq -U $ff -S $shdir.sam > $out/$f/logs/bowtie2main.log 2>&1
	  $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1
	  rm $shdir.sam
	  $shortstack --readfile $shfile --genomefile $DFA --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp > $out/$f/logs/Shortstack_2main.log 2>&1
	  rm $shfile
	  if [ `echo $f | grep -c dd` == 1 ]; then  # Reads with UMI
	    dd="_dd"
	    $gencore -i $shdir/$f.bam -h $shdir/${f}_dd.html -o $shdir/$f$dd.bam -r $DFA -s 1 > $out/$f/logs/gencore_main.log
	  fi
	  $samtools index $shdir/$f$dd.bam
	  $samtools idxstats $shdir/$f$dd.bam  > ${shdir}/mapped.txt
	  $samtools sort $shdir/$f$dd.bam -o $shdir/$f.sam
  fi
  [ `echo $f | grep -c dd` == 1 ] && cat $out/$f/logs/gencore_main.log >> $txtLog
  gawk '/reads; of these/{ print "Main_For_mapping\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/were unpaired; of these/{ print "Main_Unpaired\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/were unpaired; of these/{ print "Main_Unpaired_%\t" $2 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned 0 times/{ print "Main_Unmapped\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned 0 times/{ print "Main_Unmapped_%\t" $2 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned exactly 1 time/{ print "Main_Uniq_mapped\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned exactly 1 time/{ print "Main_Uniq_mapped_%\t" $2 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned >1 times/{ print "Main_Multimapped\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/aligned >1 times/{ print "Main_Multimapped_%\t" $2 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '/overall alignment rate/{ print "Main_Overall_alignment_rate_%\t" $1 }' $out/$f/logs/bowtie2main.log >> $txtLog
  gawk '{l+=$3; s+=$4 } END {print "Main_mapped_reads\t" l+s "\nMain_mapped_reads_less200\t" l "\nMain_mapped_reads_over200\t" s}' ${shdir}/mapped.txt >> $txtLog
fi
if   [ ! -d "${shdir}${tax}" ]; then
	if [ "${strategy}" == "successively" ]; then
		bowtie2 $bowtie2opt -k 201 -x $DB/$DV --un $out/$f/Unmapped${tax}.fq -U $out/$f/Unmapped_2main_$f.fq -S $shdir.sam > $out/$f/logs/bowtie2${tax}.log 2>&1
	else
		bowtie2 $bowtie2opt -k 201 -x $DB/$DV --un $out/$f/Unmapped${tax}.fq -U $ff                          -S $shdir.sam > $out/$f/logs/bowtie2${tax}.log 2>&1
	fi
	# $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1
	cat $shdir.sam | grep -v -F -f $out/genomeless/exclude | $samtools view -uhS -F4 | $samtools sort -@ $core - -o $shfile > /dev/null 2>&1  #  -w
	rm $shdir.sam
	rm -rf ${shdir}${tax}
	# mkdir -p ${shdir}${tax}; mv $shfile ${shdir}${tax}/$f.bam
	$shortstack --readfile $shfile --genomefile $DB/$DV.fa --outdir ${shdir}${tax} --bowtie_cores $core --mismatches 1 --bowtie_m 201 --ranmax 200 --keep_quals --inbam --nohp > $out/$f/logs/Shortstack${tax}.log 2>&1
	$samtools view -o $out/$f/$f.sam.gz $shfile && rm $shfile
	if [ `echo $f | grep -c dd` == 1 ]; then
	  dd="_dd"
####TODO Statistic from gencore
	  $gencore -i ${shdir}${tax}/$f.bam -h $shdir/${f}_dd.html -o ${shdir}${tax}/$f$dd.bam -r $DB/$DV.fa -s 1 > $out/$f/logs/gencore${tax}.log
	fi # if  UMI
	$samtools index ${shdir}${tax}/$f$dd.bam
	$samtools idxstats ${shdir}${tax}/$f$dd.bam  > ${shdir}${tax}/mapped.txt
	$samtools sort ${shdir}${tax}/$f$dd.bam -o ${shdir}${tax}/$f.sam
fi
[ `echo $f | grep -c dd` == 1 ] && cat $out/$f/logs/gencore${tax}.log >> $txtLog
gawk '/reads; of these/{ print "For_mapping\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/were unpaired; of these/{ print "Unpaired\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/were unpaired; of these/{ print "Unpaired_%\t" $2 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
echo -e "Unmapped(no_environment)\t$(sed -n '2~4p' $out/$f/Unmapped${tax}.fq | grep -v -c -F -f $out/genomeless/exclude)" >> $txtLog # -w
gawk '/aligned 0 times/{ print "Unmapped\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/aligned 0 times/{ print "Unmapped_%\t" $2 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/aligned exactly 1 time/{ print "Uniq_mapped\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/aligned exactly 1 time/{ print "Uniq_mapped_%\t" $2 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/aligned >1 times/{ print "Multimapped\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/aligned >1 times/{ print "Multimapped_%\t" $2 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '/overall alignment rate/{ print "Overall_alignment_rate_%\t" $1 }' $out/$f/logs/bowtie2${tax}.log >> $txtLog
gawk '{l+=$3; s+=$4 } END {print "mapped_reads\t" l+s "\nmapped_reads_less200\t" l "\nmapped_reads_over200\t" s}' ${shdir}${tax}/mapped.txt >> $txtLog
[ ! -e "$shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt" ] && mycount $shdir/$f.sam $shdir $DB/$DV.gtf.gz $DB $DV
mystat $shdir/$f.sam $shdir $DB/$DV.gtf.gz $DB $DV >> $txtLog
# tail -n 1 ${shdir}/htseq_${f}__no_feature.txt | gawk '{print "No_feature\t" $2}' >> $txtLog
gawk '!/^\*\t/ {print $1 "\t" $3}' ${shdir}/mapped.txt > ${shdir}/htseq-count_$DV.txt
# countOver $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog

if [ ! -e "$out/species_diagrams/$f.report${tax}.htm" ]; then
	mkdir -p $out/$f/forKrona
	gawk '!/^@/ {if($3=="*") $3=0; $3=gensub(/_.*/, "",1,$3); print $1 "\t" $3}' ${shdir}${tax}/$f.sam > $out/$f/forKrona/$f.forKrona${tax}.txt
	gawk '/^@/ {$1=gensub(/^@/,"",1,$1); print $1 "\t0"}' $out/$f/Unmapped${tax}.fq >> $out/$f/forKrona/$f.forKrona${tax}.txt
	gawk '{print $2}' $out/$f/forKrona/$f.forKrona${tax}.txt | sort -n | uniq -c > $out/$f/forKrona/$f.counts${tax}.txt
	# rm $shdir/$f.sam ${shdir}${tax}/$f.sam
	$taxonomy $out/$f/forKrona/$f.forKrona${tax}.txt -o $out/species_diagrams/$f.report${tax}.htm
fi

popd

#rm -rf $out/$f/isomiR-SEA
if [ ! -d "$out/$f/isomiR-SEA" ]; then
	mkdir -p $out/$f/isomiR-SEA
	cp www/db/genomes/mature.fa $out/$f/isomiR-SEA/mature.txt
	pushd $out/$f/isomiR-SEA
	(sed -n '1~4s/^@/>/p;2~4p' $out/$f/Unmapped${tax}.fq | gawk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$1);N++;next;} {printf("%s",$0);} END {printf("\n");}' | gawk '{print $NF}'; \
	  cat $out/$f/ShortStack/priority_$f.miR*.sam | gawk '{print $10}') \
	  | sort | uniq -c | sort -nr | gawk '{print $2 "\t" $1}' > $out/$f/isomiR-SEA/$f.txt
	$isomiRSEA -s hsa -l 16 -b 4 -i $out/$f/isomiR-SEA/ -p $out/$f/isomiR-SEA -ss 6 -h 11 -m mature -t $f > summary.txt
	popd
fi
