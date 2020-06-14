#!/bin/bash
function mycount {
      htseq_opt="htseq-count -f sam -a 0 -s no --secondary-alignments score -q" #--additional-attr=gene_name
      # htseq_opt="featureCounts -s no"
      $htseq_opt $shdir/$f.sam $DB/$DV.gtf > $shdir/htseq-count_$f.txt  #  -s reverse
      printf '%s\t%i\n' $DV `cat $shdir/htseq-count_$f.txt | grep -v "^__" | awk '{s+=$2} END {print s}'`
      stLog="$shdir/htseq_${f}__no_feature.txt $shdir/htseq_${f}__ambiguous.txt $shdir/htseq_${f}__too_low_aQual.txt $shdir/htseq_${f}__not_aligned.txt $shdir/htseq_${f}__alignment_not_unique.txt" 
      rm $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
      touch $shdir/testedID.txt $shdir/tmp.gtf $shdir/htseq-count_$f.priority.txt $stLog $shdir/testedReads.txt
      type="miRBase_hairpin_mergedFeatures"
      # for type in {"miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","miRBase_mature_mergedFeatures","GtRNAdb_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lincRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Ensembl_genes_mergedFeatures"} # ,"tRF","tRNAhalves"
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
      grep -vh -P "^__" $shdir/htseq-count_$f.lincRNA_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_hc_mergedFeatures.txt $shdir/htseq-count_$f.lncipedia_mergedFeatures.txt  > $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.miRBase_mature_mergedFeatures.txt $shdir/htseq-count_$f.miRBase_hairpin_mergedFeatures.txt $shdir/htseq-count_$f.miRNA_mergedFeatures.txt  > $shdir/htseq-count_$f.all_miRNA_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.MT_mergedFeatures.txt $shdir/htseq-count_$f.Mt_rRNA_mergedFeatures.txt $shdir/htseq-count_$f.Mt_tRNA_mergedFeatures.txt $shdir/htseq-count_$f.other\(MT\)_mergedFeatures.txt  > $shdir/htseq-count_$f.all_MT_mergedFeatures.txt
      grep -vh -P "^__" $shdir/htseq-count_$f.YRNA\(misc_RNA\)_mergedFeatures.txt $shdir/htseq-count_$f.noYorPiwi\(misc_RNA\)_mergedFeatures.txt  > $shdir/htseq-count_$f.misc_RNA_mergedFeatures.txt
      grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 1' > $shdir/htseq-count_$f.morethen1.txt
      grep -v -P "^__" $shdir/htseq-count_$f.txt | sort -nrk2 | awk '$2 > 0' > $shdir/htseq-count_$f.morethen0.txt
      cat $shdir/priority_$f.lincRNA_mergedFeatures.txt $shdir/priority_$f.lncipedia_hc_mergedFeatures.txt $shdir/priority_$f.lncipedia_mergedFeatures.txt | sort | uniq > $shdir/priority_$f.all_lncRNA_mergedFeatures.txt
      cat $shdir/htseq-count_$f.all_lncRNA_mergedFeatures.txt | gzip > $out/htseq-count_$f.all_lncRNA_mergedFeatures.txt.gz
      cat $shdir/htseq-count_$f.protein_coding_mergedFeatures.txt | gzip > $out/htseq-count_$f.protein_coding_mergedFeatures.txt.gz
      # grep -A4 -Ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt $ff | gzip > $out/$f.all_lncRNA.fq.gz
      # grep -A4 -Ff $shdir/priority_$f.protein_coding_mergedFeatures.txt $ff | gzip > $out/$f.protein_coding.fq.gz
      seqtk subseq $ff $shdir/priority_$f.all_lncRNA_mergedFeatures.txt | gzip > $out/$f.all_lncRNA.fq.gz
      seqtk subseq $ff $shdir/priority_$f.protein_coding_mergedFeatures.txt | gzip > $out/$f.protein_coding.fq.gz
}

function countOver {
  if [ "$DV" == "homo_sapiens" ]
  then 
    type="miRBase_mature_mergedFeatures"
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
  else
    grep -v "^@" $shdir/$f.sam | awk '{print $1}' > $shdir/htseq_reads_list.$DV.txt
  fi
}

function send2blast { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $fasta $out/$f 1 4; filename=${fasta##*/} ; fa=$out/$f/${filename%.*}_random${tsize}.1.faTab
  tsize="2000"
  filename=${1##*/}
  fa=faTab/${filename%.*}_random${tsize}.${3}.faTab
  mem=`expr 80000 + $4 \* 2000`
  mkdir -p "$2/faTab"
  fasta_formatter -i $1 -t | shuf | head -n $tsize > $2/$fa
  # cat $1 | tr '\n' '\t' | sed 's/\t>/\n>/g' | shuf | head -n $tsize | sort | tr '\t' '\n' > $2/$fa
  sleep 2
  bsub -M $mem -n $4 -R "rusage[mem=$mem,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c ; export BLASTDB=$HOME/data/db/blast ; R --vanilla < $HOME/bin/parse_blast_output.R $2 $fa $4 > $2/$fa.log"
}
# export BLASTDB=$HOME/data/db/blast ; R --vanilla < $HOME/bin/parse_blast_output.R /homes/pavelz/data/shared/Cristina_merged_43-91/analysis_long/homo_sapiens/38N/ /homes/pavelz/data/shared/Cristina_merged_43-91/analysis_long/homo_sapiens/38N/38N_random10kB.1.fasta 4
function send2kraken { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $ff $out/$f 1 4
  filename=${1##*/}
  fa=$2/${filename%.*}
  # mem=64000
  kraken2 --threads $4 --db $HOME/data/db/kraken2 --report $fa.kraken2.report --fasta-input $1 2> $fa.kraken2.log | cut -f2,3 > $fa.kraken2_4krona.txt
  ktImportTaxonomy $fa.kraken2_4krona.txt -o $fa.kraken2_report.html 2>&1 >> $fa.kraken2.log
  kraken-biom --max C --min G --fmt json -o $fa.json.biom $fa.kraken2.report
  # bsub -M $mem -n $4 -R "rusage[mem=$mem,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c; kraken2 --threads $4 --db $HOME/data/db/kraken2 --report $fa.kraken2.report --fasta-input $1 2> $fa.kraken2.log | cut -f2,3 > $fa.kraken2_4krona.txt ; ktImportTaxonomy $fa.kraken2_4krona.txt -o $fa.kraken2_report.html 2>&1 >> $fa.kraken2.log ; kraken-biom --max C --min G --fmt json -o $fa.json.biom $fa.kraken2.report"
}
function send2sourmash { #  set 1=$ff; set 2=$out/$f; set 3="1"; set 4=$core; send2blast $ff $out/$f 1 4
  filename=${ff##*/}
  fa=$out/${f##*/}/${filename%.*}
# trim-low-abund.py -C 3 -Z 18 -V -M 2e9 -o ${fa}.abundtrim $ff
  mem=64000
  bsub -M $mem -n 4 -R "rusage[mem=$mem,numcpus=$4.00] span[ptile=$4]" "$HOME/bin/c; sourmash compute --scaled 1000 -k 21 ${fa}.abundtrim --merge ${filename%.*} -o ${fa}-reads.sig; sourmash gather -k 21 ${fa}-reads.sig $HOME/data/db/sourmash/refseq-d2-k21.sbt.json -o ${fa}-reads.refseq21.out; sourmash gather -k 21 ${fa}-reads.sig $HOME/data/db/sourmash/genbank-d2-k21.sbt.json -o  ${fa}-reads.genbank21.out; sourmash compute --scaled 1000 -k 21 ${ff} --merge ${filename%.*} -o ${fa}.sig; sourmash gather -k 21 ${fa}.sig $HOME/data/db/sourmash/refseq-d2-k21.sbt.json -о ${fa}.refseq21.out > $fa.sourmash.log"
}
# for ff in ${files[@]}; do
#   f=${ff%.*}
#   send2sourmash $ff $out/${f##*/} 0 $core
# done

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

function send2cluster {
  # rm -rf $out
  # unset files
  # for i in {"38","56","60","215"}; do  #i="215"
  #   for j in {"preU","postU","preP","postP","N","T"}; do 
  #     files=("${files[@]}" "$i$j")
  #     echo ${files[@]}
  #   done
  # done

  files=`ls data/*.fq data/*.fa`
  echo ${files[@]}
  # files <- c(list.files("data",pattern=".fq$"),list.files("data",pattern=".fa$"))
  # files<-files[order(files)]
  # write.table(cbind(files,ad3,ad3i,ad5,pr3,pre),"files.txt",quote=F)

  i="data/46N_v1.fq"
  for i in ${files[@]}; do
    mem=`expr 24000 + $core \* 4000`
    if [ "${i}" == "data/33preU_v2.fq" ]; then 
      mem="80000"
    fi
    bsub -M $mem -n $core -R "rusage[mem=$mem,numcpus=$core.00] span[ptile=$core]" "$HOME/bin/mapping_and_RNA_catalog.sh $i"
    sleep 2
  done # i
  sleep 4
  bjobs | grep -v "JOBID" | grep "$core\*" | wc -l
  bjobs | grep -v "JOBID" | wc -l
  grep -o "Total Sequences</td><td>[0-9]</td>" qc_raw/*.html
  grep -c ^@ $out/*/*_r4.fastq > $out/raw_reads_count.txt &
  grep -c ^@ $out/*/*.fastq |  grep -v _r4.fastq > $out/trimmed_reads_count.txt &
  grep ^hg38_repeatmasker_2 $out/*/*/*_no_feature.txt | sed 's!${out}!!' > $out/no_feature_count.txt
  ls $out/*/*.fasta | xargs -l -I FILE gawk '!/^>/{ i++; seq[i] = length($1); sum += length($1) } END { asort(seq,sorted,"@val_num_asc"); print("FILE\t" sorted[int(i/2)] "\t" sum / i) }' FILE > $out/median_and_mean_reads_length.txt &
# ./analysis_long.Very_sensitive/homo_sapiens/46N_v1/46N_v1.fasta

stats=`ls $out/*/stat.txt | tr "\n" " "`  #  | grep -v 33preU_v2
count=`ls $out/*/stat.txt | wc -l`
columns="1,`seq -s "," 2 2 $((count*2+1)) | tr "\n" " "`"
paste $stats | cut -f $columns  > $out/stats.tsv

# pushd $out
# R --no-save
# stats<-as.matrix(read.delim("stats.tsv"))
# rownames(stats)<-stats[,1]
# stats<-stats[-c(6:7,10,12,14,16,17),-1]
# head(stats)
# stats[,1]
# out<-stats[,-c(1:ncol(stats))]
# for(i in seq(1,ncol(stats),2)){
#   out<-cbind(out,as.numeric(stats[,i])+as.numeric(stats[,i+1]))
#   colnames(out)[ncol(out)]<- colnames(stats)[i]
# }
# head(out)
# write.table(out,"stats_merged.tsv",quote=FALSE)
# q()
# popd

NN=RepeatMasker_rRNA_mergedFeatures.txt
NN=RepeatMasker_tRNA_mergedFeatures.txt
prefix=htseq-count_
TRNA=`find $origout/$DV | grep $prefix | grep $NN`
count=`echo $TRNA | wc -w`
columns="1,`seq -s "," 2 2 $((count*2+1)) | tr "\n" " "`"
echo -e "File\t`echo $TRNA | tr ' ' '\n' | sed -E \"s|.*/$prefix||\" | sed \"s|.$NN||\" | tr '\n' '\t'`" > $out/$NN.tsv
paste `echo $TRNA | tr '\n' ' '` | cut -f $columns >> $out/$NN.tsv
grep tRNA-Gly-GGG $DB/*/RepeatMasker_tRNA_mergedFeatures.gtf > tRNA-Gly-GGG.gtf
grep tRNA-Glu-GAG_ $DB/*/RepeatMasker_tRNA_mergedFeatures.gtf > tRNA-Glu-GAG_.gtf
grep LSU-rRNA_Hsa $DB/*/RepeatMasker_rRNA_mergedFeatures.gtf > LSU-rRNA_Hsa.gtf
grep SSU-rRNA_Hsa $DB/*/RepeatMasker_rRNA_mergedFeatures.gtf > SSU-rRNA_Hsa.gtf

mkdir $out/kraken2.reports
cp $out/*/*.kraken2.report $out/kraken2.reports/
tar czf $origout/$DV/kraken2.reports.tar.gz $origout/$DV/kraken2.reports
# echo ${files[@]} | sed "s@data/@${out}/@g"
tar czf $out/kraken2.output.tar.gz `find $origout/$DV | grep kraken2_report | tr '\n' ' '`
tar czf $out/sRNAflow.output.tar.gz `ls -d $origout/$DV/*/output`
  
conda activate qiime
FASTA=`ls -m $out/*/*.json.biom`
BIOMS=`echo $FASTA | sed 's/ //g'`
echo $BIOMS
merge_otu_tables.py -i $BIOMS -o $out/merged_samples.biom
biom convert --to-json -i $out/merged_samples.biom -o $out/merged_samples.json.biom
conda deactivate


  grep -v ^@ $out/33N_v1/ShortStack.Bowtie2_N1_Seed20_wGaps_trimm/33N_v1.sam | awk '{print $1}' | sort | uniq | wc -l
  grep -v ^@ $out/33N_v1/ShortStack.Bowtie2_N1_Seed20_wGaps_trimm/33N_v1.sam | wc -l


  grep -R "Total reads passing E-value threshold" $out/*/Unmapped_*.log
  grep -R "Total reads passing E-value threshold" $out/*/*.log | grep -v Unmapped_
  grep -R "Total reads passing E-value threshold" $out/*/*rRNA2.log

rm $out/featuresOverlap.txt
for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
do
  grep -R __ambiguous $out/*/*/htseq-nopriority_*.$type.txt | grep -Pv "\t0" | sed "s|$out/||g" | sed "s|/ShortStack.Bowtie2_N1_Seed20_wGaps_trimm/htseq-nopriority||g" | sed "s|_mergedFeatures.txt:__ambiguous||g" >> $out/featuresOverlap.txt
  echo >> $out/featuresOverlap.txt
done

  # bjobs | grep -v "$core\*"  | grep -v "JOBID"  | awk '{print $1}' | xargs -L1 bpeek
  # bjobs | grep "$core\*" | awk '{print $1}' | xargs -L1 bkill
  # bjobs | grep "PEND" | awk '{print $1}' | xargs -L1 bkill
  # bjobs | awk '{print $1}' | grep -v "JOBID" | xargs -L1 bkill
}

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
     fastqc -o $out/qc_raw $ff
     echo -e "Raw\t`grep -c ^@ $ff`" >> $txtLog
     cutadapt --quiet --quality-cutoff=20,20 -a $ad3 -m 1 -o $out/$f/${f}_r3.fastq $ff
     echo -e "QC_and_adapter3\t`grep -c ^@ $out/$f/${f}_r3.fastq`" >> $txtLog
     # cutadapt --quiet --quality-cutoff=10,10 -O 9  -e 0.1 -a $ad3 -o $out/$f/${f}_r1.fastq $ff
     # cutadapt --quiet -O 9  -e 0.1 -g $ad3                -o $out/$f/${f}_r2.fastq $out/$f/${f}_r1.fastq
     # cutadapt --quiet -O 13 -e 0   -g XATCACCGACTGCCCATAG -o $out/$f/${f}_r3.fastq $out/$f/${f}_r2.fastq
     # cutadapt --quiet -O 8 -e 0   -g ^ATCACCGACTGCC      -o $out/$f/${f}_r4.fastq $out/$f/${f}_r3.fastq
     cutadapt --quiet -m 43                         -o $out/$f/${f}_long.fastq      $out/$f/${f}_r3.fastq
     cutadapt --quiet -M 42                         -o $out/$f/${f}_r4.fastq      $out/$f/${f}_r3.fastq
     echo -e "Filter_more42\t`grep -c ^@ $out/$f/${f}_r4.fastq`" >> $txtLog
     cutadapt --quiet -m 10                         -o $out/$f/$f.fastq      $out/$f/${f}_r4.fastq
     echo -e "Filter_less10\t`grep -c ^@ $out/$f/${f}.fastq`" >> $txtLog
     ff=$out/$f/$f.fastq
     fastqc -o $out/qc $ff
# #     filename=${ff##*/}
# #     fasta=$out/$f/${filename%.*}.fasta
    fasta=$out/$f/$f.fasta
    fastq_to_fasta -i $ff -o $fasta
    # send2kraken $fasta $out/$f 0 $core
    # # send2sourmash $fasta $out/$f 0 $core
    # # send2blast $fasta $out/$f 1 $core
    # # send2blast $fasta $out/$f 2 $core
    gawk '!/^>/{ i++; seq[i] = length($1); sum += length($1) } END { asort(seq,sorted,"@val_num_asc"); print("reads_after_trimm\t" i "\nmedian_reads_length\t" sorted[int(i/2)] "\nmean_reads_length\t" sum / i) }' $out/$f/$f.fasta >> $txtLog
  fi 

  # if [ ! -e "$out/$f/Unmapped_$f.fa" ]; then 
    shdir="$out/$f/ShortStack.Bowtie2_N1_Seed20_wGaps_trimm"
    shfile="$out/$f/$f.bam"
    # rm -rf "$shdir"
    bowtie2opt="--time --end-to-end -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal"
    bowtie2 $bowtie2opt -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam >> $out/$f/bowtie2.log 2>&1
    # # bowtie2opt="--time --local -k 21 -p $core -x $DB/$DV $inFasta --un $out/$f/Unmapped_$f.fq --no-unal --score-min G,1,10"
    # # bowtie2 $bowtie2opt -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -U $ff -S $shdir.sam >> $txtLog  2>&1
    $samtools view -uhS -F4 $shdir.sam | $samtools sort -@ $core - -o $shfile
    $shortstack --readfile $shfile --genomefile $DB/$DV.fa --outdir $shdir --bowtie_cores $core --mismatches 1 --bowtie_m 21 --ranmax 20 --keep_quals --inbam --nohp
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
    countOver $shdir/$f.sam $shdir $DB/$DV.gtf $DB $DV >> $txtLog
#   fi
# # fi
ff="$out/$f/Unmapped_$f.fa"
rm $ff
if [ "$inFasta" == "" ]; then
  fastq_to_fasta -i $out/$f/Unmapped_$f.fq -o $ff
else
  ln -rfs $out/$f/Unmapped_$f.fq $ff
fi
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
mkdir -p "$out/$f/kvdb"
sortmerna --ref $rRNAdb --reads $ff --aligned $out/$f/$f\_with_rRNA2 --fastx --other $out/$f/$f\_without_rRNA2 --log -a $core -d "$out/$f/kvdb"  >> $out/$f/sortmerna.log 2>&1 &  # -m 4096 
# mv data/${i}.fq.fastq data/${i}.fq; fastqc -o qc data/${i}.fq
