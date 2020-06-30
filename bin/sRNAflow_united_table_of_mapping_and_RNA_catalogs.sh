#!/bin/bash

out="$1"
pushd $out
stats=`ls $out/*/stat.txt | tr "\n" " "`
count=`ls $out/*/stat.txt | wc -l`
columns="1,`seq -s "," 2 2 $((count*2+1)) | tr "\n" " "`"
paste $stats | cut -f $columns  > $out/stats.tsv
popd
# conda activate qiime
# FASTA=`ls -m $out/*/*.json.biom`
# BIOMS=`echo $FASTA | sed 's/ //g'`
# echo $BIOMS
# merge_otu_tables.py -i $BIOMS -o $out/merged_samples.biom
# biom convert --to-json -i $out/merged_samples.biom -o $out/merged_samples.json.biom
# conda deactivate

# rm $out/featuresOverlap.txt
# for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
# do
#   grep -R __ambiguous $out/*/*/htseq-nopriority_*.$type.txt | grep -Pv "\t0" | sed "s|$out/||g" | sed "s|/ShortStack/htseq-nopriority||g" | sed "s|_mergedFeatures.txt:__ambiguous||g" >> $out/featuresOverlap.txt
#   echo >> $out/featuresOverlap.txt
# done

function send2cluster {
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

  grep -R "Total reads passing E-value threshold" $out/*/Unmapped_*.log
  grep -R "Total reads passing E-value threshold" $out/*/*.log | grep -v Unmapped_
  grep -R "Total reads passing E-value threshold" $out/*/*rRNA2.log

rm $out/featuresOverlap.txt
for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vaultRNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
do
  grep -R __ambiguous $out/*/*/htseq-nopriority_*.$type.txt | grep -Pv "\t0" | sed "s|$out/||g" | sed "s|/ShortStack/htseq-nopriority||g" | sed "s|_mergedFeatures.txt:__ambiguous||g" >> $out/featuresOverlap.txt
  echo >> $out/featuresOverlap.txt
done
}
