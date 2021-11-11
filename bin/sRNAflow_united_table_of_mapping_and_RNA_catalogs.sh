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
# for type in {"miRBase_mature_mergedFeatures","miRBase_hairpin_mergedFeatures","miRNA_mergedFeatures","piRNAdb_mergedFeatures","piRNAbank_mergedFeatures","piRBase_mergedFeatures","YRNA(misc_RNA)_mergedFeatures","snRNA_mergedFeatures","snoRNA_mergedFeatures","tRF_mergedFeatures","tRNAhalves_mergedFeatures","GtRNAdb_mergedFeatures","rRNA_mergedFeatures","MT_mergedFeatures","Mt_rRNA_mergedFeatures","Mt_tRNA_mergedFeatures","other(MT)_mergedFeatures","vault_RNA_mergedFeatures","lncRNA_mergedFeatures","lncipedia_hc_mergedFeatures","lncipedia_mergedFeatures","noYorPiwi(misc_RNA)_mergedFeatures","protein_coding_mergedFeatures","processed_pseudogene_mergedFeatures","Other_types_mergedFeatures","RepeatMasker_tRNA_mergedFeatures","RepeatMasker_rRNA_mergedFeatures","RepeatMasker_mergedFeatures","Ensembl_genes_mergedFeatures"}
# do
#   grep -R __ambiguous $out/*/*/htseq-nopriority_*.$type.txt | grep -Pv "\t0" | sed "s|$out/||g" | sed "s|/ShortStack/htseq-nopriority||g" | sed "s|_mergedFeatures.txt:__ambiguous||g" >> $out/featuresOverlap.txt
#   echo >> $out/featuresOverlap.txt
# done

