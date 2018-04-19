mkdir fusionTrans && cd fusionTrans
git clone https://github.com/drtamermansour/2018-tamer-transcriptome.git
work_dir=$(pwd)
scripts=$(pwd)/2018-tamer-transcriptome

## start the VM on HPC
module unload Python
module swap GNU/4.4.5 GNU/6.2
module load CMake/3.1.0
module use /opt/software/ged-software/modulefiles/
module load anaconda
#conda create -n tamer_sgc python=3.5.4
source activate tamer_sgc

## Create the index
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.pc_transcripts.fa.gz
time bcalm -in gencode.v27.pc_transcripts.fa.gz -kmer-size 31 -abundance-min 1 -out bcalm.gencode31
time python $scripts/index_2.py -k 31 bcalm.gencode31.unitigs.fa gencode.v27.pc_transcripts.fa.gz -o human_cds31 &> human_cds31.index.log
time bcalm -in gencode.v27.pc_transcripts.fa.gz -kmer-size 25 -abundance-min 1 -out bcalm.gencode25
time python $scripts/index_2.py -k 25 bcalm.gencode25.unitigs.fa gencode.v27.pc_transcripts.fa.gz -o human_cds25 &> human_cds25.index.log

## stream genomic resources
## data from the 1000 genome: https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/hgsv_sv_discovery/illumina_rna.sequence.index
## A
#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/005/ERR1050075/ERR1050075_1.fastq.gz | gunzip -c | python $scripts/process_2.py -k 31 human_cds /dev/stdin
## B
#module load SRAToolkit/2.8.2
## -X select no of spots (spots mean fragments) 
## -Z Output to stdout
## --split-spot Split fragments into individual reads
## -I appened ".1" and ".2" read suffices for paired-end dat
## --fasta FASTA only, no qualities
#fastq-dump -X 5000 -Z --split-spot -I --fasta ERR1050078 | sed "s/\(^>.*\).1 /\1\/1 /; s/\(^>.*\).2 /\1\/2 /;" > pair_sample_n10000.fasta 
## https://standage.github.io/streaming-data-from-the-sra-with-fastq-dump.html

# process_2.py works only in single end mode
cat read_sample_n10000.fastq | python $scripts/process_2.py -k 31 human_cds /dev/stdin > human_cds.process.log
grep "^fusion " human_cds.process.log | awk ‘BEGIN{FS="[{}]";OFS="\n"}{print $2,$4}’ | sort | uniq -c | sort -rnk1
grep "^ambiguous fusion " human_cds.process.log

## run the process_bimode.py
cat read_sample_n10000.fastq | python $scripts/process_bimode.py -k 31 --force_single human_cds31 /dev/stdin > human_cds31.process_2PE_single.log
cat pair_sample_n10000.fasta | python $scripts/process_bimode.py -k 31 --paired human_cds31 /dev/stdin > human_cds31.process_2PE_paired.log
cat pair_sample_n10000.fasta | python $scripts/process_bimode.py -k 31 --paired -u read_sample_n10000.fastq human_cds31 /dev/stdin > human_cds31.process_2PE_both.log

## analysis of results
gunzip -c gencode.v27.pc_transcripts.fa.gz > unzippped_trans.fa
cat unzippped_trans.fa | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > unzippped_trans_unwrapped.fa
dir="/mnt/ls15/scratch/users/mansourt/Tamer/fusionTrans"
id=12175
tail -n+2 human_cds_fusion.info | awk 'BEGIN{FS="\t"}{if($4!=0 && $4!=-1)print $0;}' | grep $id
gene=$(grep -w $id $dir/human_cds.ids | awk '{print $1}');echo $gene;
grep -w $gene $dir/unzippped_trans_unwrapped.fa
grep -A1 -w $gene $dir/unzippped_trans_unwrapped.fa > $gene.fa; cat $gene.fa

###########################################################################
## testing the FusionSeq dataset
## Sample: NCIH660, Know fusion: TMPRSS2-ERG
wget http://rnaseq.gersteinlab.org/fusionseq/datasets/NCIH660.fastq.tar.gz 
tar -xzf NCIH660.fastq.tar.gz 
##
k=25 #k=31
time paste NCIH660_1.fastq NCIH660_2.fastq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | python $scripts/process_bimode.py -k "$k" --paired human_cds"$k" /dev/stdin > human_cds"$k".process.NCIH660_2PE_paired.log
tail -n8 human_cds"$k".process.NCIH660_2PE_paired.log

## replace the Ensembl ids with a pharse equal to (Ensembl ids|Gene name)
python $scripts/replace_ids.py -o human_cds"$k"_fusion.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds"$k"_fusion.calc
python $scripts/replace_ids.py -o human_cds"$k"_fusionPairs.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds"$k"_fusionPairs.calc

## for individual reads: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds"$k"_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc"$k".stat
tail -n+2 human_cds"$k"_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($4=="clear_fusion"){A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc"$k".stat
#cat fusion_clear_calc"$k".stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusion_clear.sig_calc"$k".stat.ids

## for read pairs: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds"$k"_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_calc"$k".stat
tail -n+2 human_cds"$k"_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($3=="clear_fusion"){A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_clear_calc"$k".stat
#cat fusPairs_clear_calc"$k".stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusPairs_clear.sig_calc"$k".stat.ids
#comm -12 fusion_clear.sig_calc"$k".stat.ids fusPairs_clear.sig_calc"$k".stat.ids > allfus_clear.sig_calc"$k".stat.ids
#wc -l fus*"$k"* allfus*"$k"*
awk -F"[ \t]" 'NR==FNR {h[$2 $3] = $1; next} {if(h[$2 $3]!="") print $1,h[$2 $3],$2,$3}' fusPairs_calc"$k".stat fusion_calc"$k".stat > allfus_calc"$k".stat
awk -F"[ \t]" 'NR==FNR {h[$2 $3] = $1; next} {if(h[$2 $3]!="") print $1,h[$2 $3],$2,$3}' fusPairs_clear_calc"$k".stat fusion_clear_calc"$k".stat > allfus_clear_calc"$k".stat

################
## QC analysis

## for individual reads: count the reads and fusion pairs for each gene invloved in all possible fusions, or only clear fusions
tail -n+2 human_cds"$k"_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc"$k".reads.gene
cat fusion_calc"$k".stat | awk 'BEGIN{FS="[ \t]";OFS="\n"}{print $2,$3}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc"$k".fusions.gene
join -j2 <(sort -k2 fusion_calc"$k".reads.gene) <(sort -k2 fusion_calc"$k".fusions.gene) | awk '{print $2,$3,$1}' | sort -nr > fusion_calc"$k".gene

tail -n+2 human_cds"$k"_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{if($4=="clear_fusion")print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc"$k".reads.gene
cat fusion_clear_calc"$k".stat | awk 'BEGIN{FS="[ \t]";OFS="\n"}{print $2,$3}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc"$k".fusions.gene
join -j2 <(sort -k2 fusion_clear_calc"$k".reads.gene) <(sort -k2 fusion_clear_calc"$k".fusions.gene) | awk '{print $2,$3,$1}' | sort -nr > fusion_clear_calc"$k".gene

## add rank and count of aligned reads (both actual counts and TPM)
bash $scripts/assessExp.sh "$k"

cat fusion_calc"$k".gene |  awk 'BEGIN{FS="[ :]"}{print $4}' | grep -n -Fwf - kal"$k"quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $4 in x2 {print x2[$4]}' temp fusion_calc"$k".gene > temp_sorted
paste fusion_calc"$k".gene temp_sorted > fusion_calc"$k".gene.counts
cat fusion_calc"$k".gene |  awk 'BEGIN{FS="[ :]"}{print $4}' | grep -n -Fwf - kal"$k"quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $4 in x2 {print x2[$4]}' temp fusion_calc"$k".gene > temp_sorted
paste fusion_calc"$k".gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Reads\tFusions\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$3,$5,$6,$8,$9}' > fusion_calc"$k".gene.counts_tpm

cat fusion_clear_calc"$k".gene |  awk 'BEGIN{FS="[ :]"}{print $4}' | grep -n -Fwf - kal"$k"quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $4 in x2 {print x2[$4]}' temp fusion_clear_calc"$k".gene > temp_sorted
paste fusion_clear_calc"$k".gene temp_sorted > fusion_clear_calc"$k".gene.counts
cat fusion_clear_calc"$k".gene |  awk 'BEGIN{FS="[ :]"}{print $4}' | grep -n -Fwf - kal"$k"quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $4 in x2 {print x2[$4]}' temp fusion_clear_calc"$k".gene > temp_sorted
paste fusion_clear_calc"$k".gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Reads\tFusions\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$3,$5,$6,$8,$9}' > fusion_clear_calc"$k".gene.counts_tpm

## select reads with fusion gaps
tail -n+2 human_cds"$k"_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($4=="clear_fusion" && $10!="[0]"){A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear.gap_calc"$k".stat


##
#diff --new-line-format="" --unchanged-line-format=""  fusion_clear.sig_calc.stat.ids NCIH660/fusion_clear.sig_calc.stat.ids
comm -12 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > com.ids
comm -23 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > uniq25.ids
comm -13 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > uniq31.ids

############
cat human_cds25_fusion.fa human_cds25_fusionPairs.fa > human_cds25_allfus.fa
Trinity --seqType fq --max_memory 24G --CPU 6 --output trinity_human_cds25 --run_as_paired --single human_cds25_allfus.fa
cat trinity_human_cds25/Trinity.fasta | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > Trinity.fasta

#module load SSAKE/3.8
#time cat allfus_clear_calc"$k".stat.multi | awk -F"[ :]" '{print $3,$5}' | while read geneA_id geneB_id;do 
#  echo $geneA_id $geneB_id;
#  fus=human_cds"$k"_$geneA_id"_"$geneB_id;
#  awk -F"\t" '{print $5}' human_cds"$k"_fusion.info | grep -nw "$geneA_id" | grep -w "$geneB_id" > $fus.info;
#  cat $fus.info | awk 'BEGIN{FS="[:\t]"}{print $1-1}' | python $scripts/extractFusReads.py -o $fus.fa human_cds25_fusion.fa /dev/stdin;
#  /opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $fus.fa -w 1 -m 16 -o 1;
#  rm $fus.fa.ssake_*.log $fus.fa.ssake_*.singlets $fus.fa.ssake_*.short;
#done
#real 247m5.081s  user 243m23.118s  (~4hours)

## Assemble the reads supporting one-end fusion
> human_cds"$k".top_assemblies.fa;
> human_cds"$k".top_assemblies2.fa;
> human_cds"$k".single_reads.fa;
> human_cds"$k".single_reads2.fa;
time cat allfus_clear_calc"$k".stat | awk -F"[ :|]" '{print $3,$6,$3":"$5"_"$6":"$8}' | while read geneA_id geneB_id id;do
  fus=human_cds"$k"_$geneA_id"_"$geneB_id;
  #awk -F"\t" '{if($4 ~ /_fusion/)print $5}' human_cds"$k"_fusion.info | grep -nw "$geneA_id" | grep -w "$geneB_id" > $fus.info; ## There is non-importnat bug here because this code will grep any line that has the 2 gene names which could be because of part of ambigious segment
  awk -F"\t" '{if($4 ~ /_fusion/)print $5}' human_cds"$k"_fusion.info | grep -nE "$geneA_id.*}.*$geneB_id|$geneB_id.*}.*$geneA_id" > $fus.info;
  cat $fus.info | awk 'BEGIN{FS="[:\t]"}{print $1-1}' | python $scripts/extractFusReads.py -o $fus.fa human_cds25_fusion.fa /dev/stdin;
  reads="$(cat $fus.fa | wc -l)"; echo "No of read lines:" $reads;
  if (( reads > 4 ));then ## The condition assumes fastq format
      echo "multi-reads"
      /opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $fus.fa -w 1 -m 16 -o 1;
      contigs=$(cat $fus.fa.ssake_*.contigs | wc -l);
      if (( contigs > 0 ));then
          echo "contigs found";
          rm $fus.fa.ssake_*.log $fus.fa.ssake_*.singlets $fus.fa.ssake_*.short;
          header=$(head -n1 $fus.fa.ssake_*.contigs)"|"$id
          seq=$(head -n2 $fus.fa.ssake_*.contigs | tail -n1)
          echo $header >>  human_cds"$k".top_assemblies.fa; echo $seq >>  human_cds"$k".top_assemblies.fa;
          polyA=$(echo $seq | grep "AAAAAAAAAAAAAAA");
          if [ ! -n "$polyA" ];then echo $header >>  human_cds"$k".top_assemblies2.fa; echo $seq >>  human_cds"$k".top_assemblies2.fa;fi;
      else
          echo "No contigs";
          rm $fus.fa.ssake_*.log $fus.fa.ssake_*.singlets $fus.fa.ssake_*.short $fus.fa.ssake_*.contigs; 
          awk -F"\t" '{if($4=="clear_fusion") print $5}' human_cds"$k"_fusion.info | grep -nw "$geneA_id" | grep -w "$geneB_id" > $fus.info2;
          cat $fus.info2 | awk 'BEGIN{FS="[:\t]"}{print $1-1}' | python $scripts/extractFusReads.py -o $fus.fa2 human_cds25_fusion.fa /dev/stdin;
          reads="$(cat $fus.fa2 | wc -l)";
          if (( reads > 4 ));then header=">multi_clearReads|"$id;else header=">single_clearRead|"$id;fi ## The condition assumes fastq format
          seq=$(head -n2 $fus.fa2 | tail -n1);
          echo $header >> human_cds"$k".single_reads.fa; echo $seq >> human_cds"$k".single_reads.fa;
          polyA=$(echo $seq | grep "AAAAAAAAAAAAAAA");polyT=$(echo $seq | grep "TTTTTTTTTTTTTTT");
          if [ ! -n "$polyA" ] && [ ! -n "$polyT" ];then echo $header >> human_cds"$k".single_reads2.fa; echo $seq >> human_cds"$k".single_reads2.fa;fi
      fi
  else
      echo "single"; 
      header=">single_read|"$id;
      seq=$(head -n2 $fus.fa | tail -n1);
      echo $header >>  human_cds"$k".single_reads.fa; echo $seq >>  human_cds"$k".single_reads.fa;
  fi  
done


mkdir assembly2 
mv human_cds"$k"_*_*.fa2 assembly2/.
cd assembly2
for f in *.fa2; do /opt/software/SSAKE/3.8--GCC-4.4.5/bin/SSAKE -f $f -w 1 -m 16 -o 1;done
> human_cds"$k".more_assemblies.fa;
for f in *.contigs;do if [ $(cat $f | wc -l) -gt 0 ];then 
  echo $f;
  id=$(echo $f | cut -d"_" -f3-4 | cut -d"." -f1); 
  head -n2 $f | paste - - | awk -v id=$id -v OFS="\n" '{print $1"|"id,$2 }' >> human_cds"$k".more_assemblies.fa;
fi;done
module load BLAST+/2.2.29; k=25;
DB="human_genomic_transcript" ##"human_genomic"
input=human_cds"$k".more_assemblies.fa
blastn -query "$input" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input.blastn
sort -k1,1 -k12,12nr -k11,11n  $input.blastn | sort -u -k1,1 --merge > $input.blastn.best
awk -v FS="\t" -v OFS="\t" '{if($14/$13<0.75)print $14/$13,$0}' $input.blastn.best | sort -n > $input.blastn.best.poor
cd ../

mkdir assembly
mv *ssake*.contigs human_cds"$k"_*_*.info* human_cds"$k"_*_*.fa* assembly/.

module load BLAST+/2.2.29
export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB

DB="human_genomic_transcript" ##"human_genomic"
input1=human_cds"$k".top_assemblies2.fa
blastn -query "$input1" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input1.blastn
input2=human_cds"$k".single_reads2.fa
blastn -query "$input2" -db "$DB" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input2.blastn
#qsub -v input="$input",DB="$DB" ${script_path}/blastn.sh;

sort -k1,1 -k12,12nr -k11,11n  $input1.blastn | sort -u -k1,1 --merge > $input1.blastn.best
awk -v FS="\t" -v OFS="\t" '{if($14/$13<0.75)print $14/$13,$0}' $input1.blastn.best | sort -n > $input1.blastn.best.poor
sort -k1,1 -k12,12nr -k11,11n  $input2.blastn | sort -u -k1,1 --merge > $input2.blastn.best
awk -v FS="\t" -v OFS="\t" '{if($14/$13<0.75)print $14/$13,$0}' $input2.blastn.best | sort -n > $input2.blastn.best.poor

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $input1 --output_fasta_fp $input1.NoBlast --seq_id_fp $input1.blastn.best -n 
blastn -query "$input1.NoBlast" -db "nt" -num_threads 8 -max_target_seqs 1 -outfmt "6 std qlen nident stitle salltitles" > $input1.NoBlast.blastn

## assess single event
geneA_id=11316
geneB_id=1981

cat allfus_calc"$k".stat | grep -w "$geneA_id" | grep -w "$geneB_id"
cat allfus_clear_calc"$k".stat | grep -w "$geneA_id" | grep -w "$geneB_id"
cat fusion_calc"$k".gene | grep -E " ${geneA_id}:| ${geneB_id}:"
cat fusion_clear_calc"$k".gene | grep -E " ${geneA_id}:| ${geneB_id}:"

cat human_cds"$k"_fusion.info | grep -A1 -B1 -E "$geneA_id.*}.*$geneB_id|$geneB_id.*}.*$geneA_id"
cat human_cds"$k"_fusionPairs.info | grep -nE "$geneA_id.*}.*$geneB_id|$geneB_id.*}.*$geneA_id"

cat assembly/human_cds"$k"_"$geneA_id"_"$geneB_id".*.contigs
cat assembly/human_cds"$k"_"$geneA_id"_"$geneB_id".fa


## summary stats
wc -l *.stat
awk '{if($1!=1)print $0;}' allfus_clear_calc25.stat | wc -l ## multi clear read
awk '{if($1==1)print $0;}' allfus_clear_calc25.stat | wc -l ## single clear read
echo $(cat human_cds"$k".top_assemblies.fa | wc -l)/2 | bc  ## assembled
echo $(cat human_cds"$k".top_assemblies2.fa | wc -l)/2 | bc ## Try to align 
echo "$(cat human_cds"$k".top_assemblies.fa | wc -l)/2 - $(cat human_cds"$k".top_assemblies2.fa | wc -l)/2" | bc ## polyA 
echo $(cat $input1.blastn.best | wc -l) ## candidate
echo "$(cat human_cds"$k".top_assemblies2.fa | wc -l)/2 - $(cat $input1.blastn.best | wc -l)" | bc ## failed
echo $(cat human_cds"$k".single_reads.fa | wc -l)/2 | bc ## unassembled
echo $(cat human_cds"$k".single_reads2.fa | wc -l)/2 | bc ## try to align
echo "$(cat human_cds"$k".single_reads.fa | wc -l)/2 - $(cat human_cds"$k".single_reads2.fa | wc -l)/2" | bc ## polyA 
echo $(cat $input2.blastn.best | wc -l) ## candidate
echo "$(cat human_cds"$k".single_reads2.fa | wc -l)/2 - $(cat $input2.blastn.best | wc -l)" | bc ## failed
wc -l *.blastn.best.poor
###########################################################################
## testing the the 1000 genome: https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/hgsv_sv_discovery/illumina_rna.sequence.index

mkdir ERR1050078 && cd ERR1050078
ln -s $work_dir/human_cds"$k".mphf .
ln -s $work_dir/human_cds"$k".arr .
ln -s $work_dir/human_cds"$k".ids .

ln -s $work_dir/gencode.v27.pc_transcripts.fa.gz .

#module load SRAToolkit/2.8.2
#time /opt/software/SRAToolkit/2.8.2--GCC-4.4.5/bin/fastq-dump -Z --split-spot -I --fasta ERR1050078 | sed "s/\(^>.*\)\.1 /\1\/1 /; s/\(^>.*\)\.2 /\1\/2 /;" > ERR1050078.fa ## 69073722 PE read
time /opt/software/SRAToolkit/2.8.2--GCC-4.4.5/bin/fastq-dump --defline-seq '@$sn[_$rn]/$ri' -Z --split-files --fasta ERR1050078 > ERR1050078.fa
# real    32m17.426s  user    29m28.930s  sys     0m39.289s
head -n 20000000 ERR1050078.fa | python $scripts/process_bimode.py -k "$k" --paired human_cds"$k" /dev/stdin > human_cds"$k".process.ERR1050078_2PE_paired.log ## test 5 million PE reads
#time /opt/software/SRAToolkit/2.8.2--GCC-4.4.5/bin/fastq-dump -Z --split-spot -I --fasta ERR1050078 | sed "s/\(^>.*\)\.1 /\1\/1 /; s/\(^>.*\)\.2 /\1\/2 /;" | python $scripts/process_bimode.py -k "$k" --paired ../human_cds"$k" /dev/stdin > human_cds"$k".process.ERR1050078_2PE_paired.log
# killed after: real    125m1.772s user    126m41.219s sys     4m9.716s


## get the fusion reads
# real    106m36.833s  user    101m54.025s  sys     3m32.269s



############################
mkdir fusionMap && cd fusionMap
wget http://omicsoft.com/downloads/FusionMap/FusionMap_2018-01-08.zip
unzip FusionMap_2018-01-08.zip
cd FusionMap_2018-01-08/TestDataset/input
gunzip DatasetP2_SimulatedReads_*.gz

cat DatasetP2_SimulatedReads_1.fastq | sed "s/\:1/\/1/;" > DatasetP2_SimulatedReads_1.fq
cat DatasetP2_SimulatedReads_2.fastq | sed "s/\:2/\/2/;" > DatasetP2_SimulatedReads_2.fq
module swap GNU GNU/4.8.3
module load khmer/2.0
interleave-reads.py DatasetP2_SimulatedReads_1.fq DatasetP2_SimulatedReads_2.fq -o simultatedData_interleaved.fq

time cat simultatedData_interleaved.fq | python $scripts/process_bimode.py -k "$k" --paired human_cds"$k" /dev/stdin > human_cds"$k".process.ERR1050078_2PE_paired.log
# real    0m27.582s  user    0m15.328s  sys     0m1.656s

time bash $scripts/analysis.sh "25" "$scripts";
# real    0m49.938s  user    0m36.199s  sys     0m6.839s

time bash $scripts/analysis2.sh "25" "$scripts";
# real    0m25.433s  user    0m15.676s  sys     0m3.187s

