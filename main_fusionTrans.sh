mkdir fusionTrans && cd fusionTrans
git clone https://github.com/drtamermansour/2018-tamer-transcriptome.git

module unload Python
module swap GNU/4.4.5 GNU/6.2
module load CMake/3.1.0
module use /opt/software/ged-software/modulefiles/
module load anaconda
#conda create -n tamer_sgc python=3.5.4
source activate tamer_sgc

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.pc_transcripts.fa.gz
time bcalm -in gencode.v27.pc_transcripts.fa.gz -kmer-size 31 -abundance-min 1 -out bcalm.gencode31
time python 2018-tamer-transcriptome/index_2.py -k 31 bcalm.gencode31.unitigs.fa gencode.v27.pc_transcripts.fa.gz -o human_cds31 &> human_cds31.index.log
time bcalm -in gencode.v27.pc_transcripts.fa.gz -kmer-size 25 -abundance-min 1 -out bcalm.gencode25
time python 2018-tamer-transcriptome/index_2.py -k 25 bcalm.gencode25.unitigs.fa gencode.v27.pc_transcripts.fa.gz -o human_cds25 &> human_cds25.index.log

## stream genomic resources
## A
#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/005/ERR1050075/ERR1050075_1.fastq.gz | gunzip -c | python 2018-tamer-transcriptome/process_2.py -k 31 human_cds /dev/stdin
## B
#module load SRAToolkit/2.8.2
## -X select no of spots (spots mean fragments) 
## -Z Output to stdout
## --split-spot Split fragments into individual reads
## -I appened ".1" and ".2" read suffices for paired-end dat
## --fasta FASTA only, no qualities
#fastq-dump -X 5000 -Z --split-spot -I --fasta ERR1050078 | sed "s/\(^>.*\).1 /\1\/1 /; s/\(^>.*\).2 /\1\/2 /;" > pair_sample_n10000.fasta 
## https://standage.github.io/streaming-data-from-the-sra-with-fastq-dump.html

cat read_sample_n10000.fastq | python 2018-tamer-transcriptome/process_2.py -k 31 human_cds /dev/stdin > human_cds.process.log
grep "^fusion " human_cds.process.log | awk ‘BEGIN{FS="[{}]";OFS="\n"}{print $2,$4}’ | sort | uniq -c | sort -rnk1
grep "^ambiguous fusion " human_cds.process.log

## run the process_bimode.py
cat read_sample_n10000.fastq | python 2018-tamer-transcriptome/process_bimode.py -k 31 --force_single human_cds31 /dev/stdin > human_cds31.process_2PE_single.log
cat pair_sample_n10000.fasta | python 2018-tamer-transcriptome/process_bimode.py -k 31 --paired human_cds31 /dev/stdin > human_cds31.process_2PE_paired.log
cat pair_sample_n10000.fasta | python 2018-tamer-transcriptome/process_bimode.py -k 31 --paired -u read_sample_n10000.fastq human_cds31 /dev/stdin > human_cds31.process_2PE_both.log

## analysis of results
gunzip -c gencode.v27.pc_transcripts.fa.gz > unzippped_trans.fa
cat unzippped_trans.fa | awk ‘{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}’ > unzippped_trans_unwrapped.fa
dir="/mnt/ls15/scratch/users/mansourt/Tamer/fusionTrans"
id=12175
tail -n+2 human_cds_fusion.info | awk 'BEGIN{FS="\t"}{if($4!=0 && $4!=-1)print $0;}' | grep $id
gene=$(grep -w $id $dir/human_cds.ids | awk '{print $1}');echo $gene;
grep -w $gene $dir/unzippped_trans_unwrapped.fa
grep -A1 -w $gene $dir/unzippped_trans_unwrapped.fa > $gene.fa; cat $gene.fa

## testing the FusionSeq dataset

wget http://rnaseq.gersteinlab.org/fusionseq/datasets/NCIH660.fastq.tar.gz 
tar -xzf NCIH660.fastq.tar.gz 
##
time paste NCIH660_1.fastq NCIH660_2.fastq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | python 2018-tamer-transcriptome/process_bimode.py -k 31 --paired human_cds31 /dev/stdin > human_cds31.process.NCIH660_2PE_paired.log
tail -n8 human_cds31.process.NCIH660_2PE_paired.log

## replace the Ensembl ids with a pharse equal to (Ensembl ids|Gene name)
python 2018-tamer-transcriptome/replace_ids.py -o human_cds31_fusion.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds31_fusion.calc
python 2018-tamer-transcriptome/replace_ids.py -o human_cds31_fusionPairs.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds31_fusionPairs.calc

## for individual reads: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds31_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc31.stat
tail -n+2 human_cds31_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($4=="clear_fusion"){A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc31.stat
cat fusion_clear_calc31.stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusion_clear.sig_calc31.stat.ids

## for read pairs: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds31_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_calc31.stat
tail -n+2 human_cds31_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($3=="clear_fusion"){A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_clear_calc31.stat
cat fusPairs_clear_calc31.stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusPairs_clear.sig_calc31.stat.ids
comm -12 fusion_clear.sig_calc31.stat.ids fusPairs_clear.sig_calc31.stat.ids > allfus_clear.sig_calc31.stat.ids
wc -l fus*31* allfus*31*

## for individual reads: count the genes invloved in all possible fusions, only clear fusions, and only clear fusions happening >1 
tail -n+2 human_cds31_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc31.gene
tail -n+2 human_cds31_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{if($4=="clear_fusion")print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc31.gene
cat fusion_clear_calc31.stat | awk 'BEGIN{FS="[ \t]";OFS="\n"}{if($1>1)print $2,$3}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear.sig_calc31.stat.gene
## add rank and count of aligned reads (both actual counts and TPM)
bash assessExp.sh 31

cat fusion_calc31.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_calc31.gene > temp_sorted
paste fusion_calc31.gene temp_sorted > fusion_calc31.gene.counts
cat fusion_calc31.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_calc31.gene > temp_sorted
paste fusion_calc31.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_calc31.gene.counts_tpm

cat fusion_clear_calc31.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear_calc31.gene > temp_sorted
paste fusion_clear_calc31.gene temp_sorted > fusion_clear_calc31.gene.counts
cat fusion_clear_calc31.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear_calc31.gene > temp_sorted
paste fusion_clear_calc31.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_clear_calc31.gene.counts_tpm

cat fusion_clear.sig_calc31.stat.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear.sig_calc31.stat.gene > temp_sorted
paste fusion_clear.sig_calc31.stat.gene temp_sorted > fusion_clear.sig_calc31.stat.gene.counts
cat fusion_clear.sig_calc31.stat.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal31quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear.sig_calc31.stat.gene > temp_sorted
paste fusion_clear.sig_calc31.stat.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_clear.sig_calc31.stat.gene.counts_tpm

##
time paste NCIH660_1.fastq NCIH660_2.fastq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' | python 2018-tamer-transcriptome/process_bimode.py -k 25 --paired human_cds25 /dev/stdin > human_cds25.process.NCIH660_2PE_paired.log
tail -n8 human_cds25.process.NCIH660_2PE_paired.log

## replace the Ensembl ids with a pharse equal to (Ensembl ids|Gene name)
python 2018-tamer-transcriptome/replace_ids.py -o human_cds25_fusion.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds25_fusion.calc
python 2018-tamer-transcriptome/replace_ids.py -o human_cds25_fusionPairs.v2.calc gencode.v27.pc_transcripts.fa.gz human_cds25_fusionPairs.calc

## for individual reads: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds25_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc25.stat
tail -n+2 human_cds25_fusion.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($4=="clear_fusion"){A[1]=$5;A[2]=$6;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc25.stat
cat fusion_clear_calc25.stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusion_clear.sig_calc25.stat.ids

## for read pairs: count the gene fusion pairs (all possible fusions, only clear fusions, and only clear fusions happening >1) 
tail -n+2 human_cds25_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_calc25.stat
tail -n+2 human_cds25_fusionPairs.v2.calc | awk 'BEGIN{FS=OFS="\t"}{if($3=="clear_fusion"){A[1]=$4;A[2]=$5;asort(A);print A[1],A[2];}}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusPairs_clear_calc25.stat
cat fusPairs_clear_calc25.stat | awk 'BEGIN{FS="[ \t]";OFS="\t"}{if($1>1)print $2,$3}' | sort > fusPairs_clear.sig_calc25.stat.ids
comm -12 fusion_clear.sig_calc25.stat.ids fusPairs_clear.sig_calc25.stat.ids > allfus_clear.sig_calc25.stat.ids
wc -l fus*25* allfus*25*

## for individual reads: count the genes invloved in all possible fusions, only clear fusions, and only clear fusions happening >1 
tail -n+2 human_cds25_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_calc25.gene
tail -n+2 human_cds25_fusion.v2.calc | awk 'BEGIN{FS="\t";OFS="\n"}{if($4=="clear_fusion")print $5,$6;}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear_calc25.gene
cat fusion_clear_calc25.stat | awk 'BEGIN{FS="[ \t]";OFS="\n"}{if($1>1)print $2,$3}' | sort | uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > fusion_clear.sig_calc25.stat.gene
## add rank and count of aligned reads (both actual counts and TPM)
bash assessExp.sh 25

cat fusion_calc25.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_calc25.gene > temp_sorted
paste fusion_calc25.gene temp_sorted > fusion_calc25.gene.counts
cat fusion_calc25.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_calc25.gene > temp_sorted
paste fusion_calc25.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_calc25.gene.counts_tpm

cat fusion_clear_calc25.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear_calc25.gene > temp_sorted
paste fusion_clear_calc25.gene temp_sorted > fusion_clear_calc25.gene.counts
cat fusion_clear_calc25.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear_calc25.gene > temp_sorted
paste fusion_clear_calc25.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_clear_calc25.gene.counts_tpm

cat fusion_clear.sig_calc25.stat.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_count.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear.sig_calc25.stat.gene > temp_sorted
paste fusion_clear.sig_calc25.stat.gene temp_sorted > fusion_clear.sig_calc25.stat.gene.counts
cat fusion_clear.sig_calc25.stat.gene |  awk 'BEGIN{FS="[ :]"}{print $3}' | grep -n -Fwf - kal25quant/agg_id_tpm.tsv | awk 'BEGIN{FS="[:\t]";OFS="\t";}{print $2,$1,$3}' > temp
awk 'BEGIN{FS="[ :\t]"} FNR==NR {x2[$1] = $0; next} $3 in x2 {print x2[$3]}' temp fusion_clear.sig_calc25.stat.gene > temp_sorted
paste fusion_clear.sig_calc25.stat.gene.counts temp_sorted | awk 'BEGIN{FS="[ \t]";OFS="\t";print "Freq\tGene\tcount.rank\tcount\ttpm.rank\ttpm";}{print $1,$2,$4,$5,$7,$8}' > fusion_clear.sig_calc25.stat.gene.counts_tpm


##
#diff --new-line-format="" --unchanged-line-format=""  fusion_clear.sig_calc.stat.ids NCIH660/fusion_clear.sig_calc.stat.ids
comm -12 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > com.ids
comm -23 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > uniq25.ids
comm -13 allfus_clear.sig_calc31.stat.ids allfus_clear.sig_calc25.stat.ids > uniq31.ids


