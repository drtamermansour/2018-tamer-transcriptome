k=$1

module swap GNU GNU/4.9
module load kallisto/0.43.0
kallisto index -k $k -i kal${k} gencode.v27.pc_transcripts.fa.gz 
kallisto quant -i kal${k} -o kal${k}quant -t4 NCIH660_1.fastq NCIH660_2.fastq 

cd kal${k}quant
#cat abundance.tsv | awk -F "\t" '{if($5!=0)print $0}' > eff_abundance.tsv
cat abundance.tsv | awk -F "\t" '{print $1}' | awk 'BEGIN{FS=OFS="|"}{print $2,$6}' > id.tsv

cat abundance.tsv | awk -F "\t" '{print $5}' > tpm.tsv
paste id.tsv tpm.tsv > id_tpm.tsv
tail -n+2 id_tpm.tsv | awk 'BEGIN{FS=OFS="\t"}{a[$1]+=$2}END{for(i in a) print i,a[i]}' | sort -nrk2,2 > agg_id_tpm.tsv

cat abundance.tsv | awk -F "\t" '{print $4}' > count.tsv
paste id.tsv count.tsv > id_count.tsv
tail -n+2 id_count.tsv | awk 'BEGIN{FS=OFS="\t"}{a[$1]+=$2}END{for(i in a) print i,a[i]}' | sort -nrk2,2 > agg_id_count.tsv       

tail -n+2 id_tpm.tsv | awk 'BEGIN{FS=OFS="\t"}{sum+=$2}END{print sum}'
tail -n+2 id_count.tsv | awk 'BEGIN{FS=OFS="\t"}{sum+=$2}END{print sum}'
