set -e

cd /mnt/output
md5sum $1
md5sum $2

merge-mates $1 $2 --interleave --gzip -o pro.joined.fq

iva --fr pro.joined.fq -t 8 pro
mv pro/contigs.fasta pro.fasta

bwa bwasw /opt/micall/micall/utils/hcv_geno/hxb2.fasta pro.fasta | samtools view -bS - | samtools sort - > prov.bam
samtools index prov.bam

bwa index pro.fasta
bwa mem pro.fasta $1 $2 | samtools view -bS - | samtools sort - > remap.bam
samtools index remap.bam

bwa mem /opt/micall/micall/utils/hcv_geno/hxb2.fasta $1 $2 | samtools view -bS - | samtools sort - > hxb2_remap.bam
samtools index hxb2_remap.bam


python /opt/micall/micall/utils/plot_proviral.py prov.bam remap.bam hxb2_remap.bam `cat $3`
mv /tmp/summary.csv $5
mv /tmp/alignment.svg $4
mv /tmp/alignment.png $7
mv pro.fasta $6

