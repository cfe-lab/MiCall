cd /mnt/output
md5sum $1
md5sum $2

merge-mates $1 $2 --interleave --gzip -o pro.joined.fq

iva --fr pro.joined.fq -t 8 pro
mv pro/contigs.fasta pro.fasta

bwa bwasw /opt/micall/micall/utils/hcv_geno/hxb2.fasta pro.fasta | samtools view -bS - | samtools sort - > prov.bam
samtools index prov.bam

python /opt/micall/micall/utils/plot_proviral.py prov.bam
#mv /tmp/summary.csv $5
mv /tmp/alignment.svg $3
