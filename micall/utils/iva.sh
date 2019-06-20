cd /mnt/output
md5sum $1
md5sum $2
md5sum $3
md5sum $4

merge-mates $1 $2 --interleave --gzip -o wg.joined.fq
merge-mates $3 $4 --interleave --gzip -o mid.joined.fq

iva --fr wg.joined.fq -t 8 wg
mv wg/contigs.fasta $5

iva --fr mid.joined.fq -t 8 mid
mv mid/contigs.fasta $6

bwa bwasw /opt/micall/micall/utils/hcv_geno/hcv.fasta $5 | samtools view -bS - | samtools sort - > wg.bam
samtools index wg.bam

bwa bwasw /opt/micall/micall/utils/hcv_geno/hcv.fasta $6 | samtools view -bS - | samtools sort - > mid.bam
samtools index mid.bam

python /opt/micall/micall/utils/plot.py wg.bam mid.bam
mv /tmp/alignment.svg $7
