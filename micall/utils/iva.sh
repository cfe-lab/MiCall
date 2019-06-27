set -e

cd /mnt/output
md5sum $1
md5sum $2
md5sum $3
md5sum $4

merge-mates $1 $2 --interleave --gzip -o wg.joined.fq
merge-mates $3 $4 --interleave --gzip -o mid.joined.fq

iva --fr wg.joined.fq -t 8 wg
mv wg/contigs.fasta $6

iva --fr mid.joined.fq -t 8 mid
mv mid/contigs.fasta $7

bwa bwasw /opt/micall/micall/utils/hcv_geno/hcv.fasta $6 | samtools view -bS - | samtools sort - > wg.bam
samtools index wg.bam

bwa bwasw /opt/micall/micall/utils/hcv_geno/hcv.fasta $7 | samtools view -bS - | samtools sort - > mid.bam
samtools index mid.bam

python /opt/micall/micall/utils/plot.py wg.bam mid.bam `cat $5`
mv /tmp/alignment.svg $8
mv /tmp/subtyping.csv $9
