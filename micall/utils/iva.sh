cd /mnt/output
merge-mates $1 $2 --interleave -o wg.joined.fq
merge-mates $3 $4 --interleave -o mid.joined.fq

iva --fr wg.joined.fq -t 8 wg
mv wg/contigs.fasta wg.fasta

iva --fr mid.joined.fq -t 8 mid
mv mid/contigs.fasta mid.fasta

