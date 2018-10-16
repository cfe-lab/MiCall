import argparse
import sys
from subprocess import Popen, PIPE
from tempfile import TemporaryDirectory
from Bio import SeqIO
from micall.utils.gt import genotype

PEAR = "/opt/bin/pear"
SAVAGE = "/opt/savage_wrapper.sh"
IVA = "iva"
IS_SAVAGE_ENABLED = False

ASSEMBLY_TIMEOUT=1800 # 30 minutes

def geno(f, contigs):
    rec = 0
    gts = genotype(f)
    for record in SeqIO.parse(f, "fasta"):
        if rec in gts:
            print("{},{}".format(gts[rec], record.seq), file=contigs)
        else:
            print("unknown,{}".format(record.seq), file=contigs)
        rec += 1

def denovo(fq1, fq2, contigs):
    prefix = "sample"

    with TemporaryDirectory() as tmp_dir:
        iva_proc = Popen([\
                IVA, "-f", fq1, "-r", fq2, "-t", "1", "iva"],
                cwd=tmp_dir)

        if IS_SAVAGE_ENABLED:
            pear_proc = Popen([PEAR, "-f", fq1, "-r", fq2, "-o", prefix],
                    cwd=tmp_dir)

            if pear_proc.wait():
                raise Exception

            savage_proc = Popen([\
                        SAVAGE, "--split", "1", "-s",\
                        "{}/sample.assembled.fastq".format(tmp_dir), "-p1",\
                        "{}/sample.unassembled.forward.fastq".format(tmp_dir),
                        "-p2",
                        "{}/sample.unassembled.reverse.fastq".format(tmp_dir),
                        "-t", "1",\
                        "--merge_contigs", "0.01", "--overlap_len_stage_c",\
                        "100",\
                ], cwd=tmp_dir, shell=True)

            savage_proc.wait(timeout=ASSEMBLY_TIMEOUT)
            try:
                if not savage_proc.wait():
                    geno("{}/contigs_stage_c.fasta".format(tmp_dir), contigs)
                else:
                    print("savage exits with error", file=sys.stderr)
            except:
                print("savage did not finish", file=sys.stderr)

        if not iva_proc.wait(timeout=ASSEMBLY_TIMEOUT):
            geno("{}/iva/contigs.fasta".format(tmp_dir), contigs)
        else:
            print("iva exits with error", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('fastq2')
    parser.add_argument('contigs', type=argparse.FileType('w'))

    args = parser.parse_args()
    denovo(args.fastq1, args.fastq2, args.contigs)
