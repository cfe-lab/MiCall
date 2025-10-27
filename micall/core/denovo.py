import argparse
import logging
from pathlib import Path
from typing import Optional
from csv import DictReader
from datetime import datetime
from shutil import rmtree, copyfileobj
from subprocess import PIPE, CalledProcessError, STDOUT
import subprocess
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.utils.cache import cached
from micall.utils.work_dir import WorkDir


IVA = "iva"
logger = logging.getLogger(__name__)


def count_fasta_sequences(file_path: Path) -> int:
    with file_path.open() as file:
        return sum(1 for line in file if line.startswith(">"))


@cached("denovo[iva]")
def run_subprocess(
    fastq1: Path,
    fastq2: Path,
    merged_contigs_csv: Optional[Path],
) -> Path:
    """
    Run IVA de novo assembly subprocess.

    Uses work_dir from WorkDir dynamic scoping for temporary file storage.
    """
    # Get work_dir from dynamic scope - required to be set by caller
    work_dir = WorkDir.get()

    old_tmp_dirs = work_dir.glob("assembly_*")
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = Path(mkdtemp(dir=work_dir, prefix="assembly_"))

    joined_path = tmp_dir / "joined.fastq"
    subprocess.run(
        ["merge-mates", fastq1, fastq2, "--interleave", "-o", str(joined_path)],
        check=True,
        cwd=tmp_dir,
    )
    iva_out_path = tmp_dir / "iva_out"
    contigs_fasta_path = iva_out_path / "contigs.fasta"
    iva_args = [IVA, "--fr", str(joined_path), "-t", "2"]
    if merged_contigs_csv is not None:
        seeds_fasta_path = tmp_dir / "seeds.fasta"
        with open(seeds_fasta_path, "w") as seeds_fasta, \
             merged_contigs_csv.open() as merged_contigs_reader:
            SeqIO.write(
                (
                    SeqRecord(Seq(row["contig"]), f"seed-{i}", "", "")
                    for i, row in enumerate(DictReader(merged_contigs_reader))
                ),
                seeds_fasta,
                "fasta",
            )
            seeds_size = seeds_fasta.tell()
        if seeds_size > 0:
            iva_args.extend(["--contigs", str(seeds_fasta_path), "--make_new_seeds"])
    iva_args.append(str(iva_out_path))
    try:
        subprocess.run(iva_args, check=True, stdout=PIPE, stderr=STDOUT)
    except CalledProcessError as ex:
        output = ex.output and ex.output.decode("UTF8")
        if output != "Failed to make first seed. Cannot continue\n":
            logger.warning("iva failed to assemble.", exc_info=True)
            logger.warning(output)
        contigs_fasta_path.touch()

    return contigs_fasta_path


def denovo(
    fastq1: Path,
    fastq2: Path,
    fasta: Path,
    merged_contigs_csv: Optional[Path] = None,
):
    """Use de novo assembly to build contigs from reads.

    :param fastq1: FASTQ file for read 1 reads
    :param fastq2: FASTQ file for read 2 reads
    :param fasta: file to write assembled contigs to
    :param work_dir: path for writing temporary files (optional, uses dynamic scope if not provided)
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    """

    start_time = datetime.now()

    contigs_fasta_path = run_subprocess(fastq1, fastq2, merged_contigs_csv)
    with contigs_fasta_path.open() as reader, \
         fasta.open("w") as writer:
        copyfileobj(reader, writer)

    duration = datetime.now() - start_time
    contig_count = count_fasta_sequences(contigs_fasta_path)
    logger.info(
        "Assembled %d contigs in %s (%ds) on %s.",
        contig_count,
        duration,
        duration.total_seconds(),
        fastq1,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(
        description="A script to perform de novo assembly of reads to build contigs."
    )
    parser.add_argument(
        "fastq1",
        type=Path,
        help="Path to the FASTQ file containing read 1 of paired-end sequencing data.",
    )
    parser.add_argument(
        "fastq2",
        type=Path,
        help="Path to the FASTQ file containing read 2 of paired-end sequencing data.",
    )
    parser.add_argument(
        "fasta",
        type=Path,
        help="Path to the output FASTA file where assembled contigs will be written.",
    )

    args = parser.parse_args()
    with WorkDir.using(Path.cwd()):
        denovo(args.fastq1, args.fastq2, args.fasta)
