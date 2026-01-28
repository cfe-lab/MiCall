import argparse
import logging
from typing import Optional
from datetime import datetime
import shutil
from pathlib import Path
from subprocess import CalledProcessError
import subprocess
from tempfile import mkdtemp

from micall.utils.cache import cached
from micall.utils.work_dir import WorkDir


HAPLOFLOW = "haploflow"
logger = logging.getLogger(__name__)


def count_fasta_sequences(file_path: Path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))


@cached("denovo[haploflow]")
def run_subprocess(
    fastq1: Path,
    fastq2: Path,
    merged_contigs_csv: Optional[Path],
) -> Path:

    if merged_contigs_csv is not None:
        # TODO: implement this.
        logger.error("Haploflow implementation does not support contig extensions yet.")

    work_dir = WorkDir.get()
    old_tmp_dirs = work_dir.glob('assembly_*')
    for old_tmp_dir in old_tmp_dirs:
        shutil.rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = Path(mkdtemp(dir=work_dir, prefix='assembly_'))

    joined_path = tmp_dir / 'joined.fastq'
    subprocess.run(['merge-mates',
                    str(fastq1),
                    str(fastq2),
                    '--interleave',
                    '-o', str(joined_path)],
                   check=True)

    haplo_args = {'long': 0,
                  'filter': 80,
                  'thres': -1,
                  'strict': 5,
                  'error': 0.02,
                  'kmer': 23,
                  'merge': False,
                  'scaffold': False,
                  'patch': False,
                  'ref': None,
                  'RP': False,
                  }

    assembly_out_path = tmp_dir / 'haplo_out'
    contigs_fasta_path = assembly_out_path / 'contigs.fa'

    assembly_out_path.mkdir(exist_ok=True, parents=True)
    contigs_fasta_path.touch()

    haplo_cmd = [HAPLOFLOW,
                 '--read-file', str(joined_path),
                 '--out', str(assembly_out_path),
                 '--k', str(haplo_args['kmer']),
                 '--error-rate', str(haplo_args['error']),
                 '--strict', str(haplo_args['strict']),
                 '--filter', str(haplo_args['filter']),
                 '--thres', str(haplo_args['thres']),
                 '--long', str(haplo_args['long'])]
    try:
        subprocess.run(haplo_cmd, check=True)
    except CalledProcessError:
        logger.warning('Haploflow failed to assemble.', exc_info=True)

    return contigs_fasta_path


def denovo(fastq1_path: Path,
           fastq2_path: Path,
           fasta: Path,
           merged_contigs_csv: Optional[Path] = None,
           ):
    """ Use de novo assembly to build contigs from reads.

    :param fastq1: FASTQ file for read 1 reads
    :param fastq2: FASTQ file for read 2 reads
    :param fasta: file to write assembled contigs to
    :param work_dir: path for writing temporary files
    :param merged_contigs_csv: file to read contigs that were merged from
        amplicon reads
    """

    start_time = datetime.now()
    contigs_fasta_path = run_subprocess(fastq1_path,
                                        fastq2_path,
                                        merged_contigs_csv)
    shutil.copy(contigs_fasta_path, fasta)

    duration = datetime.now() - start_time
    contig_count = count_fasta_sequences(contigs_fasta_path)
    logger.info('Assembled %d contigs in %s (%ds) on %s.',
                contig_count,
                duration,
                duration.total_seconds(),
                fastq1_path)


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
