from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import shutil
import sys
import tempfile

from micall.core.trim_fastqs import TrimSteps, trim


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('micall')


def concatenate_files(input_file1, input_file2, output_file):
    with open(input_file1, 'rb') as src1, \
         open(input_file2, 'rb') as src2, \
         open(output_file, 'wb') as dst:

        shutil.copyfileobj(src1, dst)
        shutil.copyfileobj(src2, dst)


def merge_fastqs(args):
    with tempfile.NamedTemporaryFile() as trimmed_fastq1_a, \
         tempfile.NamedTemporaryFile() as trimmed_fastq2_a, \
         tempfile.NamedTemporaryFile() as trimmed_fastq1_b, \
         tempfile.NamedTemporaryFile() as trimmed_fastq2_b:

        logger.info('Processing reads of Sample A.')

        trim((args.fastq1_a, args.fastq2_a),
             args.bad_cycles_a_csv,
             (trimmed_fastq1_a.name, trimmed_fastq2_a.name),
             use_gzip=not args.unzipped)

        logger.info('Processing reads of Sample B.')

        trim((args.fastq1_b, args.fastq2_b),
             args.bad_cycles_b_csv,
             (trimmed_fastq1_b.name, trimmed_fastq2_b.name),
             use_gzip=not args.unzipped)

        logger.info('Merging resuling reads files.')

        concatenate_files(trimmed_fastq1_a.name, trimmed_fastq1_b.name,
                          args.fastq1_result)
        concatenate_files(trimmed_fastq2_a.name, trimmed_fastq2_b.name,
                          args.fastq2_result)

    logger.info('Done.')


def main(argv) -> int:
    parser = ArgumentParser(
        description="Combine and filter the FASTQ files from two samples into a single output file.",
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "fastq1_a",
        help="FASTQ file containing forward reads of sample A",
    )
    parser.add_argument(
        "fastq2_a",
        help="FASTQ file containing reverse reads of sample A",
    )
    parser.add_argument(
        "fastq1_b",
        help="FASTQ file containing forward reads of sample B",
    )
    parser.add_argument(
        "fastq2_b",
        help="FASTQ file containing reverse reads of sample B",
    )
    parser.add_argument(
        "fastq1_result",
        help="Resulting combined FASTQ file containing forward reads",
    )
    parser.add_argument(
        "fastq2_result",
        help="Resulting combined FASTQ file containing reverse reads",
    )
    parser.add_argument(
        "--bad_cycles_a_csv",
        help="list of tiles and cycles rejected for poor quality in sample A",
    )
    parser.add_argument(
        "--bad_cycles_b_csv",
        help="list of tiles and cycles rejected for poor quality in sample B",
    )
    parser.add_argument(
        '--unzipped', '-u',
        action='store_true',
        help='Set if the original FASTQ files are not compressed',
    )

    args = parser.parse_args(argv)
    merge_fastqs(args)
    return 0


if __name__ == '__main__': exit(main(sys.argv[1:]))
