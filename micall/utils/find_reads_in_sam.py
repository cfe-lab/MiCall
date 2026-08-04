from argparse import ArgumentParser, FileType
from io import StringIO
from operator import attrgetter
from subprocess import run
import sys

from pysam import AlignmentFile


def parse_args():
    parser = ArgumentParser(description="Find reads that map to a position.")
    parser.add_argument("bam_file", type=FileType('rb'))
    parser.add_argument('target_region', help='Positions to look for (e.g., HCV-1a:100-200).')
    return parser.parse_args()


def main():
    args = parse_args()
    with args.bam_file:
        bam_reader = AlignmentFile(args.bam_file)
        
        if not bam_reader.has_index():
            print('Adding index...')
            index_args = ['samtools', 'index', args.bam_file.name]
            run(index_args, check=True)
            args.bam_file.seek(0)  # Go back to start of header.
            bam_reader = AlignmentFile(args.bam_file)
            bam_reader.check_index()
        x = bam_reader.parse_region(region=args.target_region)
        sequences = sorted(bam_reader.fetch(region=args.target_region),
                           key=attrgetter('qname'))
        print(len(sequences))
        for seq in sequences:
            print(seq.qname)
 
if __name__ == '__main__':
    main()
