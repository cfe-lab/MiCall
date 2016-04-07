import os
import argparse
import itertools
from csv import DictReader
from matplotlib import pyplot as plt
from operator import itemgetter

from micall.core import project_config, aln2counts


def coverage_plot(amino_csv):
    """ Generate coverage plots.

    @param amino_csv: an open file object that holds amino acid frequencies
    @return: a list of full paths to the image files.
    """
    # imports project information from JSON
    path, filename = os.path.split(amino_csv.name)
    prefix = filename.split('.')[0]
    projects = project_config.ProjectConfig.loadDefault()
    reader = DictReader(amino_csv)
    paths = []

    for region, group in itertools.groupby(reader, itemgetter('region')):
        # use region to retrieve coordinate reference
        coord_ref = projects.getReference(region)
        x = [i+1 for i in range(len(coord_ref))]
        y = [0]*len(x)  # coverage
        for row in group:
            pos = int(row['refseq.aa.pos'])-1  # adjust from 1-index
            total = sum([int(row[aa]) for aa in aln2counts.AMINO_ALPHABET])
            y[pos] = total

        plt.step(x, y, linewidth=3, where='mid')
        plt.xlim([1, len(coord_ref)])
        plt.ylim([0.1, 1000000])
        plt.yscale('log')
        plt.xlabel('Reference coordinates (AA)', fontsize=18)
        plt.ylabel('Coverage', fontsize=18)
        figname = '%s.%s.png' % (prefix, region)
        dest = os.path.join(path, figname)
        paths.append(dest)
        plt.savefig(dest)  # write image to file
        plt.clf()

    return paths  # locations of image files


def parse_args():
    parser = argparse.ArgumentParser(description='Generate coverage plots from MiCall outputs.')
    parser.add_argument('amino_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing amino acid frequency outputs.')
    return parser.parse_args()


def main():
    args = parse_args()
    coverage_plot(amino_csv=args.amino_csv)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/sam_g2p.py -h
    main()
