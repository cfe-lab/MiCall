from csv import DictReader
from itertools import groupby
from operator import itemgetter

from ..core.aln2counts import AMINO_ALPHABET


def read_aminos(amino_csv, min_fraction, reported_regions=None):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    for region, rows in groupby(DictReader(amino_csv),
                                itemgetter('region')):
        if reported_regions is not None and region not in reported_regions:
            continue
        aminos = []
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = sum(counts)
            min_count = max(1, coverage * min_fraction)  # needs at least 1
            pos_aminos = [report_names[i]
                          for i, count in enumerate(counts)
                          if count >= min_count and report_names[i] != '*']
            ins_count = int(row['ins'])
            if ins_count >= min_count:
                pos_aminos.append('i')
            aminos.append(pos_aminos)
        yield region, aminos
