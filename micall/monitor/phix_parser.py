from struct import unpack
import csv
import os
from operator import itemgetter
import sys
import math
from itertools import groupby


def read_phix(data_file):
    """ Read error rate data from a phiX data file.

    :param file data_file: an open file-like object. Needs to have a two-byte
    header with the file version and the length of each record, followed by the
    records.
    :return: an iterator over the records of data in the file. Each record is a
    dictionary with the following keys:
    - lane [uint16]
    - tile [uint16]
    - cycle [uint16]
    - error_rate [float]
    - num_0_errors [uint32]
    - num_1_error [uint32]
    - num_2_errors [uint32]
    - num_3_errors [uint32]
    - num_4_errors [uint32]
    """
    header = data_file.read(2)  # ignore header
    version, record_length = unpack('!bb', header)
    PARSED_LENGTH = 30
    if version < 3:
        raise IOError('Old phiX error file version: {}'.format(version))
    while True:
        data = data_file.read(record_length)
        read_length = len(data)
        if read_length == 0:
            break
        if read_length < PARSED_LENGTH:
            message = 'Partial record of length {} found in phiX error file.'.format(
                read_length)
            raise IOError(message)
        fields = unpack('<HHHfLLLLL', data[:PARSED_LENGTH])
        yield dict(lane=fields[0],
                   tile=fields[1],
                   cycle=fields[2],
                   error_rate=fields[3],
                   num_0_errors=fields[4],
                   num_1_error=fields[5],
                   num_2_errors=fields[6],
                   num_3_errors=fields[7],
                   num_4_errors=fields[8])


def _yield_cycles(records, read_lengths):
    sorted_records = map(itemgetter('tile', 'cycle', 'error_rate'), records)
    sorted_records.sort()
    max_forward_cycle = read_lengths and read_lengths[0] or sys.maxint
    min_reverse_cycle = read_lengths and sum(read_lengths[:-1])+1 or sys.maxint
    for record in sorted_records:
        cycle = record[1]
        if cycle >= min_reverse_cycle:
            cycle = min_reverse_cycle - cycle - 1
            yield record[0], cycle, record[2]
        elif cycle <= max_forward_cycle:
            yield record


def _record_grouper(record):
    # Group by tile and sign of cycle (forward or reverse).
    return (record[0], int(math.copysign(1, record[1])))


def write_phix_csv(out_file, records, read_lengths=None):
    """ Write phiX error rate data to a comma-separated-values file.

    Missing cycles are written with blank error rates, index reads are not
    written, and reverse reads are written with negative cycles.
    :param out_file: an open file to write to
    :param records: a sequence of dictionaries like those yielded from
    read_phix().
    :param read_lengths: a list of lengths for each type of read: forward,
    indexes, and reverse
    """
    writer = csv.writer(out_file, lineterminator=os.linesep)
    writer.writerow(['tile', 'cycle', 'errorrate'])

    for (_tile, sign), group in groupby(_yield_cycles(records, read_lengths),
                                        _record_grouper):
        previous_cycle = 0
        for record in group:
            cycle = record[1]
            previous_cycle += sign
            while previous_cycle*sign < cycle*sign:
                writer.writerow((record[0], previous_cycle))
                previous_cycle += sign
            writer.writerow(record)
        if read_lengths:
            read_length = read_lengths[0] if sign == 1 else -read_lengths[-1]
            while previous_cycle*sign < read_length*sign:
                previous_cycle += sign
                writer.writerow((record[0], previous_cycle))

if __name__ == '__main__':
    in_path = ('/home/don/data/RAW_DATA/MiSeq/runs/140123_M01841_microtest/'
               'InterOp/ErrorMetricsOut.bin')
    out_path = '/home/don/data/miseq/140123_M01841_microtest/quality2.csv'
    with open(in_path, 'rb') as in_file, open(out_path, 'wb') as out_file:
        records = read_phix(in_file)
        write_phix_csv(out_file, records, read_lengths=[251, 8, 8, 251])
    print('Done.')
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.phix_parser_test import PhixParserTest

    suite = unittest.TestSuite()
    suite.addTest(PhixParserTest("test_write_missing_end"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
