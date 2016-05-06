from struct import unpack

from micall.monitor.error_metrics_parser import read_records


def read_quality(data_file):
    """ Read a quality metrics data file.

    :param file data_file: an open file-like object. Needs to have a two-byte
    header with the file version and the length of each record, followed by the
    records.
    :return: an iterator over the records of data in the file. Each record is a
    dictionary with the following keys:
    - lane [uint16]
    - tile [uint16]
    - cycle [uint16]
    - quality_bins [list of 50 uint32, representing quality 1 to 50]
    """
    PARSED_LENGTH = 206
    format_string = '<HHH' + 'L'*50
    for data in read_records(data_file, min_version=4):
        fields = unpack(format_string, data[:PARSED_LENGTH])
        yield dict(lane=fields[0],
                   tile=fields[1],
                   cycle=fields[2],
                   quality_bins=fields[3:])


def summarize_quality_records(records, summary, read_lengths=None):
    """ Calculate the portion of clusters and cycles with quality >= 30.

    :param records: a sequence of dictionaries like those yielded from
    read_quality().
    :param dict summary: a dictionary to hold the summary values:
    q30_fwd and q30_rev. If read_lengths is None, only fwd_q30 will be set.
    :param list read_lengths: a list of lengths for each type of read: forward,
    indexes, and reverse
    """
    good_count = total_count = 0
    good_reverse = total_reverse = 0
    if read_lengths is None:
        last_forward_cycle = first_reverse_cycle = None
    else:
        last_forward_cycle = read_lengths[0]
        first_reverse_cycle = sum(read_lengths[:-1]) + 1
    for record in records:
        cycle = record['cycle']
        cycle_clusters = sum(record['quality_bins'])
        cycle_good = sum(record['quality_bins'][29:])

        if last_forward_cycle is None or cycle <= last_forward_cycle:
            total_count += cycle_clusters
            good_count += cycle_good
        elif cycle >= first_reverse_cycle:
            total_reverse += cycle_clusters
            good_reverse += cycle_good

    if total_count > 0:
        summary['q30_fwd'] = good_count/float(total_count)
    if total_reverse > 0:
        summary['q30_rev'] = good_reverse/float(total_reverse)


def summarize_quality(filename, summary, read_lengths=None):
    """ Summarize the records from a quality metrics file. """
    with open(filename, 'rb') as data_file:
        records = read_quality(data_file)
        summarize_quality_records(records, summary, read_lengths)

if __name__ == '__live_coding__':
    import unittest
    from micall.tests.quality_metrics_parser_test import QualityMetricsParserTest

    suite = unittest.TestSuite()
    suite.addTest(QualityMetricsParserTest("test_summarize_blank"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
