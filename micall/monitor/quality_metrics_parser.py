from struct import unpack


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
    header = data_file.read(2)  # ignore header
    version, record_length = unpack('!BB', header)
    PARSED_LENGTH = 206
    format_string = '<HHH' + 'L'*50
    if version < 4:
        raise IOError('Old quality metrics file version: {}'.format(version))
    while True:
        data = data_file.read(record_length)
        read_length = len(data)
        if read_length == 0:
            break
        if read_length < PARSED_LENGTH:
            message = 'Partial record of length {} found in quality metrics file.'.format(
                read_length)
            raise IOError(message)
        fields = unpack(format_string, data[:PARSED_LENGTH])
        yield dict(lane=fields[0],
                   tile=fields[1],
                   cycle=fields[2],
                   quality_bins=fields[3:])


def summarize_quality(records, read_lengths=None):
    """ Calculate the portion of clusters and cycles with quality >= 30.

    :param records: a sequence of dictionaries like those yielded from
    read_quality().
    :param list read_lengths: a list of lengths for each type of read: forward,
    indexes, and reverse
    :return: (fwd_q30, rev_q30) if read_lengths is not None, otherwise just
    returns a float
    """
    good_count = total_count = 0
    if read_lengths is None:
        last_forward_cycle = first_reverse_cycle = None
    else:
        good_reverse = total_reverse = 0
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

    def summarize(good, total):
        return good/float(total) if total != 0 else 0
    if read_lengths is None:
        return summarize(good_count, total_count)
    return (summarize(good_count, total_count), summarize(good_reverse, total_reverse))

if __name__ == '__live_coding__':
    import unittest
    from micall.tests.quality_metrics_parser_test import QualityMetricsParserTest

    suite = unittest.TestSuite()
    suite.addTest(QualityMetricsParserTest("test_summarize_blank"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
