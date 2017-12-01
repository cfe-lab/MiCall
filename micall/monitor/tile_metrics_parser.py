from struct import unpack

from micall.monitor.error_metrics_parser import read_records


class MetricCodes(object):
    """ Constants for metric codes used in a tile metrics data file.

    Other codes:
    (200 + (N - 1) * 2): phasing for read N
    (201 + (N - 1) * 2): prephasing for read N
    (300 + N - 1): percent aligned for read N
    """
    CLUSTER_DENSITY = 100  # K/mm2
    CLUSTER_DENSITY_PASSING_FILTERS = 101  # K/mm2
    CLUSTER_COUNT = 102
    CLUSTER_COUNT_PASSING_FILTERS = 103


def read_tiles(data_file):
    """ Read a tile metrics data file.

    :param file data_file: an open file-like object. Needs to have a two-byte
    header with the file version and the length of each record, followed by the
    records.
    :return: an iterator over the records of data in the file. Each record is a
    dictionary with the following keys:
    - lane [uint16]
    - tile [uint16]
    - metric_code [uint16]
    - metric_value [float32]
    """
    PARSED_LENGTH = 10
    format_string = '<HHHf'
    for data in read_records(data_file, min_version=2):
        fields = unpack(format_string, data[:PARSED_LENGTH])
        yield dict(lane=fields[0],
                   tile=fields[1],
                   metric_code=fields[2],
                   metric_value=fields[3])


def summarize_tile_records(records, summary):
    """ Summarize the records from a tile metrics file.

    :param records: a sequence of dictionaries from read_tiles()
    :param dict summary: a dictionary to hold the summary values:
    cluster_density and pass_rate.
    """
    density_sum = 0.0
    density_count = 0
    total_clusters = 0.0
    passing_clusters = 0.0
    for record in records:
        if record['metric_code'] == MetricCodes.CLUSTER_DENSITY:
            density_sum += record['metric_value']
            density_count += 1
        elif record['metric_code'] == MetricCodes.CLUSTER_COUNT:
            total_clusters += record['metric_value']
        elif record['metric_code'] == MetricCodes.CLUSTER_COUNT_PASSING_FILTERS:
            passing_clusters += record['metric_value']
    if density_count > 0:
        summary['cluster_density'] = density_sum/density_count
    if total_clusters > 0.0:
        summary['pass_rate'] = passing_clusters / total_clusters


def summarize_tiles(filename, summary):
    """ Summarize the records from a tile metrics file. """
    with open(filename, 'rb') as data_file:
        records = read_tiles(data_file)
        summarize_tile_records(records, summary)
