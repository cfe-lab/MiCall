import pytest

from micall.core.consensus_builder import ConsensusBuilder


def prepare_reads(*merged_sequences):
    ignored = None
    for merged_sequence in merged_sequences:
        yield (ignored,  # name
               ignored,  # read 1
               ignored,  # read 2
               merged_sequence)


def test_simple():
    merged_reads = list(prepare_reads("AAACCCTTTGGGAAACCC",
                                      "ATACCCTTTGGGAAACCC",
                                      "AAACCCTTTGGGAAACCC"))
    expected_consensus = "AAACCCTTTGGGAAACCC"
    builder = ConsensusBuilder()

    returned_reads = list(builder.build(merged_reads))

    assert merged_reads == returned_reads
    assert expected_consensus == builder.get_consensus()


def test_different_lengths():
    merged_reads = prepare_reads("AAACCCTTTGGGAAACCC",
                                 "ATACCCTTTGGGAAACCC",
                                 "AAACCCTTTGGGAAACCC",
                                 "CATGAGACATCACACAC",
                                 "CATGAGACATCACACA",
                                 "CATGAGACATCACAC",
                                 "CATGAGACATCACA")
    expected_consensus = "AAACCCTTTGGGAAACCC"
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    assert expected_consensus == builder.get_consensus()


def test_unmerged_pair():
    merged_reads = list(prepare_reads("AAACCCTTTGGGAAACCC",
                                      "ATACCCTTTGGGAAACCC",
                                      None,  # unmerged pair
                                      "AAACCCTTTGGGAAACCC"))
    expected_consensus = "AAACCCTTTGGGAAACCC"
    builder = ConsensusBuilder()

    returned_reads = list(builder.build(merged_reads))

    assert merged_reads == returned_reads
    assert expected_consensus == builder.get_consensus()


def test_no_reads():
    merged_reads = prepare_reads()
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    with pytest.raises(IndexError):
        builder.get_consensus()


def test_build_by_lengths():
    raw_reads = []
    for n in range(100, 141):
        if n == 120:
            raw_reads += ['A'*n] * 100
        else:
            raw_reads += ['C'*n] * 5
    for n in range(200, 241):
        if n == 220:
            raw_reads += ['G'*n] * 1000
        else:
            raw_reads += ['T'*n] * 50
    merged_reads = prepare_reads(*raw_reads)
    expected_consensuses = ['A'*120, 'G'*220]
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    consensuses = list(builder.get_consensus_by_lengths())
    assert expected_consensuses == consensuses


def test_build_by_lengths_when_spike_too_short():
    raw_reads = []
    for n in range(100, 141):
        if n == 120:
            raw_reads += ['A'*n] * 99
        else:
            raw_reads += ['C'*n] * 5
    merged_reads = prepare_reads(*raw_reads)
    expected_consensuses = []
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    consensuses = list(builder.get_consensus_by_lengths())
    assert expected_consensuses == consensuses


def test_build_by_lengths_with_no_neighbours():
    raw_reads = ['A'*120] * 19 + ['G'*220] * 20
    merged_reads = prepare_reads(*raw_reads)
    expected_consensuses = ['G'*220]
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    consensuses = list(builder.get_consensus_by_lengths())
    assert expected_consensuses == consensuses


def test_build_by_lengths_no_reads():
    merged_reads = prepare_reads()
    expected_consensuses = []
    builder = ConsensusBuilder()

    list(builder.build(merged_reads))

    consensuses = list(builder.get_consensus_by_lengths())
    assert expected_consensuses == consensuses
