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
