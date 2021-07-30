from micall.utils.report_amino import ReportNucleotide, SeedAmino


def test_report_nucleotide_repr():
    report_nuc = ReportNucleotide(10)

    assert repr(report_nuc) == 'ReportNucleotide(10, SeedNucleotide({}))'


def test_nuc_consensus_index_during_add():
    seed_amino1 = SeedAmino(consensus_nuc_index=30)
    seed_amino1.count_aminos('GAG', 9)

    seed_amino2 = SeedAmino(None)
    seed_amino2.add(seed_amino1, start_nuc=1)

    assert seed_amino1.nucleotides[1].consensus_index == 31
    assert seed_amino2.nucleotides[1].consensus_index == 31