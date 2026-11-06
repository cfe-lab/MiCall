"""Regression tests for non-stealing denovo remapping and depth accounting."""
import csv
from io import StringIO

import pytest

from micall.core.denovo_multimap import select_best_pairs
from micall.core.sam2aln import sam2aln

CSV_HEADER = 'qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n'
SAM_HEADER = '@SQ\tSN:A\tLN:100\n@SQ\tSN:B\tLN:100\n'


def placement(qname, reference, pos=1, score=20, secondary=False,
              cigar='5M', seq='ACGTA', quality=None, flag1=99, flag2=147):
    quality = quality or 'I' * len(seq)
    offset = 256 if secondary else 0

    def record(flag):
        return [qname, str(flag + offset), reference, str(pos), '0', cigar,
                '=', str(pos), '0', seq, quality, f'AS:i:{score}']
    return [record(flag1), record(flag2)]


def sam(*pairs):
    return StringIO(SAM_HEADER + ''.join(
        '\t'.join(row) + '\n' for pair in pairs for row in pair))


def normalize(*pairs, limit=8):
    output = StringIO()
    counts = select_best_pairs(sam(*pairs), output, limit)
    return output.getvalue(), counts


def csv_rows(sam_text):
    return [row for row in csv.DictReader(StringIO(
        CSV_HEADER + '\n'.join(
            ','.join(line.split('\t')[:11])
            for line in sam_text.splitlines() if not line.startswith('@'))))]


def remap_csv(sam_text):
    return StringIO(CSV_HEADER + ''.join(
        ','.join(line.split('\t')[:11]) + '\n'
        for line in sam_text.splitlines() if not line.startswith('@')))


def test_two_equally_best_contigs_keep_original_names_and_mates():
    text, counts = normalize(
        placement('read1', 'A'), placement('read1', 'B', secondary=True))
    rows = csv_rows(text)
    assert [(r['qname'], r['rname']) for r in rows] == [
        ('read1', 'A'), ('read1', 'A'), ('read1', 'B'), ('read1', 'B')]
    assert counts == {'A': 2, 'B': 2}
    assert int(rows[2]['flag']) & 0x100


def test_suboptimal_secondary_placement_does_not_inflate_depth():
    text, counts = normalize(
        placement('read1', 'A', score=20),
        placement('read1', 'B', score=19, secondary=True))
    assert counts == {'A': 2}
    assert {r['rname'] for r in csv_rows(text)} == {'A'}


def test_duplicates_and_same_locus_alternatives_cannot_double_count():
    text, counts = normalize(
        placement('read1', 'A'),
        placement('read1', 'A', secondary=True),
        placement('read1', 'A', pos=2, secondary=True),
        placement('read1', 'B', secondary=True))
    assert counts == {'A': 2, 'B': 2}
    assert len(csv_rows(text)) == 4


def test_distinct_loci_on_same_contig_are_distinct_hypotheses():
    text, counts = normalize(placement('read1', 'A'),
                             placement('read1', 'A', pos=50, secondary=True))
    assert counts == {'A': 4}
    assert [r['pos'] for r in csv_rows(text)] == ['1', '1', '50', '50']


def test_invalid_supplementary_pair_excluded():
    text, counts = normalize(
        placement('read1', 'A'),
        placement('read1', 'B', secondary=True, flag1=99 + 2048,
                  flag2=147 + 2048))
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2



def test_mismatched_secondary_status_is_not_cross_paired():
    malformed = placement('read1', 'B', secondary=True)
    malformed[1][1] = '147'  # first mate secondary, second primary
    text, counts = normalize(placement('read1', 'A'), malformed)
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2


def test_mixed_primary_pair_preserves_rescue_input():
    mixed = placement('read1', 'A')
    mixed[1][2] = 'B'
    mixed[0][6] = 'B'
    mixed[1][6] = 'A'
    text, _ = normalize(mixed)
    assert [r['rname'] for r in csv_rows(text)] == ['A', 'B']



def test_discordant_same_reference_primary_preserves_legacy_alignment():
    # Bowtie2 can report an R1/R2 mapping without the proper-pair bit.
    # The previous MiCall implementation did not discard these outright.
    discordant = placement('read1', 'A', flag1=97, flag2=145)
    text, counts = normalize(discordant)
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2


def test_bad_supplementary_mate_never_resurrected_by_fallback():
    suspect = placement('read1', 'A', flag1=97, flag2=145 + 2048)
    text, counts = normalize(suspect)
    assert not counts
    assert not csv_rows(text)


def test_malformed_single_mapping_is_discarded():
    malformed = placement('read1', 'A')
    malformed[1][7] = '999'  # mate coordinate does not agree
    text, counts = normalize(malformed)
    assert not counts
    assert not csv_rows(text)



@pytest.mark.parametrize('field,value', [
    (3, '0'),          # SAM coordinates are 1-based
    (3, 'not-an-integer'),
    (7, 'not-an-integer'),
    (1, 'not-an-integer'),
    (5, '3M2Z'),       # malformed CIGAR
    (5, '2M1S'),       # CIGAR consumes fewer bases than the sequence
    (10, '!!!'),       # quality length mismatch
])
def test_invalid_primary_sam_records_fail_closed(field, value):
    invalid = placement('read1', 'A')
    invalid[0][field] = value
    text, counts = normalize(invalid)
    assert not counts
    assert not csv_rows(text)


def test_mate_orientation_flags_must_be_reciprocal():
    invalid = placement('read1', 'A')
    invalid[0][1] = '67'  # fails to flag the reverse-complement R2
    text, counts = normalize(invalid)
    assert not counts
    assert not csv_rows(text)

def test_alignment_cap_is_explicit():
    text, counts = normalize(
        placement('read1', 'A'),
        placement('read1', 'B', secondary=True), limit=1)
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2


def test_avg_divides_each_read_pair_by_its_own_number_of_placements():
    text, _ = normalize(
        placement('read1', 'A'), placement('read1', 'B', secondary=True),
        placement('read2', 'A'))
    depth = StringIO()
    aligned = StringIO()
    sam2aln(remap_csv(text), aligned, depth_csv=depth)
    rows = list(csv.DictReader(StringIO(depth.getvalue())))
    assert len(rows) == 10
    by_ref_pos = {(r['refname'], int(r['position'])): r for r in rows}
    assert by_ref_pos['A', 1]['max_depth'] == '2'
    assert float(by_ref_pos['A', 1]['avg_depth']) == 1.5
    assert by_ref_pos['B', 1]['max_depth'] == '1'
    assert float(by_ref_pos['B', 1]['avg_depth']) == 0.5
    assert all(by_ref_pos['A', i]['max_depth'] == '2' for i in range(1, 6))


def test_ambiguous_alignment_does_not_change_single_mapping_format():
    text, _ = normalize(placement('read1', 'A'))
    aligned = StringIO()
    depth = StringIO()
    sam2aln(remap_csv(text), aligned, depth_csv=depth)
    assert aligned.getvalue().splitlines()[1] == 'A,15,0,1,0,ACGTA'
    assert all(r['max_depth'] == '1' and r['avg_depth'] == '1'
               for r in csv.DictReader(StringIO(depth.getvalue())))


def test_bad_quality_not_in_any_count_or_avg_denominator():
    text, _ = normalize(
        placement('read1', 'A', quality='!!!!!'),
        placement('read1', 'B', secondary=True, quality='IIIII'))
    depth = StringIO()
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth)
    rows = list(csv.DictReader(StringIO(depth.getvalue())))
    assert {(r['refname'], r['avg_depth']) for r in rows} == {('B', '1')}


@pytest.mark.parametrize('cigar,seq,first,last', [
    ('2M1I2M', 'ACGTA', 3, 6),
    ('1S4M', 'ACGTA', 3, 6),
])
def test_insertions_and_softclips_do_not_create_reference_positions(
        cigar, seq, first, last):
    text, _ = normalize(placement('read1', 'A', pos=first,
                                  cigar=cigar, seq=seq))
    depth = StringIO()
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth)
    rows = list(csv.DictReader(StringIO(depth.getvalue())))
    assert [int(r['position']) for r in rows] == list(range(first, last + 1))


def test_reference_deletion_not_counted_as_covered_base():
    text, _ = normalize(placement('read1', 'A', cigar='3M1D2M'))
    depth = StringIO()
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth)
    positions = [int(r['position']) for r in
                 csv.DictReader(StringIO(depth.getvalue()))]
    assert positions == [1, 2, 3, 5, 6]


def test_high_offset_depth_uses_absolute_reference_coordinates():
    # Avoid regression if temporary paired-read strings use relative offsets.
    text, _ = normalize(placement('read1', 'A', pos=5000))
    depth = StringIO()
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth)
    assert [int(r['position']) for r in
            csv.DictReader(StringIO(depth.getvalue()))] == list(range(5000, 5005))


def test_zero_coverage_reference_positions_are_reported():
    text, _ = normalize(placement('read1', 'A', pos=3))
    depth = StringIO()
    reference = StringIO('region,sequence\nA,ACGTACGTAC\nB,ACGTACGTAC\n')
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth,
            depth_reference_csv=reference)
    rows = list(csv.DictReader(StringIO(depth.getvalue())))
    assert len(rows) == 20
    assert [(r['max_depth'], r['avg_depth']) for r in rows[:3]] == [
        ('0', '0'), ('0', '0'), ('1', '1')]
    assert all(r['max_depth'] == '0' and r['avg_depth'] == '0'
               for r in rows[10:])


def test_non_numeric_alignment_score_is_not_silently_accepted():
    malformed = placement('read1', 'B', secondary=True)
    malformed[0][-1] = 'AS:i:bad'
    text, counts = normalize(placement('read1', 'A'), malformed)
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2


def test_only_softclipping_is_not_a_valid_aligned_placement():
    invalid = placement('read1', 'A', cigar='5S')
    text, counts = normalize(invalid)
    assert not counts
    assert not csv_rows(text)


def test_empty_or_incomplete_sam_lines_cannot_crash_selector():
    text = StringIO(SAM_HEADER + 'broken\t99\tA\n'
                    'broken\t147\tA\n')
    result = StringIO()
    assert not select_best_pairs(text, result)
    assert not csv_rows(result.getvalue())


@pytest.mark.parametrize('alternatives', range(1, 9))
def test_fractional_depth_conserves_one_pair_across_all_placements(alternatives):
    from micall.core.denovo_multimap import MAX_DENOVO_ALIGNMENTS

    assert alternatives <= MAX_DENOVO_ALIGNMENTS
    references = [f'R{i}' for i in range(alternatives)]
    pairs = [placement('read1', reference, secondary=i > 0)
             for i, reference in enumerate(references)]
    text, counts = normalize(*pairs)
    depth = StringIO()
    sam2aln(remap_csv(text), StringIO(), depth_csv=depth)
    rows = list(csv.DictReader(StringIO(depth.getvalue())))
    assert counts == {reference: 2 for reference in references}
    assert len(rows) == 5 * alternatives
    assert all(int(r['max_depth']) == 1 for r in rows)
    for pos in range(1, 6):
        weights = [float(r['avg_depth']) for r in rows
                   if int(r['position']) == pos]
        assert sum(weights) == pytest.approx(1.0)


def test_reporting_cap_is_auditable_not_mistaken_for_exhaustiveness():
    from micall.core.denovo_multimap import MAX_DENOVO_ALIGNMENTS

    pairs = [placement('read1', f'R{i}', secondary=i != 0)
             for i in range(MAX_DENOVO_ALIGNMENTS)]
    diagnostics = {}
    out = StringIO()
    select_best_pairs(sam(*pairs), out, diagnostics=diagnostics)
    assert diagnostics['at_report_limit'] == 1


def test_discarded_mapped_pair_is_reported_in_diagnostics():
    suspect = placement('read1', 'A')
    suspect[0][3] = '0'
    diagnostics = {}
    out = StringIO()
    select_best_pairs(sam(suspect), out, diagnostics=diagnostics)
    assert diagnostics['discarded_mapped_pairs'] == 1


def test_unmapped_pair_does_not_trigger_malformed_mapped_warning():
    # Both mates are unmapped and should be dropped without alarming users.
    unmapped = placement('read1', '*', flag1=77, flag2=141)
    unmapped[0][5] = unmapped[1][5] = '*'
    diagnostics = {}
    out = StringIO()
    select_best_pairs(sam(unmapped), out, diagnostics=diagnostics)
    assert 'discarded_mapped_pairs' not in diagnostics
