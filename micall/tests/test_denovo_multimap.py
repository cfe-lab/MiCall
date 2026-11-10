"""Regression tests for non-stealing denovo remapping and depth accounting."""
import csv
import random
from io import StringIO
from itertools import pairwise

import pytest

from micall.core.denovo_multimap import (
    select_best_pairs)
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
    # Identical duplicate placements are removed, but a shifted placement
    # must retain its previously uncovered final reference position.
    assert counts == {'A': 4, 'B': 2}
    assert len(csv_rows(text)) == 6
    aligned, depth = StringIO(), StringIO()
    sam2aln(remap_csv(text), aligned, depth_csv=depth)
    by_position = {(r['refname'], int(r['position'])): r
                   for r in csv.DictReader(StringIO(depth.getvalue()))}
    assert {pos for ref, pos in by_position if ref == 'A'} == set(range(1, 7))
    assert all(r['max_depth'] == '1' for r in by_position.values())
    assert float(by_position['A', 2]['avg_depth']) == pytest.approx(2/3)


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


def test_primary_singleton_survives_without_claiming_a_complete_pair():
    primary = placement('orphan', 'A', flag1=73, flag2=133)
    primary[1][5] = '*'  # R2 is unmapped; R1 is valid mapped evidence.
    normalized, counts = normalize(primary)
    assert counts == {'A': 1}
    rows = csv_rows(normalized)
    assert [int(r['flag']) for r in rows] == [73, 133]
    aligned = StringIO()
    depth = StringIO()
    failures = StringIO()
    sam2aln(remap_csv(normalized), aligned, failed_csv=failures,
            depth_csv=depth)
    # Existing MiCall requires both CIGARs to merge; keeping singleton
    # records maintains remap/cascade parity, not false Q15 depth.
    assert list(csv.DictReader(StringIO(aligned.getvalue()))) == []
    assert list(csv.DictReader(StringIO(depth.getvalue()))) == []
    assert list(csv.DictReader(StringIO(failures.getvalue()))) == [
        {'qname': 'orphan', 'cause': 'badCigar'}]


def test_unmapped_first_mate_primary_is_sorted_before_mapped_second():
    # This ordering occurred in real MiCall HIV/HCV source reads.
    primary = placement('reverse-orphan', 'A', flag1=69, flag2=153)
    primary[0][5] = '*'
    normalized, counts = normalize(list(reversed(primary)))
    assert counts == {'A': 1}
    assert [int(r['flag']) for r in csv_rows(normalized)] == [69, 153]


def test_primary_singleton_rejects_malformed_evidence():
    primary = placement('orphan', 'A', flag1=73, flag2=133)
    primary[1][5] = '*'
    invalid = [list(r) for r in primary]
    invalid[0][5] = '2M1S'  # Claims only three bases consumed; read has five.
    assert normalize(invalid)[1] == {}
    invalid = [list(r) for r in primary]
    invalid[0][1] = str(73 | 0x800)  # Supplementary is not primary evidence.
    assert normalize(invalid)[1] == {}
    invalid = [list(r) for r in primary]
    invalid[1][5] = '5M'  # An unmapped mate cannot have an aligned CIGAR.
    assert normalize(invalid)[1] == {}


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


def test_overlap_choice_does_not_depend_on_bowtie_reporting_order():
    left = placement('read1', 'A', pos=1)
    right = placement('read1', 'A', pos=3, secondary=True)
    first, _ = normalize(left, right)
    second, _ = normalize(right, left)
    assert [row['pos'] for row in csv_rows(first)] == ['1', '1', '3', '3']
    assert [row['pos'] for row in csv_rows(second)] == ['1', '1', '3', '3']


def test_callback_counts_unique_pairs_not_bowtie_multi_alignment_rows(tmp_path):
    from collections import Counter
    from unittest.mock import Mock
    from micall.core.remap import map_to_reference

    records = [placement(f'read{i}', f'contig{j}', secondary=j > 0)
               for i in range(80) for j in range(8)]
    lines = ['\t'.join(row) + '\n' for pair in records for row in pair]
    mapper = Mock()
    mapper.yield_output.return_value = lines
    build = Mock()
    progress = Mock()
    raw_reads = 160
    map_to_reference('R1.fastq', 'R2.fastq', {'contig0': 'ACGTA'},
                     str(tmp_path / 'ref.fa'), str(tmp_path / 'reads.sam'),
                     StringIO(), StringIO(), mapper, build, raw_reads,
                     15, 15, Counter(), StringIO(), progress,
                     max_alignments=8)
    values = [call.kwargs['progress'] for call in progress.call_args_list]
    assert values[-1] == raw_reads
    assert max(values) <= raw_reads
    assert all(a <= b for a, b in pairwise(values))
    assert values[1] < raw_reads


def test_missing_as_tag_does_not_beat_negative_scored_alignment():
    scored = placement('read1', 'A', score=-10)
    unscored = placement('read1', 'B', secondary=True)
    for row in unscored:
        row.pop()  # no AS tag: must not silently outrank -10
    text, counts = normalize(scored, unscored)
    assert counts == {'A': 2}
    assert len(csv_rows(text)) == 2


def test_concordant_pairs_without_alignment_scores_are_not_selected():
    first = placement('read1', 'A')
    second = placement('read1', 'B', secondary=True)
    for row in (*first, *second):
        row.pop()
    text, counts = normalize(first, second)
    assert not counts
    assert not csv_rows(text)


def test_maximum_depth_exactly_matches_accepted_q15_consensus_input():
    """No independent rematching logic may diverge from the aligned CSV."""
    pairs = [
        placement('read1', 'A'),
        placement('read1', 'B', pos=4, secondary=True),
        placement('read2', 'A', pos=3, quality='II!II'),
        placement('read3', 'A', pos=8, cigar='3M1D2M'),
        placement('read4', 'B', pos=14, cigar='1S4M'),
    ]
    text, _ = normalize(*pairs)
    aligned, depths = StringIO(), StringIO()
    sam2aln(remap_csv(text), aligned, depth_csv=depths)

    expected = {}
    for row in csv.DictReader(StringIO(aligned.getvalue())):
        if row['qcut'] != '15':
            continue
        reference, offset = row['refname'], int(row['offset'])
        for i, base in enumerate(row['seq'], start=offset + 1):
            if base in 'ACGT':
                key = (reference, i)
                expected[key] = expected.get(key, 0) + int(row['count'])
    actual = {(row['refname'], int(row['position'])): int(row['max_depth'])
              for row in csv.DictReader(StringIO(depths.getvalue()))}
    assert actual == expected


@pytest.mark.parametrize('seed', range(24))
def test_deterministic_fuzzed_placements_match_q15_depth(seed):
    rng = random.Random(seed)
    candidates = []
    for i in range(rng.randint(1, 8)):
        region = rng.choice(['A', 'B', 'C'])
        position = rng.randint(1, 35)
        score = rng.choice([20, 20, 20, 19])
        candidates.append(placement(
            'fragment', region, pos=position, score=score,
            secondary=i > 0))
    original = tuple(candidates)
    shuffled = list(candidates)
    rng.shuffle(shuffled)

    def analyze(pairs):
        text, counts = normalize(*pairs)
        aligned, depths = StringIO(), StringIO()
        sam2aln(remap_csv(text), aligned, depth_csv=depths)
        physical = {}
        for row in csv.DictReader(StringIO(aligned.getvalue())):
            for pos, base in enumerate(row['seq'],
                                       start=int(row['offset']) + 1):
                if base in 'ACGT':
                    key = (row['refname'], pos)
                    physical[key] = physical.get(key, 0) + 1
        max_depths = {(r['refname'], int(r['position'])):
                      int(r['max_depth'])
                      for r in csv.DictReader(StringIO(depths.getvalue()))}
        return counts, physical, max_depths

    expected = analyze(original)
    assert expected == analyze(shuffled)
    assert expected[1] == expected[2]



def test_orphaned_secondary_sam_record_is_reported():
    from micall.core.denovo_multimap import select_best_pairs

    valid = placement('read1', 'A')
    orphan = placement('read1', 'B', secondary=True)[0]
    content = sam(valid)
    content.seek(0, 2)
    content.write('\t'.join(orphan) + '\n')
    content.seek(0)
    output = StringIO()
    diagnostics = {}
    select_best_pairs(content, output, diagnostics=diagnostics)
    assert diagnostics['incomplete_pair_groups'] == 1
    assert len(csv_rows(output.getvalue())) == 2


def test_unbounded_selector_retains_all_nine_equivalent_contigs():
    pairs = [placement('read1', f'R{i}', secondary=i > 0)
             for i in range(9)]
    output = StringIO()
    diagnostics = {}
    counts = select_best_pairs(sam(*pairs), output,
                               limit=None, diagnostics=diagnostics)
    assert len(csv_rows(output.getvalue())) == 18
    assert len(counts) == 9
    assert diagnostics.get('at_report_limit', 0) == 0


def test_unbounded_sam_reader_rejects_pathologically_repetitive_pair():
    from micall.core.denovo_multimap import MAX_DENOVO_RECORDS_PER_PAIR

    source = sam(*(placement('huge', f'R{i}', secondary=i != 0)
                   for i in range(MAX_DENOVO_RECORDS_PER_PAIR // 2 + 1)))
    out = StringIO()
    diagnostics = {}
    counts = select_best_pairs(source, out, limit=None,
                               diagnostics=diagnostics)
    assert not counts
    assert not csv_rows(out.getvalue())
    assert diagnostics['oversized_pair_groups'] == 1


def test_unequal_bowtie2_discordant_mate_alternatives_are_not_malformed():
    # In Bowtie2 -a mode, a discordant read may have three R1 placements
    # and only one R2 placement; adjacent record pairing is not meaningful.
    primary = placement('discordant', 'A', flag1=97, flag2=145)
    primary[1][2] = 'B'
    primary[0][6] = 'B'
    primary[1][6] = 'A'
    secondary_a = primary[0].copy()
    secondary_b = primary[0].copy()
    secondary_a[1] = str(int(secondary_a[1]) + 256)
    secondary_a[2] = 'C'
    secondary_b[1] = str(int(secondary_b[1]) + 256)
    secondary_b[2] = 'D'
    raw = sam([primary[0], secondary_a, secondary_b, primary[1]])
    out = StringIO()
    diagnostics = {}
    select_best_pairs(raw, out, limit=None, diagnostics=diagnostics)
    observed = csv_rows(out.getvalue())
    assert [(r['rname'], r['flag']) for r in observed] == [
        ('A', '97'), ('B', '145')]
    assert diagnostics.get('incomplete_pair_groups', 0) == 0


def test_partially_overlapping_alternatives_keep_distant_mate_coverage():
    # R1 overlaps between alternatives, but R2 is moved to a distant locus.
    # The molecule supports the shared R1 reference positions only once and
    # still contributes to the distant R2 locus in the second hypothesis.
    near = placement('fragment', 'A', pos=100)
    distant = placement('fragment', 'A', pos=101, secondary=True)
    distant[0][7] = '800'
    distant[1][3] = '800'
    distant[1][7] = '101'
    text, counts = normalize(near, distant)
    assert counts == {'A': 4}
    aligned, depth = StringIO(), StringIO()
    sam2aln(remap_csv(text), aligned, depth_csv=depth)
    coverage = {(int(r['position'])):
                (int(r['max_depth']), float(r['avg_depth']))
                for r in csv.DictReader(StringIO(depth.getvalue()))}
    assert coverage[101] == (1, 1.0)
    assert coverage[800] == (1, 0.5)
    assert coverage[804] == (1, 0.5)
    assert set(coverage) == set(range(100, 106)) | set(range(800, 805))
    from_aligned = {}
    for row in csv.DictReader(StringIO(aligned.getvalue())):
        for pos, base in enumerate(row['seq'], start=int(row['offset']) + 1):
            if base in 'ACGT':
                from_aligned[pos] = from_aligned.get(pos, 0) + 1
    assert from_aligned == {pos: max_depth
                            for pos, (max_depth, _) in coverage.items()}
    # Internal masked overlaps are absent coverage, NOT deletion calls.
    # SeedNucleotide counts '-' as deletion, but ignores lowercase 'n'.
    from micall.utils.report_amino import SeedNucleotide
    nuc = SeedNucleotide()
    seqs = [r['seq'] for r in
            csv.DictReader(StringIO(aligned.getvalue()))]
    assert any('n' in seq for seq in seqs)
    for seq in seqs:
        for base in seq:
            nuc.count_nucleotides(base)
    assert nuc.counts['-'] == 0


def test_every_denovo_search_uses_finite_bowtie_reporting_limit(tmp_path):
    from collections import Counter
    from unittest.mock import Mock
    from micall.core.remap import map_to_reference

    records = [placement('one', 'A')]
    mapper = Mock()
    mapper.yield_output.return_value = [
        '\t'.join(row) + '\n' for pair in records for row in pair]
    build = Mock()
    from micall.core.denovo_multimap import MAX_DENOVO_ALIGNMENTS
    limit = MAX_DENOVO_ALIGNMENTS
    map_to_reference('R1.fastq', 'R2.fastq', {'A': 'ACGTA'},
                     str(tmp_path / 'ref.fa'), str(tmp_path / 'reads.sam'),
                     StringIO(), StringIO(), mapper, build, 2,
                     15, 15, Counter(), StringIO(), None,
                     max_alignments=limit)
    args = mapper.yield_output.call_args.args[0]
    assert '-a' not in args
    assert args[args.index('-k') + 1] == '12'


@pytest.mark.parametrize('deletion_first', [True, False])
def test_alternative_deletion_does_not_override_nucleotide_support(
        deletion_first):
    # The second placement reaches new positions, so the overlap must be
    # retained without introducing a deletion at position 104.
    nucleotide = placement('fragment', 'A', pos=100)
    deletion = placement('fragment', 'A', pos=101, secondary=True,
                         cigar='3M1D2M')
    pairs = (deletion, nucleotide) if deletion_first else (nucleotide, deletion)
    aligned, depth = StringIO(), StringIO()
    sam2aln(remap_csv(sam(*pairs).getvalue()), aligned, depth_csv=depth)

    rows = list(csv.DictReader(StringIO(aligned.getvalue())))
    assert len(rows) == 2
    assert all('-' not in row['seq'] for row in rows)
    covered = set()
    for row in rows:
        for pos, base in enumerate(row['seq'], int(row['offset']) + 1):
            if base in 'ACGT':
                assert pos not in covered
                covered.add(pos)
    assert covered == set(range(100, 107))
    depths = {(row['refname'], int(row['position'])): row
              for row in csv.DictReader(StringIO(depth.getvalue()))}
    assert depths['A', 104]['max_depth'] == '1'


def test_shared_deletion_position_is_counted_once_across_placements():
    first = placement('fragment', 'A', pos=100, cigar='3M1D2M')
    second = placement('fragment', 'A', pos=101, cigar='2M1D3M',
                       secondary=True)
    aligned = StringIO()
    sam2aln(remap_csv(sam(first, second).getvalue()), aligned)
    rows = list(csv.DictReader(StringIO(aligned.getvalue())))
    assert len(rows) == 2
    assert sum(row['seq'].count('-') for row in rows) == 1
