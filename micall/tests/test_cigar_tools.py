import pytest
from typing import List, Tuple
from math import floor
import itertools

from micall.utils.consensus_aligner import CigarActions
from micall.utils.cigar_tools import Cigar, CigarHit, parse_cigar_operation, CIGAR_OP_MAPPING


cigar_mapping_cases: List[Tuple[Cigar, 'mapping', 'closest_mapping']] = [
    # Simple cases
    ('3M',     {0: 0, 1: 1, 2: 2},
               {0: 0, 1: 1, 2: 2}),
    ('1M1D1M', {0: 0, 2: 1},
               {0: 0, 1: 0, 2: 1}),
    ('1M1I1M', {0: 0, 1: 2},
               {0: 0, 1: 2}),
    ('2M2D2M', {0: 0, 1: 1, 4: 2, 5: 3},
               {0: 0, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3}),
    ('2M2I2M', {0: 0, 1: 1, 2: 4, 3: 5},
               {0: 0, 1: 1, 2: 4, 3: 5}),
    ('3M1D3M', {0: 0, 1: 1, 2: 2, 4: 3, 5: 4, 6: 5},
               {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4, 6: 5}),
    ('3M1I3M', {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 6},
               {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 6}),
    ('7M1I3M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 9, 9: 10},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 9, 9: 10}),
    ('5M2D4M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 7: 5, 8: 6, 9: 7, 10: 8},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 4, 6: 5, 7: 5, 8: 6, 9: 7, 10: 8}),
    ('5M3I4M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 8, 6: 9, 7: 10, 8: 11},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 8, 6: 9, 7: 10, 8: 11}),
    ('1M1D',   {0: 0},
               {0: 0, 1: 0}),
    ('1M1I',   {0: 0},
               {0: 0}),
    ('1I1M',   {0: 1},
               {0: 1}),
    ('1D1M',   {1: 0},
               {1: 0, 0: 0}),

    # Multiple deletions and insertions
    ('2M2D2M2I2M', {0: 0, 1: 1, 4: 2, 5: 3, 6: 6, 7: 7},
               {0: 0, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 6, 7: 7}),
    ('2M2I2M2D2M', {0: 0, 1: 1, 2: 4, 3: 5, 6: 6, 7: 7},
               {0: 0, 1: 1, 2: 4, 3: 5, 4: 5, 5: 6, 6: 6, 7: 7}),
    ('2=1X2N1N2=1H2S', {0: 0, 1: 1, 2: 2, 6: 3, 7: 4},
               {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 3, 6: 3, 7: 4}),
    ('2M2D2M2I2M', {0: 0, 1: 1, 4: 2, 5: 3, 6: 6, 7: 7},
               {0: 0, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 6, 7: 7}),
    ('3=1X2N1N2=1H2S', {0: 0, 1: 1, 2: 2, 3: 3, 7: 4, 8: 5},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 3, 6: 4, 7: 4, 8: 5}),

    # Edge cases
    ('', {}, {}),
    ('12I', {}, {}),
    ('12D', {}, ValueError()),
]


@pytest.mark.parametrize("cigar_str, expected_mapping", [(x[0], x[1]) for x in cigar_mapping_cases])
def test_cigar_to_coordinate_mapping(cigar_str, expected_mapping):
    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    assert expected_mapping == mapping.ref_to_query_d
    assert expected_mapping == {i: mapping.ref_to_query(i)
                                for i in mapping.mapped_reference_coordinates()}


@pytest.mark.parametrize("cigar_str", [x[0] for x in cigar_mapping_cases])
def test_cigar_to_coordinate_bijection_property(cigar_str):
    inverse = lambda d: {v: k for k, v in d.items()}

    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    assert mapping.query_to_ref_d == inverse(mapping.ref_to_query_d)
    assert mapping.ref_to_query_d == inverse(mapping.query_to_ref_d)
    assert mapping.ref_to_query_d == inverse(inverse(mapping.ref_to_query_d))
    assert mapping.query_to_ref_d == inverse(inverse(mapping.query_to_ref_d))


@pytest.mark.parametrize("cigar_str, expected_closest_mapping", [(x[0], x[2]) for x in cigar_mapping_cases])
def test_cigar_to_closest_coordinate_mapping(cigar_str, expected_closest_mapping):
    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    def test():
        fullrange = {i: mapping.ref_to_closest_query(i)
                     for i in mapping.all_reference_coordinates()}
        assert expected_closest_mapping == fullrange

    if isinstance(expected_closest_mapping, Exception):
        with pytest.raises(type(expected_closest_mapping)):
            test()
    else:
        test()


@pytest.mark.parametrize("cigar_str, expected_mapping", [(x[0], x[1]) for x in cigar_mapping_cases])
def test_cigar_hit_to_coordinate_mapping(cigar_str, expected_mapping):
    cigar = Cigar.coerce(cigar_str)
    hit = CigarHit(cigar, r_st=5, r_ei=(5 + cigar.ref_length - 1), q_st=7, q_ei=(7 + cigar.query_length - 1))
    mapping = hit.coordinate_mapping

    # Coordinates are translated by q_st and r_st.
    expected_mapping = {k + hit.r_st: v + hit.q_st for (k, v) in expected_mapping.items()}
    assert mapping.ref_to_query(0) == None
    assert mapping.query_to_ref(0) == None
    assert expected_mapping \
        == {i: mapping.ref_to_query(i)
            for i in mapping.mapped_reference_coordinates()}


@pytest.mark.parametrize("cigar_str, expected_closest_mapping", [(x[0], x[2]) for x in cigar_mapping_cases])
def test_cigar_hit_to_coordinate_closest_mapping(cigar_str, expected_closest_mapping):
    cigar = Cigar.coerce(cigar_str)
    hit = CigarHit(cigar, r_st=5, r_ei=(5 + cigar.ref_length - 1), q_st=7, q_ei=(7 + cigar.query_length - 1))
    mapping = hit.coordinate_mapping

    def test(expected):
        # Coordinates are translated by q_st and r_st.
        fullrange = {i: mapping.ref_to_closest_query(i)
                     for i in mapping.all_reference_coordinates()}
        assert expected == fullrange

    if isinstance(expected_closest_mapping, Exception):
        with pytest.raises(type(expected_closest_mapping)):
            test(expected_closest_mapping)
    else:
        test({k + hit.r_st: v + hit.q_st for (k, v) in expected_closest_mapping.items()})


def test_invalid_operation_in_cigar_string():
    with pytest.raises(ValueError):
        Cigar.coerce('3M1Z3M') # Z operation is not implemented


def test_invalid_operation_in_cigar_list():
    with pytest.raises(ValueError):
        Cigar.coerce([(3, 42)]) # Operation code "42" does not exist


def test_invalid_cigar_string():
    with pytest.raises(ValueError):
        Cigar.coerce('3MMMMMM3M') # Too many Ms
    with pytest.raises(ValueError):
        Cigar.coerce('3') # Not enough Ms


cigar_hit_ref_cut_cases = [
    # # Trivial cases
    (CigarHit('4M', r_st=1, r_ei=4, q_st=1, q_ei=4), 2.5,
     [CigarHit('2M', r_st=1, r_ei=2, q_st=1, q_ei=2),
      CigarHit('2M', r_st=3, r_ei=4, q_st=3, q_ei=4)]),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 3.5,
     [CigarHit('3M', r_st=1, r_ei=3, q_st=1, q_ei=3),
      CigarHit('6M', r_st=4, r_ei=9, q_st=4, q_ei=9)]),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 4.5,
     [CigarHit('4M', r_st=1, r_ei=4, q_st=1, q_ei=4),
      CigarHit('5M', r_st=5, r_ei=9, q_st=5, q_ei=9)]),

    (CigarHit('9M', r_st=0, r_ei=8, q_st=0, q_ei=8), 3.5,
     [CigarHit('4M', r_st=0, r_ei=3, q_st=0, q_ei=3),
      CigarHit('5M', r_st=4, r_ei=8, q_st=4, q_ei=8)]),

    # Simple cases
    (CigarHit('9M9D9M', r_st=1, r_ei=27, q_st=1, q_ei=18), 3.5,
     [CigarHit('3M', r_st=1, r_ei=3, q_st=1, q_ei=3),
      CigarHit('6M9D9M', r_st=4, r_ei=27, q_st=4, q_ei=18)]),

    (CigarHit('9M9D9M', r_st=1, r_ei=27, q_st=1, q_ei=18), 20.5,
     [CigarHit('9M9D2M', r_st=1, r_ei=20, q_st=1, q_ei=11),
      CigarHit('7M', r_st=21, r_ei=27, q_st=12, q_ei=18)]),

    (CigarHit('9M9I9M', r_st=1, r_ei=18, q_st=1, q_ei=27), 3.5,
     [CigarHit('3M', r_st=1, r_ei=3, q_st=1, q_ei=3),
      CigarHit('6M9I9M', r_st=4, r_ei=18, q_st=4, q_ei=27)]),

    (CigarHit('9M9I9M', r_st=1, r_ei=18, q_st=1, q_ei=27), 13.5 or 27/2,
     [CigarHit('9M9I4M', r_st=1, r_ei=13, q_st=1, q_ei=22),
      CigarHit('5M', r_st=14, r_ei=18, q_st=23, q_ei=27)]),

    # Ambigous cases
    (CigarHit('9M9D9M', r_st=1, r_ei=27, q_st=1, q_ei=18), 13.5 or 27/2,
     [CigarHit('9M4D', r_st=1, r_ei=13, q_st=1, q_ei=9),
      CigarHit('5D9M', r_st=14, r_ei=27, q_st=10, q_ei=18)]),

    (CigarHit('9M9I9M', r_st=1, r_ei=18, q_st=1, q_ei=27), 9.2,
     [CigarHit('9M1I', r_st=1, r_ei=9, q_st=1, q_ei=10),
      CigarHit('8I9M', r_st=10, r_ei=18, q_st=11, q_ei=27)]),

    # Edge cases
    (CigarHit('9M9I9M', r_st=1, r_ei=18, q_st=1, q_ei=27), 9.5, # no middlepoint
     [CigarHit('9M5I', r_st=1, r_ei=9, q_st=1, q_ei=14),
      CigarHit('4I9M', r_st=10, r_ei=18, q_st=15, q_ei=27)]),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 8.5, # one is singleton
     [CigarHit('8M', r_st=1, r_ei=8, q_st=1, q_ei=8),
      CigarHit('1M', r_st=9, r_ei=9, q_st=9, q_ei=9)]),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 9.5, # one is empty
     [CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9),
      CigarHit('', r_st=10, r_ei=9, q_st=10, q_ei=9)]),

    (CigarHit('7M', r_st=3, r_ei=9, q_st=3, q_ei=9), 2.5, # one is empty
     [CigarHit('', r_st=3, r_ei=2, q_st=3, q_ei=2),
      CigarHit('7M', r_st=3, r_ei=9, q_st=3, q_ei=9)]),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 0.5, # one is empty around 0
     [CigarHit('', r_st=1, r_ei=0, q_st=1, q_ei=0),
      CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9)]),

    (CigarHit('9M', r_st=0, r_ei=8, q_st=0, q_ei=8), -0.5, # another one is empty and negative
     [CigarHit('', r_st=0, r_ei=-1, q_st=0, q_ei=-1),
      CigarHit('9M', r_st=0, r_ei=8, q_st=0, q_ei=8)]),

    (CigarHit('9D', r_st=1, r_ei=9, q_st=1, q_ei=0), 3.5,
     [CigarHit('3D', r_st=1, r_ei=3, q_st=1, q_ei=0),
      CigarHit('6D', r_st=4, r_ei=9, q_st=1, q_ei=0)]),

    (CigarHit('9D', r_st=0, r_ei=8, q_st=0, q_ei=-1), -0.5,
     [CigarHit('', r_st=0, r_ei=-1, q_st=0, q_ei=-1),
      CigarHit('9D', r_st=0, r_ei=8, q_st=0, q_ei=-1)]),

    (CigarHit('2=1X2N1N2=1H2S', r_st=1, r_ei=8, q_st=1, q_ei=7), 3.5,
     [CigarHit('2=1X', r_st=1, r_ei=3, q_st=1, q_ei=3),
      CigarHit('3N2=1H2S', r_st=4, r_ei=8, q_st=4, q_ei=7)]),

    # Negative cases
    (CigarHit('9M9I9M', r_st=1, r_ei=18, q_st=1, q_ei=27), 20.5,
     IndexError("20.5 is bigger than reference (18)")),

    (CigarHit('', r_st=2, r_ei=1, q_st=2, q_ei=1), 2.5,
     IndexError("Empty string cannot be cut")),

    (CigarHit('', r_st=2, r_ei=1, q_st=2, q_ei=1), 1.5,
     IndexError("Empty string cannot be cut")),

    (CigarHit('9I', r_st=1, r_ei=0, q_st=1, q_ei=9), 3.5,
     IndexError("Out of reference bounds")),

    (CigarHit('9M', r_st=1, r_ei=9, q_st=1, q_ei=9), 4,
     ValueError("Cut point must not be an integer")),

]

@pytest.mark.parametrize('hit, cut_point, expected_result', cigar_hit_ref_cut_cases)
def test_cigar_hit_ref_cut(hit, cut_point, expected_result):
    if isinstance(expected_result, Exception):
        with pytest.raises(type(expected_result)):
            hit.cut_reference(cut_point)

    else:
        expected_left, expected_right = expected_result
        left, right = hit.cut_reference(cut_point)
        assert expected_left == left
        assert expected_right == right


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in cigar_hit_ref_cut_cases
                                            if not isinstance(x[2], Exception)])
def test_cigar_hit_ref_cut_add_prop(hit, cut_point):
    left, right = hit.cut_reference(cut_point)
    assert left + right == hit


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in cigar_hit_ref_cut_cases
                                            if not isinstance(x[2], Exception)])
def test_cigar_hit_ref_cut_add_prop_exhaustive(hit, cut_point):
    percentage = cut_point - floor(cut_point)

    for cut_point in range(hit.r_st, hit.r_ei + 2):
        left, right = hit.cut_reference(cut_point - percentage)
        assert left + right == hit


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in cigar_hit_ref_cut_cases
                                            if not isinstance(x[2], Exception)])
def test_cigar_hit_ref_cut_add_associativity(hit, cut_point):
    percentage = cut_point - floor(cut_point)

    for ax_cut in range(hit.r_st, hit.r_ei + 2):
        a, x = hit.cut_reference(ax_cut - percentage)

        for bc_cut in range(a.r_ei + 1, hit.r_ei + 2):
            if x.ref_length == 0: continue

            b, c = x.cut_reference(bc_cut - percentage)

            assert (a + b) + c == a + (b + c)


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases])
def test_cigar_hit_lstrip_is_stringlike(hit):
    all_chars = CIGAR_OP_MAPPING.keys()

    actions_of = lambda s: (x for x in s if x in all_chars)

    for r in range(len(all_chars) + 1):
        for char_set in itertools.combinations(all_chars, r):
            actions = set(map(parse_cigar_operation, char_set))
            chars = ''.join(char_set)

            p = lambda x: ''.join(actions_of(str(x.cigar)))
            g = lambda x: x.lstrip(actions)
            h = lambda x: x.lstrip(chars)

            assert p(g(hit)) == h(p(hit))


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases])
def test_cigar_hit_rstrip_is_stringlike(hit):
    all_chars = CIGAR_OP_MAPPING.keys()

    actions_of = lambda s: (x for x in s if x in all_chars)

    for r in range(len(all_chars) + 1):
        for char_set in itertools.combinations(all_chars, r):
            actions = set(map(parse_cigar_operation, char_set))
            chars = ''.join(char_set)

            p = lambda x: ''.join(actions_of(str(x.cigar)))
            g = lambda x: x.rstrip(actions)
            h = lambda x: x.rstrip(chars)

            assert p(g(hit)) == h(p(hit))


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases
                                 if not isinstance(x[2], Exception)])
def test_cigar_hit_gaps_no_m_or_i(hit):
    gaps = list(hit.gaps())

    if 'D' in str(hit.cigar):
        assert len(gaps) > 0

    for gap in gaps:
        assert 'M' not in str(gap.cigar)
        assert 'I' not in str(gap.cigar)


@pytest.mark.parametrize("reference_seq, query_seq, cigar, expected_reference, expected_query", [
    ('ACTG',   'ACTG',   '4M',     'ACTG',   'ACTG'),
    ('ACTG',   '',       '4D',     'ACTG',   '----'),
    ('',       'ACTG',   '4I',     '----',   'ACTG'),
    ('ACTGAC', 'ACAC',   '2M2D2M', 'ACTGAC', 'AC--AC'),
    ('ACAC',   'ACTGAC', '2M2I2M', 'AC--AC', 'ACTGAC'),
    ('GCTATGGGAA', 'GCTATGGGAA', '5M3D2M', 'GCTATGGGAA', 'GCTAT---GG'),
    ('ACTG',   'ACTG', '2M99H77P2M', 'ACTG', 'ACTG'), # Ignores non-consuming operations.
])
def test_cigar_to_msa(reference_seq, query_seq, cigar, expected_reference, expected_query):
    assert Cigar.coerce(cigar).to_msa(reference_seq, query_seq) \
        == (expected_reference, expected_query)


@pytest.mark.parametrize("cigar, reference_seq, query_seq", [
    ('10M', 'A' * 3, 'A' * 10), # reference is shorter than CIGAR
    ('10M', 'A' * 10, 'A' * 3), # query is shorter than CIGAR
    ('10D', 'A' * 3, 'A' * 3),
    ('10I', 'A' * 3, 'A' * 3),
])
def test_illigal_cigar_to_msa(cigar, reference_seq, query_seq):
    with pytest.raises(ValueError):
        Cigar.coerce(cigar).to_msa(reference_seq, query_seq)
