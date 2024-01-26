import pytest
from typing import List, Tuple
from math import floor
import itertools
import re

from micall.utils.consensus_aligner import CigarActions
from micall.utils.cigar_tools import Cigar, CigarHit, connect_cigar_hits, CoordinateMapping


cigar_mapping_cases = [
    # Simple cases
    ('3M',     {0: 0, 1: 1, 2: 2},  # exact mapping
               {0: 0, 1: 1, 2: 2}), # closest mapping
    ('1M1D1M', {0: 0, 2: 1},        # exact mapping
               {0: 0, 1: 0, 2: 1}), # closest mapping
    ('1M1I1M', {0: 0, 1: 2},
               {0: 0, 1: 2}),
    ('2M2D2M', {0: 0, 1: 1, 4: 2, 5: 3},
               {0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 3}),
    ('2M2I2M', {0: 0, 1: 1, 2: 4, 3: 5},
               {0: 0, 1: 1, 2: 4, 3: 5}),
    ('3M1D3M', {0: 0, 1: 1, 2: 2, 4: 3, 5: 4, 6: 5},
               {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4, 6: 5}),
    ('3M1I3M', {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 6},
               {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 6}),
    ('7M1I3M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 9, 9: 10},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 9, 9: 10}),
    ('5M2D4M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 7: 5, 8: 6, 9: 7, 10: 8},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 4, 6: 4, 7: 5, 8: 6, 9: 7, 10: 8}),
    ('5M3I4M', {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 8, 6: 9, 7: 10, 8: 11},
               {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 8, 6: 9, 7: 10, 8: 11}),
    ('1M1D',   {0: 0},
               {0: 0, 1: 0}),
    ('1M1I',   {0: 0},
               {0: 0}),
    ('1I1M',   {0: 1},
               {0: 1}),
    ('1D1M',   {1: 0},
               {1: 0, 0: None}),

    # Multiple deletions and insertions
    ('2M2D2M2I2M', {0: 0, 1: 1, 4: 2, 5: 3, 6: 6, 7: 7},
                   {0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 7}),
    ('2M2I2M2D2M', {0: 0, 1: 1, 2: 4, 3: 5, 6: 6, 7: 7},
                   {0: 0, 1: 1, 2: 4, 3: 5, 4: 5, 5: 5, 6: 6, 7: 7}),
    ('2=1X2N1N2=1H2S', {0: 0, 1: 1, 2: 2, 6: 3, 7: 4},
                   {0: 0, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 3, 7: 4}),
    ('2M2D2M2I2M', {0: 0, 1: 1, 4: 2, 5: 3, 6: 6, 7: 7},
                   {0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 7}),
    ('3=1X2N1N2=1H2S', {0: 0, 1: 1, 2: 2, 3: 3, 7: 4, 8: 5},
                   {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 3, 6: 3, 7: 4, 8: 5}),

    # Edge cases
    ('', {}, {}),
    ('3I3D',   {},
               {0: None, 1: None, 2: None}),
    ('3D3I',   {},
               {0: None, 1: None, 2: None}),
    ('12I', {}, {}),
    ('12D', {}, {k: None for k in range(12)}),
]


@pytest.mark.parametrize("cigar_str, expected_mapping", [(x[0], x[1]) for x in cigar_mapping_cases])
def test_cigar_to_coordinate_mapping(cigar_str, expected_mapping):
    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    assert expected_mapping == mapping.ref_to_query
    assert expected_mapping == {i: mapping.ref_to_query[i]
                                for i in mapping.ref_to_query.keys()}


@pytest.mark.parametrize("cigar_str", [x[0] for x in cigar_mapping_cases])
def test_cigar_to_coordinate_bijection_property(cigar_str):
    inverse = lambda d: {v: k for k, v in d.items()}

    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    assert mapping.query_to_ref == inverse(mapping.ref_to_query)
    assert mapping.ref_to_query == inverse(mapping.query_to_ref)
    assert mapping.ref_to_query == inverse(inverse(mapping.ref_to_query))
    assert mapping.query_to_ref == inverse(inverse(mapping.query_to_ref))


@pytest.mark.parametrize("cigar_str, expected_leftmax_mapping", [(x[0], x[2]) for x in cigar_mapping_cases])
def test_cigar_to_coordinate_mapping_leftmax(cigar_str, expected_leftmax_mapping):
    mapping = Cigar.coerce(cigar_str).coordinate_mapping

    def test():
        fullrange = {i: mapping.ref_to_query.left_max(i)
                     for i in mapping.ref_to_query.domain}
        assert expected_leftmax_mapping == fullrange

    if isinstance(expected_leftmax_mapping, Exception):
        with pytest.raises(type(expected_leftmax_mapping)):
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
    assert mapping.ref_to_query.get(0, None) == None
    assert mapping.query_to_ref.get(0, None) == None
    assert expected_mapping \
        == {i: mapping.ref_to_query[i]
            for i in mapping.ref_to_query.keys()}


@pytest.mark.parametrize("cigar_str, expected_leftmax_mapping", [(x[0], x[2]) for x in cigar_mapping_cases])
def test_cigar_hit_to_coordinate_mapping_leftmax(cigar_str, expected_leftmax_mapping):
    cigar = Cigar.coerce(cigar_str)
    hit = CigarHit(cigar, r_st=5, r_ei=(5 + cigar.ref_length - 1), q_st=7, q_ei=(7 + cigar.query_length - 1))
    mapping = hit.coordinate_mapping

    def test(expected):
        # Coordinates are translated by q_st and r_st.
        fullrange = {i: mapping.ref_to_query.left_max(i)
                     for i in mapping.ref_to_query.domain}
        assert expected == fullrange

    if isinstance(expected_leftmax_mapping, Exception):
        with pytest.raises(type(expected_leftmax_mapping)):
            test(expected_leftmax_mapping)
    else:
        test({k + hit.r_st: v + hit.q_st if v is not None else v for (k, v) in expected_leftmax_mapping.items()})


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


CIGAR_REGEX = re.compile(r"(.*)@([0-9]+)->([0-9]+)")
def parsed_hit(string):
    match = CIGAR_REGEX.match(string)
    assert match, f"Cannot parse {string}"
    cigar_str, q_st, r_st = match.groups()
    cigar = Cigar.coerce(cigar_str)
    r_ei = int(r_st) + cigar.ref_length - 1
    q_ei = int(q_st) + cigar.query_length - 1
    return CigarHit(cigar, int(r_st), int(r_ei), int(q_st), int(q_ei))


cigar_hit_ref_cut_cases = [
    # Trivial cases
    ('4M@1->1', 2.5,
     ['2M@1->1', '2M@3->3']),

    ('9M@1->1', 3.5,
     ['3M@1->1', '6M@4->4']),

    ('9M@1->1', 4.5,
     ['4M@1->1', '5M@5->5']),

    ('9M@0->0', 3.5,
     ['4M@0->0', '5M@4->4']),

    # Simple cases
    ('9M9D9M@1->1', 3.5,
     ['3M@1->1', '6M9D9M@4->4']),

    ('9M9D9M@1->1', 20.5,
     ['9M9D2M@1->1', '7M@12->21']),

    ('9M9I9M@1->1', 3.5,
     ['3M@1->1', '6M9I9M@4->4']),

    ('9M9I9M@1->1', 13.5 or 27/2,
     ['9M9I4M@1->1', '5M@23->14']),

    ('5M6I@1->1', 3.5,
     ['3M@1->1', '2M6I@4->4']),

    ('6I5M@1->1', 3.5,
     ['6I3M@1->1', '2M@10->4']),

    ('5M6D@1->1', 3.5,
     ['3M@1->1', '2M6D@4->4']),

    ('6D5M@1->1', 3.5,
     ['3D@1->1', '3D5M@1->4']),

    ('5M6D@1->1', 7.5,
     ['5M2D@1->1', '4D@6->8']),

    ('6D5M@1->1', 7.5,
     ['6D1M@1->1', '4M@2->8']),

    ('6D5M@1->1', 6.5,
     ['6D@1->1', '5M@1->7']),

    # Ambigous cases
    ('9M9D9M@1->1', 13.5 or 27/2,
     ['9M4D@1->1', '5D9M@10->14']),

    ('9M9I9M@1->1', 9.2,
     ['9M1I@1->1', '8I9M@11->10']),

    ('9M9D9I9M@1->1', 13.5 or 27/2,
     ['9M4D@1->1', '5D9I9M@10->14']),

    ('9M9I9D9M@1->1', 13.5 or 27/2,
     ['9M9I4D@1->1', '5D9M@19->14']),

    ('9M9D9I9D9I9D9M@1->1', 13.5 or 27/2,
     ['9M4D@1->1', '5D9I9D9I9D9M@10->14']),

    ('9M9I9D9I9D9I9M@1->1', 13.5 or 27/2,
     ['9M9I4D@1->1', '5D9I9D9I9M@19->14']),

    ('1M1I1D1M@1->1', 1.5, # same as previous 2 cases but smaller
     ['1M1I@1->1', '1D1M@3->2']),

    ('1M1D1I1M@1->1', 1.5, # same as previous 2 cases but smaller
     ['1M@1->1', '1D1I1M@2->2']),

    # Edge cases
    ('9M9I9M@1->1', 9.5, # no middlepoint
     ['9M5I@1->1', '4I9M@15->10']),

    ('9M@1->1', 8.5,
     ['8M@1->1', '1M@9->9']),

    ('9M@1->1', 9.5,
     ['9M@1->1', '@10->10']),

    ('7M@3->3', 2.5,
     ['@3->3', '7M@3->3']),

    ('9M@1->1', 0.5,
     ['@1->1', '9M@1->1']),

    ('9M@0->0', -0.5,
     ['@0->0', '9M@0->0']),

    ('9D@1->1', 3.5,
     ['3D@1->1', '6D@1->4']),

    ('9D@0->0', -0.5,
     ['@0->0', '9D@0->0']),

    ('1M7I1M@1->1', 1.5,
     ['1M4I@1->1', '3I1M@6->2']),

    ('1M6I1M@1->1', 1.5,
     ['1M3I@1->1', '3I1M@5->2']),

    ('1M7I1M@1->1', 1.999,
     ['1M7I@1->1', '1M@9->2']),

    ('1M7I1M@1->1', 1.001,
     ['1M@1->1', '7I1M@2->2']),

    ('2=1X2N1N2=1H2S@1->1', 3.5,
     ['2=1X@1->1', '3N2=1H2S@4->4']),

    # Negative cases
    ('9M9I9M@1->1', 20.5,
     IndexError("20.5 is bigger than reference (18)")),

    ('@2->2', 2.5,
     IndexError("Empty string cannot be cut")),

    ('@2->2', 1.5,
     IndexError("Empty string cannot be cut")),

    ('9I@1->1', 3.5,
     IndexError("Out of reference bounds")),

    ('9M@1->1', 4,
     ValueError("Cut point must not be an integer")),

]


@pytest.mark.parametrize('hit, cut_point, expected_result', cigar_hit_ref_cut_cases)
def test_cigar_hit_ref_cut(hit, cut_point, expected_result):
    hit = parsed_hit(hit)

    if isinstance(expected_result, Exception):
        with pytest.raises(type(expected_result)):
            hit.cut_reference(cut_point)

    else:
        expected_result = list(map(parsed_hit, expected_result))
        expected_left, expected_right = expected_result
        left, right = hit.cut_reference(cut_point)
        assert expected_left == left
        assert expected_right == right


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in cigar_hit_ref_cut_cases
                                            if not isinstance(x[2], Exception)])
def test_cigar_hit_ref_cut_add_prop(hit, cut_point):
    hit = parsed_hit(hit)
    left, right = hit.cut_reference(cut_point)
    assert left + right == hit


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in
                                            [x for x in cigar_hit_ref_cut_cases
                                             if not isinstance(x[2], Exception)]])
def test_cigar_hit_ref_cut_add_prop_exhaustive(hit, cut_point):
    hit = parsed_hit(hit)
    percentage = cut_point - floor(cut_point)

    for cut_point in range(hit.r_st, hit.r_ei + 2):
        left, right = hit.cut_reference(cut_point - percentage)
        assert left + right == hit


lstrip_reference_cases = [
    ('9M@1->1', '9M@1->1'),
    ('5M6D@1->1', '5M6D@1->1'),
    ('6D5M@1->1', '6D5M@1->1'),
    ('6I5M@1->1', '5M@7->1'),
    ('6I4D5M@1->1', '4D5M@7->1'),
    ('6D4I5M@1->1', '6D5M@5->1'),
    ('3D3D4I5M@1->1', '6D5M@5->1'),
    ('3I3I4D5M@1->1', '4D5M@7->1'),
    ('3D2I3D2I5M@1->1', '6D5M@5->1'),
    ('3I2D3I2D5M@1->1', '4D5M@7->1'),
    ('4D6I5M@1->1', '4D5M@7->1'),
    ('4I6D5M@1->1', '6D5M@5->1'),
    ('6I4D@1->1', '4D@7->1'),
    ('6D4I@1->1', '6D4I@1->1'),
    ('4D6I@1->1', '4D6I@1->1'),
    ('4I6D@1->1', '6D@5->1'),
    ('@1->1', '@1->1'),
]


rstrip_reference_cases = [
    ('9M@1->1', '9M@1->1'),
    ('5M6D@1->1', '5M6D@1->1'),
    ('5M6I@1->1', '5M@1->1'),
    ('6D5M@1->1', '6D5M@1->1'),
    ('5M4I6D@1->1', '5M6D@1->1'),
    ('5M4D6I@1->1', '5M4D@1->1'),
    ('5M4I3D3D@1->1', '5M6D@1->1'),
    ('5M4D3I3I@1->1', '5M4D@1->1'),
    ('5M2I3D2I3D@1->1', '5M6D@1->1'),
    ('5M2D3I2D3I@1->1', '5M4D@1->1'),
    ('5M6D4I@1->1', '5M6D@1->1'),
    ('5M6I4D@1->1', '5M4D@1->1'),
    ('6D4I@1->1', '6D@1->1'),
    ('6I4D@1->1', '6I4D@1->1'),
    ('4I6D@1->1', '4I6D@1->1'),
    ('4D6I@1->1', '4D@1->1'),
    ('@1->1', '@1->1'),
]


lstrip_query_cases = [
    ('9M@1->1', '9M@1->1'),
    ('5M6D@1->1', '5M6D@1->1'),
    ('6D5M@1->1', '5M@1->7'),
    ('6I5M@1->1', '6I5M@1->1'),
    ('6I4D5M@1->1', '6I5M@1->5'),
    ('6D4I5M@1->1', '4I5M@1->7'),
    ('6D4I5M@1->1', '4I5M@1->7'),
    ('6I4D5M@1->1', '6I5M@1->5'),
    ('3D2I3D2I5M@1->1', '4I5M@1->7'),
    ('3I2D3I2D5M@1->1', '6I5M@1->5'),
    ('4D6I5M@1->1', '6I5M@1->5'),
    ('4I6D5M@1->1', '4I5M@1->7'),
    ('6I4D@1->1', '6I4D@1->1'),
    ('6D4I@1->1', '4I@1->7'),
    ('4D6I@1->1', '6I@1->5'),
    ('4I6D@1->1', '4I6D@1->1'),
    ('@1->1', '@1->1'),
]


rstrip_query_cases = [
    ('9M@1->1', '9M@1->1'),
    ('5M6D@1->1', '5M@1->1'),
    ('5M6I@1->1', '5M6I@1->1'),
    ('6D5M@1->1', '6D5M@1->1'),
    ('5M4I6D@1->1', '5M4I@1->1'),
    ('5M4D6I@1->1', '5M6I@1->1'),
    ('5M4I6D@1->1', '5M4I@1->1'),
    ('5M4D6I@1->1', '5M6I@1->1'),
    ('5M2I3D2I3D@1->1', '5M4I@1->1'),
    ('5M2D3I2D3I@1->1', '5M6I@1->1'),
    ('5M6D4I@1->1', '5M4I@1->1'),
    ('5M6I4D@1->1', '5M6I@1->1'),
    ('6D4I@1->1', '6D4I@1->1'),
    ('6I4D@1->1', '6I@1->1'),
    ('4I6D@1->1', '4I@1->1'),
    ('4D6I@1->1', '4D6I@1->1'),
    ('@1->1', '@1->1'),
]


strip_prop_cases_all = \
    [x[0] for x in cigar_hit_ref_cut_cases] + \
    [x[0] for x in lstrip_reference_cases] + \
    [x[0] for x in rstrip_reference_cases]


@pytest.mark.parametrize('hit, expected', lstrip_reference_cases)
def test_cigar_hit_lstrip_reference(hit, expected):
    hit = parsed_hit(hit)
    expected = parsed_hit(expected)
    assert expected == hit.lstrip_reference()


@pytest.mark.parametrize('hit, expected', rstrip_reference_cases)
def test_cigar_hit_rstrip_reference(hit, expected):
    hit = parsed_hit(hit)
    expected = parsed_hit(expected)
    assert expected == hit.rstrip_reference()



@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_reference_strip_combines_with_connect(hit):
    hit = parsed_hit(hit)

    for cut_point in range(hit.r_st - 1, hit.r_ei):
        left, right = hit.cut_reference(cut_point + hit.epsilon)

        left = left.rstrip_reference()
        right = right.lstrip_reference()

        assert left.connect(right).coordinate_mapping.ref_to_query \
            == hit.coordinate_mapping.ref_to_query


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_reference_strip_combines_with_add(hit):
    hit = parsed_hit(hit)

    for cut_point in range(hit.r_st - 1, hit.r_ei):
        left, right = hit.cut_reference(cut_point + hit.epsilon)

        left = left.rstrip_reference()
        right = right.lstrip_reference()

        if left.touches(right):
            assert left + right == hit


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_reference_strip_never_crashes(hit):
    hit = parsed_hit(hit)

    hit.rstrip_reference().lstrip_reference()
    hit.lstrip_reference().rstrip_reference()
    hit.lstrip_reference().lstrip_reference()
    hit.rstrip_reference().rstrip_reference()


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_reference_strip_is_idempotent(hit):
    hit = parsed_hit(hit)

    h1 = hit.rstrip_reference()
    assert h1 == h1.rstrip_reference() == h1.rstrip_reference().rstrip_reference()

    h1 = hit.lstrip_reference()
    assert h1 == h1.lstrip_reference() == h1.lstrip_reference().lstrip_reference()

    h1 = hit.lstrip_reference().rstrip_reference()
    assert h1 == h1.lstrip_reference() == h1.rstrip_reference()

    h1 = hit.rstrip_reference().lstrip_reference()
    assert h1 == h1.rstrip_reference() == h1.lstrip_reference()


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_reference_strips_are_commutative(hit):
    hit = parsed_hit(hit)

    assert hit.rstrip_reference().lstrip_reference() \
        == hit.lstrip_reference().rstrip_reference()


@pytest.mark.parametrize('hit, expected', lstrip_query_cases)
def test_cigar_hit_lstrip_query(hit, expected):
    hit = parsed_hit(hit)
    expected = parsed_hit(expected)
    assert expected == hit.lstrip_query()


@pytest.mark.parametrize('hit, expected', rstrip_query_cases)
def test_cigar_hit_rstrip_query(hit, expected):
    hit = parsed_hit(hit)
    expected = parsed_hit(expected)
    assert expected == hit.rstrip_query()


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_query_strip_combines_with_connect(hit):
    hit = parsed_hit(hit)

    for cut_point in range(hit.r_st - 1, hit.r_ei):
        left, right = hit.cut_reference(cut_point + hit.epsilon)

        left = left.rstrip_query()
        right = right.lstrip_query()

        assert left.connect(right).coordinate_mapping.ref_to_query \
            == hit.coordinate_mapping.ref_to_query


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_query_strip_combines_with_add(hit):
    hit = parsed_hit(hit)

    for cut_point in range(hit.r_st - 1, hit.r_ei):
        left, right = hit.cut_reference(cut_point + hit.epsilon)

        left = left.rstrip_query()
        right = right.lstrip_query()

        if left.touches(right):
            assert left + right == hit


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_query_strip_never_crashes(hit):
    hit = parsed_hit(hit)

    hit.rstrip_query().lstrip_query()
    hit.lstrip_query().rstrip_query()
    hit.lstrip_query().lstrip_query()
    hit.rstrip_query().rstrip_query()


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_query_strip_is_idempotent(hit):
    hit = parsed_hit(hit)

    h1 = hit.rstrip_query()
    assert h1 == h1.rstrip_query() == h1.rstrip_query().rstrip_query()

    h1 = hit.lstrip_query()
    assert h1 == h1.lstrip_query() == h1.lstrip_query().lstrip_query()

    h1 = hit.lstrip_query().rstrip_query()
    assert h1 == h1.lstrip_query() == h1.rstrip_query()

    h1 = hit.rstrip_query().lstrip_query()
    assert h1 == h1.rstrip_query() == h1.lstrip_query()


@pytest.mark.parametrize('hit', strip_prop_cases_all)
def test_cigar_hit_query_strips_are_commutative(hit):
    hit = parsed_hit(hit)

    assert hit.rstrip_query().lstrip_query() \
        == hit.lstrip_query().rstrip_query()


@pytest.mark.parametrize('hit, cut_point', [(x[0], x[1]) for x in cigar_hit_ref_cut_cases
                                            if not isinstance(x[2], Exception)])
def test_cigar_hit_ref_cut_add_associativity(hit, cut_point):
    hit = parsed_hit(hit)
    percentage = cut_point - floor(cut_point)

    for ax_cut in range(hit.r_st, hit.r_ei + 2):
        a, x = hit.cut_reference(ax_cut - percentage)

        for bc_cut in range(a.r_ei + 1, hit.r_ei + 2):
            if x.ref_length == 0: continue

            b, c = x.cut_reference(bc_cut - percentage)

            assert (a + b) + c == a + (b + c)


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases
                                 if not isinstance(x[2], Exception)])
def test_cigar_hit_deletions_no_m_or_i(hit):
    hit = parsed_hit(hit)
    deletions = list(hit.deletions())

    assert len(deletions) == len([op for op in str(hit.cigar) if op in 'DN'])

    for deletion in deletions:
        assert 'M' not in str(deletion.cigar)
        assert 'I' not in str(deletion.cigar)


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases
                                 if not isinstance(x[2], Exception)])
def test_cigar_hit_deletions_lengths(hit):
    hit = parsed_hit(hit)
    deletions = list(hit.deletions())

    for deletion in deletions:
        assert deletion.query_length == 0
        assert deletion.ref_length > 0
        assert deletion.coordinate_mapping.ref_to_query == {}


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases
                                 if not isinstance(x[2], Exception)])
def test_cigar_hit_insertions_no_m_or_i(hit):
    hit = parsed_hit(hit)
    insertions = list(hit.insertions())

    if 'I' in str(hit.cigar):
        assert len(insertions) > 0

    for insertion in insertions:
        assert 'M' not in str(insertion.cigar)
        assert 'D' not in str(insertion.cigar)


@pytest.mark.parametrize('hit', [x[0] for x in cigar_hit_ref_cut_cases
                                 if not isinstance(x[2], Exception)])
def test_cigar_hit_insertions_lengths(hit):
    hit = parsed_hit(hit)
    insertions = list(hit.insertions())

    for insertion in insertions:
        assert insertion.ref_length == 0
        assert insertion.query_length > 0
        assert insertion.coordinate_mapping.ref_to_query == {}


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


connect_cigar_hits_cases = [
    # Non-overlapping hits should be connected with deletions/insertions
    (
        ['4M@1->1', '4M@8->10'],
        ['4M5D3I4M@1->1']
    ),
    # Overlapping hits should ignore later ones
    (
        ['4M@1->1', '5M@3->3'],
        ['4M@1->1']
    ),
    # Touching hits should be simply concatenated
    (
        ['4M@1->1', '4M@5->5'],
        ['8M@1->1']
    ),
    # Hits that touch at only one boundary should combine just fine
    (
        ['3M@1->1', '6M@6->4'],
        ['3M2I6M@1->1']
    ),
    # Hits that are subsets of earlier hits should be ignored
    (
        ['8M@1->1', '3M@3->3'],
        ['8M@1->1']
    ),
    # Hits that are out of order should be connected if no overlap
    (
        ['3M@6->10', '3M@1->1'],
        ['3M6D2I3M@1->1']
    ),
    # Hits that overlap by a single base should prioritize the first hit and not combine
    (
        ['3M@1->1', '3M@3->3'],
        ['3M@1->1']
    ),
    # Non-overlapping hits in the query space but overlapping in reference space
    (
        ['5M@1->1', '1M@10->3'],
        ['5M@1->1']
    ),
    # Combining more than two hits
    (
        ['3M@1->1', '3M@7->7', '3M@16->12'],
        ['3M3D3I3M2D6I3M@1->1']
    ),
    # Combining hits including hard-clipping, which should be ignored in alignments
    (
        ['2H5M1H@3->1', '2H5M1H@13->11'],
        ['2H5M1H5D5I2H5M1H@3->1']
    ),
    # An empty list of hits should raise a ValueError
    (
        [],
        ValueError("Expected a non-empty list of cigar hits")
    ),
    # Before by reference, after by query
    (
        ['4M@8->1', '4M@1->10'],
        ['4M@8->1', '4M@1->10']
    ),
]
@pytest.mark.parametrize('hits, expected_result', connect_cigar_hits_cases)
def test_connect_cigar_hits(hits, expected_result):
    hits = list(map(parsed_hit, hits))

    if isinstance(expected_result, Exception):
        with pytest.raises(type(expected_result)):
            connect_cigar_hits(hits)
    else:
        expected_result = list(map(parsed_hit, expected_result))
        result = connect_cigar_hits(hits)
        assert expected_result == result
