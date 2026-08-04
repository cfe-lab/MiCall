import pytest
import random

from micall.tests.utils import fixed_random_seed
from micall.utils.overlap_stitcher import calculate_concordance_norm, disambiguate_concordance, exp_dropoff_array


def test_concordance_same_length_inputs():
    with pytest.raises(ValueError):
        calculate_concordance_norm("abc", "ab")


def test_concordance_completely_different_strings():
    result = calculate_concordance_norm("a" * 30, "b" * 30)
    assert any(n > 0 for n in result)


def generate_random_string_pair(length):
    left = "".join(random.choice("ACGT") for _ in range(length))
    right = "".join(random.choice("ACGT") for _ in range(length))
    return left, right


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ("aaaaa", "aaaaa", [0.0, 0.78, 1.0, 0.78, 0.0]),
        ("abcdd", "abcdd", [0.0, 0.78, 1.0, 0.78, 0.0]),
        ("aaaaaaaa", "baaaaaab", [0.0, 0.95, 0.99, 1.0, 1.0, 0.99, 0.95, 0.0]),
        ("aaaaaaaa", "aaaaaaab", [0.94, 0.98, 0.99, 1.0, 0.99, 0.98, 0.94, 0.0]),
        ("aaaaaaaa", "aaaaabbb", [0.96, 0.99, 1.0, 0.99, 0.96, 0.02, 0.0, 0.02]),
        ("aaaaaaaaaaa", "aaaaabaaaaa", [0.96, 0.99, 1.0, 0.99, 0.96, 0.0, 0.96, 0.99, 1.0, 0.99, 0.96]),
        ("aaaaaaaaaaa", "aaaabbbaaaa", [0.98, 1.0, 1.0, 0.98, 0.02, 0.0, 0.02, 0.98, 1.0, 1.0, 0.98]),
        ("aaaaaaaaaaa", "aaabbbbbaaa", [0.98, 1.0, 0.98, 0.04, 0.01, 0.0, 0.01, 0.04, 0.98, 1.0, 0.98]),
        ("aaaaaaaa", "aaabbaaa", [0.98, 1.0, 0.98, 0.0, 0.0, 0.98, 1.0, 0.98]),
        ("aaaaa", "bbbbb", [1.0, 0.22, 0.0, 0.22, 1.0]),
        ("", "", []),
    ],
)
def test_concordance_simple(left, right, expected):
    result = [round(float(x), 2) for x in calculate_concordance_norm(left, right)]
    assert result == expected


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ("a" * 128, "a" * 128, 64),
        ("a" * 128, "a" * 64 + "b" * 64, 32),
        ("a" * 128, "a" * 64 + "ba" * 32, 32),
        ("a" * 128, "a" * 54 + "b" * 20 + "a" * 54, 28),  # two peaks
        ("a" * 128, "a" * 63 + "b" * 2 + "a" * 63, 32),  # two peaks
        ("a" * 1280, "b" * 640 + "a" * 640, round(1280 * 3 / 4)),
        ("a" * 128, "b" * 48 + "a" * 32 + "b" * 48, 64),
        (
            "a" * 128,
            "b" * 48 + "a" * 15 + "ab" + "a" * 15 + "b" * 48,
            48 + 15 // 2,
        ),  # two peaks - choosing 1st
        (
            "a" * 128,
            "b" * 48 + "a" * 15 + "ba" + "a" * 15 + "b" * 48,
            48 + 15 + 16 // 2,
        ),  # two peaks - choosing 2nd
        (
            "a" * 128,
            "b" * 48 + "a" * 15 + "bb" + "a" * 15 + "b" * 48,
            48 + 15 // 2,
        ),  # two peaks - choosing 1st
    ],
)
def test_concordance_simple_index(left, right, expected):
    concordance = calculate_concordance_norm(left, right)
    concordance_d = list(disambiguate_concordance(concordance))
    index = max(range(len(concordance)), key=lambda i: concordance_d[i])
    if abs(index - expected) > 1:
        assert index == expected


def generate_test_cases(num_cases):
    with fixed_random_seed(42):
        length = random.randint(1, 80)
        return [generate_random_string_pair(length) for _ in range(num_cases)]


concordance_cases = generate_test_cases(num_cases=100)


@pytest.mark.parametrize("left, right", concordance_cases)
def test_concordance_output_range(left, right):
    result = calculate_concordance_norm(left, right)
    assert all(
        0 <= n <= 1 for n in result
    ), "All values in result should be between 0 and 1"


@pytest.mark.parametrize("left, right", concordance_cases)
def test_concordance_higher_if_more_matches_added(left, right):
    # Insert exact matches in the middle
    matching_sequence = "A" * 30
    insert_position = len(left) // 2
    new_left = (
        left[:insert_position]
        + matching_sequence
        + left[insert_position + len(matching_sequence):]
    )
    new_right = (
        right[:insert_position]
        + matching_sequence
        + right[insert_position + len(matching_sequence):]
    )

    old_conc = calculate_concordance_norm(left, right)
    new_conc = calculate_concordance_norm(new_left, new_right)
    old_average = sum(old_conc) / len(old_conc)
    new_average = sum(new_conc) / len(new_conc)
    assert old_average <= new_average


@pytest.mark.parametrize("left, right", concordance_cases)
def test_concordance_higher_in_matching_areas(left, right):
    # Insert exact matches in the middle
    matching_sequence = "A" * 30
    insert_position = len(left) // 2
    new_left = (
        left[:insert_position]
        + matching_sequence
        + left[insert_position + len(matching_sequence):]
    )
    new_right = (
        right[:insert_position]
        + matching_sequence
        + right[insert_position + len(matching_sequence):]
    )

    concordance_scores = calculate_concordance_norm(new_left, new_right)

    # Check concordance in the matching area
    matching_area_concordance = concordance_scores[
        insert_position:insert_position + len(matching_sequence)
    ]

    # Calculate average concordance inside and outside the matching area
    average_inside = sum(matching_area_concordance) / len(matching_sequence)
    average_outside = (sum(concordance_scores) - sum(matching_area_concordance)) / (
        len(concordance_scores) - len(matching_sequence)
    )

    # Assert that the concordance is indeed higher in the matching area
    assert (
        average_inside > average_outside
    ), "Concordance in matching areas should be higher than in non-matching areas"


@pytest.mark.parametrize(
    "input, output",
    [
        ([0, 0, 0, 1, 0, 0, 0], [0.125, 0.25, 0.5, 1, 0.5, 0.25, 0.125]),
        ([0, 0, 1, 1, 0, 0, 0], [0.25, 0.5, 1, 1, 0.5, 0.25, 0.125]),
        ([0, 0, 1, 1, 0, 1, 0], [0.25, 0.5, 1, 1, 0.5, 1, 0.5]),
        ([0, 0, 1, 0, 0, 0, 0], [0.25, 0.5, 1, 0.5, 0.25, 0.125, 0.0625]),
        ([0, 0, 0, 8, 0, 0, 0], [1, 2, 4, 8, 4, 2, 1]),
        ([0, 0, 0, -8, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]),
        ([0, 0, 1, -8, 0, 0, 0], [0.25, 0.5, 1, 0.5, 0.25, 0.125, 0.0625]),
        ([0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]),
        ([1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]),
        ([], []),
    ],
)
def test_exp_dropoff_simple(input, output):
    exp_dropoff_array(input)
    assert input == output
