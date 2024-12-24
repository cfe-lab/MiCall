import pytest
import micall.utils.merge_contigs as merge_contigs


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ('tttttxxxxx', 'uuuuuxxx', -10),
        ('tttttxxx', 'uuuuuxxxxx', -10),
        ('xxxttttt', 'xxxxxuuuuu', -10),
        ('xxxxxttttt', 'xxxuuuuu', -10),
        ('xxxxttttt', 'uuuuuxxx', -14),
        ('aaaaxxxxttttt', 'uxxxk', -10),
        ('tttttxxxx', 'xxxuuuuu', -3),
        ('tttttxxxx', 'xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu', -3),
        ('tttttxxxxpppppppppppppppppppp', 'xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu', -23),
        ('', 'tttttxxxx', ValueError),
        ('tttttxxxx', '', ValueError),
        ('', '', ValueError),
    ],
)
def test_maximum_overlap_cases(left, right, expected):
    if isinstance(expected, int):
        assert merge_contigs.find_maximum_overlap(left, right) == expected
    else:
        with pytest.raises(expected):
            merge_contigs.find_maximum_overlap(left, right)
