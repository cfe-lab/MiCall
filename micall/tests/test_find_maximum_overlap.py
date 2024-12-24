import pytest
from micall.utils.find_maximum_overlap import find_maximum_overlap, show_maximum_overlap, cli_main


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
        assert find_maximum_overlap(left, right) == expected
    else:
        with pytest.raises(expected):
            find_maximum_overlap(left, right)


print_cases = \
    [
        ('tttttxxxxx', 'uuuuuxxx', '''\
tttttxxxxx
uuuuuxxx--
'''),
        ('tttttxxx', 'uuuuuxxxxx', '''\
--tttttxxx
uuuuuxxxxx
'''),
        ('xxxttttt', 'xxxxxuuuuu', '''\
--xxxttttt
xxxxxuuuuu
'''),
        ('xxxxxttttt', 'xxxuuuuu', '''\
xxxxxttttt
xxxuuuuu--
'''),
        ('xxxxttttt', 'uuuuuxxx', '''\
-----xxxxttttt
uuuuuxxx------
'''),
        ('aaaaxxxxttttt', 'uxxxk', '''\
aaaaxxxxttttt
---uxxxk-----
'''
         ),
        ('tttttxxxx', 'xxxuuuuu', '''\
tttttxxxx-----
------xxxuuuuu
'''),
        ('tttttxxxx', 'xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu', '''\
tttttxxxx-------------------------------
------xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
'''),
        ('tttttxxxxpppppppppppppppppppp', 'xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu', '''\
tttttxxxxpppppppppppppppppppp-----------
------xxxuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
'''),
]


@pytest.mark.parametrize("left, right, expected", print_cases)
def test_print(left, right, expected):
    shift = find_maximum_overlap(left, right)
    ret = show_maximum_overlap(left, right, shift)
    assert ret == expected


@pytest.mark.parametrize("left, right, expected", print_cases)
def test_main(left, right, expected, capsys):
    cli_main([left, right])
    ret = capsys.readouterr().out
    assert ret == expected
