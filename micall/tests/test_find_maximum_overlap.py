import pytest
from micall.utils.find_maximum_overlap import find_maximum_overlap, show_maximum_overlap, cli_main


@pytest.mark.parametrize(
    "left, right, expected",
    [
        ('aaaaaxxxx', 'xxxbbbbb', -3),
        ('aaaaaxxxx', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', -3),
        ('aaaaaxxxxcccccccccccccccccccc', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', -23),
        ('aaaaaxxxxx', 'bbbbbxxx', -10),
        ('bbbbbxxx', 'aaaaaxxxxx', -10),
        ('aaaaaxxx', 'bbbbbxxxxx', -10),
        ('xxxaaaaa', 'xxxxxbbbbb', -10),
        ('xxxxxaaaaa', 'xxxbbbbb', -10),
        ('xxxxaaaaa', 'bbbbbxxx', -14),
        ('aaaaxxxxaaaaa', 'bxxxk', -10),
        ('aaaa', 'bbbb', 0),
        ('aaaax', 'xbbbb', -1),
        ('', 'aaaaaxxxx', ValueError),
        ('aaaaaxxxx', '', ValueError),
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
        ('bbbbbxxx', 'aaaaaxxxxx', '''\
--bbbbbxxx
aaaaaxxxxx
'''),
        ('aaaaaxxxxx', 'bbbbbxxx', '''\
aaaaaxxxxx
bbbbbxxx--
'''),
        ('aaaaaxxx', 'bbbbbxxxxx', '''\
--aaaaaxxx
bbbbbxxxxx
'''),
        ('xxxaaaaa', 'xxxxxbbbbb', '''\
--xxxaaaaa
xxxxxbbbbb
'''),
        ('xxxxxaaaaa', 'xxxbbbbb', '''\
xxxxxaaaaa
xxxbbbbb--
'''),
        ('xxxxaaaaa', 'bbbbbxxx', '''\
-----xxxxaaaaa
bbbbbxxx------
'''),
        ('aaaaxxxxaaaaa', 'bxxxk', '''\
aaaaxxxxaaaaa
---bxxxk-----
'''
         ),
        ('aaaaaxxxx', 'xxxbbbbb', '''\
aaaaaxxxx-----
------xxxbbbbb
'''),
        ('aaaaaxxxx', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', '''\
aaaaaxxxx-------------------------------
------xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
'''),
        ('aaaaaxxxxcccccccccccccccccccc', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb', '''\
aaaaaxxxxcccccccccccccccccccc-----------
------xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
'''),
        ('aaaa', 'bbbb', '''\
aaaa----
----bbbb
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
