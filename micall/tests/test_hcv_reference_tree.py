from io import StringIO

from micall.utils.hcv_reference_tree import filter_hcv_fasta, check_tree


def test_add_subtype():
    raw_hcv = StringIO("""\
>Ref.1a.Foo
ACTAGGA
GAGATTT
>Ref.2b.Bar
TAGACT
""")
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.1a.Foo-1a
ACTAGGA
GAGATTT
>Ref.2b.Bar-2b
TAGACT
"""

    filter_hcv_fasta(raw_hcv, filtered_hcv)

    assert expected_filtered_hcv == filtered_hcv.getvalue()


def test_subtype_uppercase():
    raw_hcv = StringIO("""\
>Ref.1a.Foo
ACTAGGA
GAGATTT
>Ref.2B.Bar
TAGACT
""")
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.1a.Foo-1a
ACTAGGA
GAGATTT
>Ref.2B.Bar-2b
TAGACT
"""

    filter_hcv_fasta(raw_hcv, filtered_hcv)

    assert expected_filtered_hcv == filtered_hcv.getvalue()


def test_consensus():
    """ Consensus sequences have a different name format.
     
    They also need to have the parentheses replaced by brackets so they don't
    get dropped from the tree file.
    """
    raw_hcv = StringIO("""\
>Ref.CON_1a(142)
ACTAGGA
GAGATTT
>Ref.2b.Bar
TAGACT
""")
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.CON_1a[142]-1a
ACTAGGA
GAGATTT
>Ref.2b.Bar-2b
TAGACT
"""

    filter_hcv_fasta(raw_hcv, filtered_hcv)

    assert expected_filtered_hcv == filtered_hcv.getvalue()


def test_filter_missing_subtype():
    raw_hcv = StringIO("""\
>Ref.1a.Foo
ACTAGGA
GAGATTT
>Ref.2.Bar
TAGACT
""")
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.1a.Foo-1a
ACTAGGA
GAGATTT
"""
    expected_invalid_headers = ['>Ref.2.Bar']

    invalid_headers = filter_hcv_fasta(raw_hcv, filtered_hcv)

    assert expected_filtered_hcv == filtered_hcv.getvalue()
    assert expected_invalid_headers == invalid_headers


def test_filter_duplicates():
    raw_hcv = StringIO("""\
>Ref.1a.Foo
ACTAGGA
GAGATTT
>Ref.2b.Bar
TAGACT
>Ref.1a.Foo
AAA
""")
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.1a.Foo-1a
ACTAGGA
GAGATTT
>Ref.2b.Bar-2b
TAGACT
"""

    expected_invalid_headers = ['>Ref.1a.Foo (duplicate)']

    invalid_headers = filter_hcv_fasta(raw_hcv, filtered_hcv)

    assert expected_filtered_hcv == filtered_hcv.getvalue()
    assert expected_invalid_headers == invalid_headers


def test_filter_exclusions():
    raw_hcv = StringIO("""\
>Ref.1a.Foo
ACTAGGA
GAGATTT
>Ref.2b.Bar
TAGACT
""")
    excluded = ['Bar']
    filtered_hcv = StringIO()
    expected_filtered_hcv = """\
>Ref.1a.Foo-1a
ACTAGGA
GAGATTT
"""

    expected_invalid_headers = ['>Ref.2b.Bar (excluded)']

    invalid_headers = filter_hcv_fasta(raw_hcv, filtered_hcv, excluded)

    assert expected_filtered_hcv == filtered_hcv.getvalue()
    assert expected_invalid_headers == invalid_headers


def test_monophyletic():
    tree_source = "(foo-1a, (bar-2b, baz-2b));"
    report = StringIO()
    expected_report = ''

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_polyphyletic():
    tree_source = "((foo-1a, bar-2b), (floop-1a, baz-2c));"
    report = StringIO()
    expected_report = """\
1a polyphyletic bar-2b, baz-2c
1 polyphyletic bar-2b, baz-2c
2 polyphyletic floop-1a, foo-1a

      /-foo-1a
   /-|
  |   \-bar-2b
--|
  |   /-floop-1a
   \-|
      \-baz-2c
"""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_paraphyletic():
    tree_source = "(bar-2b, (baz-2c, (foo-1a, floop-1a)));"
    report = StringIO()
    expected_report = """\
2 paraphyletic floop-1a, foo-1a

   /-bar-2b
--|
  |   /-baz-2c
   \-|
     |   /-foo-1a
      \-|
         \-floop-1a
"""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_paraphyletic_genotype6():
    tree_source = "(bar-6b, (baz-6c, (foo-1a, floop-1a)));"
    report = StringIO()
    expected_report = ""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()
