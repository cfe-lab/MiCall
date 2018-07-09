from io import StringIO

from micall.utils.hcv_reference_tree import filter_hcv_fasta, check_tree, combine_samples


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
    tree_source = "(Ref.foo-1a, (Sample.bar-2b, Ref.baz-2b));"
    report = StringIO()
    expected_report = ''

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_polyphyletic():
    tree_source = "((Ref.foo-1a, Sample.bar-2b), (Sample.floop-1a, Ref.baz-2c));"
    report = StringIO()
    expected_report = """\
1a polyphyletic Ref.baz-2c, Sample.bar-2b
1 polyphyletic Ref.baz-2c, Sample.bar-2b
2 polyphyletic Ref.foo-1a, Sample.floop-1a
Sample.bar-2b treed with 1a
Sample.floop-1a treed with 2c

      /-Ref.foo-1a
   /-|
  |   \-Sample.bar-2b
--|
  |   /-Sample.floop-1a
   \-|
      \-Ref.baz-2c
"""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_paraphyletic():
    tree_source = "((Sample.bar-2b, Ref.byre-2b), (Ref.baz-2c, (Ref.foo-1a, Sample.floop-1a)));"
    report = StringIO()
    expected_report = """\
2 paraphyletic Ref.foo-1a, Sample.floop-1a

      /-Sample.bar-2b
   /-|
  |   \-Ref.byre-2b
--|
  |   /-Ref.baz-2c
   \-|
     |   /-Ref.foo-1a
      \-|
         \-Sample.floop-1a
"""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_paraphyletic2():
    tree_source = "((Ref.byre-2b, Sample.bar-2b), (Ref.baz-2c, (Ref.foo-1a, Sample.floop-1a)));"
    report = StringIO()
    expected_report = """\
2 paraphyletic Ref.foo-1a, Sample.floop-1a

      /-Ref.byre-2b
   /-|
  |   \-Sample.bar-2b
--|
  |   /-Ref.baz-2c
   \-|
     |   /-Ref.foo-1a
      \-|
         \-Sample.floop-1a
"""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_paraphyletic_genotype6():
    tree_source = "((Sample.bar-6b, Ref.byre-6b), (Ref.baz-6c, (Ref.foo-1a, Sample.floop-1a)));"
    report = StringIO()
    expected_report = ""

    check_tree(tree_source, report)

    assert expected_report == report.getvalue()


def test_combine_samples():
    filtered_hcv = StringIO("""\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
""")
    consensus_file = StringIO("""\
sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
2130A-HCV_S15,HCV-2a,15,MAX,7669,TCCATGTCATACTCCTGGACAGGGGCTCTG
2130A-HCV_S15,HCV-2a,15,0.010,7669,TCHATGTCATACTCCTGGACAGGGGCTCTG
""")
    coverage_scores = StringIO("""\
sample,seed,on.score
2130A-HCV_S15,HCV-2a,4
""")
    combined_file = StringIO()
    expected_combined_file = """\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
>Sample.2a.2130A-HCV_S15-2a
TCCATGTCATACTCCTGGACAGGGGCTCTG
"""

    combine_samples(filtered_hcv, consensus_file, coverage_scores, combined_file)

    assert expected_combined_file == combined_file.getvalue()


def test_combine_samples_low_coverage():
    filtered_hcv = StringIO("""\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
""")
    consensus_file = StringIO("""\
sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
2000A-HCV_S14,HCV-3c,15,MAX,7669,TCCATGTCATACTCCTGGACAGGGGCTCTG
2130A-HCV_S15,HCV-2a,15,MAX,7669,TCCATGTCATACTCCTGGACAGGGGCTCTG
2130A-HCV_S15,HCV-3c,15,MAX,7669,TCH-----ATACTCCTGGACAGGGGCTCTG
""")
    coverage_scores = StringIO("""\
sample,seed,on.score
2000A-HCV_S14,HCV-3c,4
2130A-HCV_S15,HCV-2a,4
2130A-HCV_S15,HCV-3c,1
""")
    combined_file = StringIO()
    expected_combined_file = """\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
>Sample.3c.2000A-HCV_S14-3c
TCCATGTCATACTCCTGGACAGGGGCTCTG
>Sample.2a.2130A-HCV_S15-2a
TCCATGTCATACTCCTGGACAGGGGCTCTG
"""

    combine_samples(filtered_hcv, consensus_file, coverage_scores, combined_file)

    assert expected_combined_file == combined_file.getvalue()


def test_combine_samples_hcv_only():
    filtered_hcv = StringIO("""\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
""")
    consensus_file = StringIO("""\
sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
2130A-HIV_S15,HIV1-B-FR-KF716496-seed,15,MAX,7669,TCCATGTCATACTCCTGGACAGGGGCTCTG
2130A-HIV_S15,HIV1-B-FR-KF716496-seed,15,0.010,7669,TCHATGTCATACTCCTGGACAGGGGCTCTG
""")
    coverage_scores = StringIO("""\
sample,seed,on.score
2130A-HCV_S15,HCV-2a,4
""")
    combined_file = StringIO()
    expected_combined_file = """\
>Ref.1a.Foo-1a
ACTACCTGA
TGACG
"""

    combine_samples(filtered_hcv, consensus_file, coverage_scores, combined_file)

    assert expected_combined_file == combined_file.getvalue()
