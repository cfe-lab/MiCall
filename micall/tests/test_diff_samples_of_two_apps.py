import csv
import pandas as pd
import pytest
from pathlib import Path
import logging

from micall.utils.user_error import UserError
from micall.utils.analyze_kive_batches.diff_samples_of_two_apps import \
    diff_samples_of_two_apps


def write_csv(path: Path, header, rows):
    """Utility: write a CSV by header+rows."""
    with open(path, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(header)
        w.writerows(rows)


def read_df(path: Path):
    """Utility: read with pandas for easy assertions."""
    return pd.read_csv(path, dtype=str)


def assert_df_equal_ignore_dtype(actual: pd.DataFrame, expected: pd.DataFrame):
    """Compare two small DataFrames by content only."""
    pd.testing.assert_frame_equal(actual.sort_index(axis=1),
                                  expected.sort_index(axis=1),
                                  check_dtype=False,
                                  check_like=True)


def test_simple_pair(tmp_path):
    # Example 1:
    # sample,app,size,type
    # 1,bob,50,x
    # 1,alice,70,y
    inp = tmp_path / "in.csv"
    out = tmp_path / "out.csv"
    header = ["sample", "app", "size", "type"]
    rows = [
        ("1", "bob",   "50", "x"),
        ("1", "alice", "70", "y"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)

    # Expect one row: 1,alice/bob,70-50=20,y/x + base columns
    want = pd.DataFrame([{
        "sample": "1",
        "app":    "alice/bob",
        "size":   "20",
        "type":   "y/x",
        "size_base": "70",
        "type_base": "y",
    }])
    got = read_df(out)
    assert_df_equal_ignore_dtype(got, want)


def test_three_apps_only_two_selected(tmp_path):
    # Example 2: three apps but we pick only bob vs alice
    inp = tmp_path / "in2.csv"
    out = tmp_path / "out2.csv"
    header = ["sample", "app", "size", "type"]
    rows = [
        ("1", "bob",   "50", "x"),
        ("1", "celia", "60", "y"),
        ("1", "alice", "70", "z"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)

    want = pd.DataFrame([{
        "sample": "1",
        "app":    "alice/bob",
        "size":   "20",
        "type":   "z/x",
        "size_base": "70",
        "type_base": "z",
    }])
    got = read_df(out)
    assert_df_equal_ignore_dtype(got, want)


def test_multiple_rows_per_sample(tmp_path):
    # Example 3: two rows for bob, two rows for alice
    inp = tmp_path / "in3.csv"
    out = tmp_path / "out3.csv"
    header = ["sample", "app", "size", "type"]
    rows = [
        ("1", "bob",   "50", "x"),
        ("1", "bob",   "50", "m"),
        ("1", "alice", "70", "z"),
        ("1", "alice", "70", "k"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)

    # We get two output rows in the same order pairs were seen:
    want = pd.DataFrame([
        {"sample": "1", "app": "alice/bob", "size": "20", "type": "z+k/x+m", "size_base": "70", "type_base": "z+k"},
    ])
    got = read_df(out)
    assert_df_equal_ignore_dtype(got, want)


def test_two_samples(tmp_path):
    # Two distinct samples, both shared
    inp = tmp_path / "in4.csv"
    out = tmp_path / "out4.csv"
    header = ["sample", "app", "size"]
    rows = [
        ("1", "bob",   "10"),
        ("1", "alice", "20"),
        ("2", "bob",   "30"),
        ("2", "alice", "40"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)

    want = pd.DataFrame([
        {"sample": "1", "app": "alice/bob", "size": "10", "size_base": "20"},
        {"sample": "2", "app": "alice/bob", "size": "10", "size_base": "40"},
    ])
    got = read_df(out)
    assert_df_equal_ignore_dtype(got, want)


def test_no_common_samples_yields_only_header(tmp_path):
    inp = tmp_path / "in_empty.csv"
    out = tmp_path / "out_empty.csv"
    header = ["sample", "app", "size"]
    rows = [
        ("1", "bob", "5"),
        ("2", "alice", "6"),
    ]
    write_csv(inp, header, rows)

    # bob vs alice share no sample → empty data
    diff_samples_of_two_apps(inp, app1="bob", app2="alice", output=out)

    # read with pandas: zero rows but correct header
    df = pd.read_csv(out)
    expected_cols = header + ["size_base"]  # base columns added
    assert list(df.columns) == expected_cols
    assert df.shape[0] == 0


def test_missing_app_raises(tmp_path):
    inp = tmp_path / "in_missing.csv"
    out = tmp_path / "out_missing.csv"
    header = ["sample", "app", "size"]
    rows = [
        ("1", "bob", "5"),
    ]
    write_csv(inp, header, rows)

    with pytest.raises(UserError) as ei:
        diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    msg = str(ei.value)
    assert "App(s) not found" in msg
    assert "alice" in msg


def test_coerce_and_nan(tmp_path):
    # size is supposed to be numeric, but we give one bad entry
    inp = tmp_path / "in_nan.csv"
    out = tmp_path / "out_nan.csv"
    header = ["sample", "app", "size"]
    rows = [
        ("1", "bob",   "foo"),   # non‐numeric
        ("1", "alice", "10"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = pd.read_csv(out)

    # foo → NaN, so result is 10/foo, base is 10
    assert df.loc[0, "size"] == "10/foo"
    assert str(df.loc[0, "size_base"]) == "10"
    # still has the correct header
    expected_cols = header + ["size_base"]
    assert list(df.columns) == expected_cols


def test_identical_non_numeric(tmp_path):
    # test that if a non‐numeric column is the same,
    # we do not insert slash.
    inp = tmp_path / "in_same.csv"
    out = tmp_path / "out_same.csv"
    header = ["sample", "app", "tag"]
    rows = [
        ("1", "bob",   "X"),
        ("1", "alice", "X"),   # same tag
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = pd.read_csv(out, dtype=str)
    # tag is the same → we expect "X" not "X/X", base is "X"
    assert df.loc[0, "tag"] == "X"
    assert df.loc[0, "tag_base"] == "X"
    assert df.loc[0, "app"] == "alice/bob"


def test_unbalanced_rows_pairs_and_logs_warning(tmp_path, caplog):
    """
    If app1 has 3 rows for a sample but app2 only 2,
    we should emit a warning and only produce 2 output rows.
    """
    caplog.set_level(logging.WARNING)
    inp = tmp_path / "unbal.csv"
    out = tmp_path / "unbal_out.csv"
    header = ["sample", "app", "size", "type"]
    rows = [
        ("1", "bob",   "10", "X"),
        ("1", "bob",   "20", "Y"),
        ("1", "bob",   "30", "Z"),
        ("1", "alice", "15", "A"),
        ("1", "alice", "25", "B"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = read_df(out)

    # only the first pair shows up
    assert len(df) == 1
    # sizes: 15-10=5, 25-20=5, but aggregated so mean of (15,25) - mean of (10,20,30) = 20-20=0
    assert set(df["size"]) == {'0'}
    # types: A/X and B/Y, but aggregated so A+B/X+Y+Z
    assert set(df["type"]) == {"A+B/X+Y+Z"}
    # base columns should have alice's aggregated values
    assert set(df["size_base"]) == {'20'}
    assert set(df["type_base"]) == {"A+B"}


def test_float_subtraction_and_negative(tmp_path):
    """
    Test with floats and ensure we get correct decimal subtraction
    (and negative results if app2 > app1).
    """
    inp = tmp_path / "floats.csv"
    out = tmp_path / "floats_out.csv"
    header = ["sample", "app", "val"]
    rows = [
        ("1", "bob",   "3.5"),
        ("1", "alice", "2.1"),
        ("2", "bob",   "1.0"),
        ("2", "alice", "5.0"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = read_df(out)

    # For sample=1: 2.1-3.5 = -1.4, base=2.1
    # For sample=2: 5.0-1.0 = 4.0, base=5.0
    got = dict(zip(df["sample"], df["val"]))
    got_base = dict(zip(df["sample"], df["val_base"]))
    assert got["1"] == "-1.4"
    assert got["2"] == "4"
    assert got_base["1"] == "2.1"
    assert got_base["2"] == "5"


def test_column_order_is_preserved(tmp_path):
    """
    If the input columns occur in a custom order, the output
    must preserve exactly that order in the CSV.
    """
    inp = tmp_path / "order.csv"
    out = tmp_path / "order_out.csv"
    header = ["type", "sample", "foo", "app", "size"]
    rows = [
        ("X", "1", "alpha", "bob",   "100"),
        ("Y", "1", "beta",  "alice", "110"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    # read only header line - should include base columns
    actual_header = out.read_text().splitlines()[0].split(',')
    expected_header = header + ["type_base", "foo_base", "size_base"]
    assert actual_header == expected_header

    # check the single data row
    df = read_df(out)
    rec = df.iloc[0].to_dict()

    # NOTE: type: Y/X, foo: beta/alpha, app: alice/bob, size:10
    # base columns have alice's values
    assert rec == {
        "type":   "Y/X",
        "sample": "1",
        "foo":    "beta/alpha",
        "app":    "alice/bob",
        "size":   "10",
        "type_base": "Y",
        "foo_base": "beta",
        "size_base": "110",
    }


def test_boolean_and_nonstring_columns(tmp_path):
    """
    A boolean column should be treated like a non‐numeric
    and joined with slash if different.
    """
    inp = tmp_path / "bool.csv"
    out = tmp_path / "bool_out.csv"
    header = ["sample", "app", "flag"]
    rows = [
        ("1", "bob",   "True"),
        ("1", "alice", "False"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = read_df(out)
    assert df.loc[0, "flag"] == "False/True"
    assert df.loc[0, "flag_base"] == "False"
    assert df.loc[0, "app"]  == "alice/bob"


def test_sample_codes_are_preserved_as_strings(tmp_path):
    """
    If the sample column is alphanumeric, we still
    use it verbatim in the output.
    """
    inp = tmp_path / "alpha.csv"
    out = tmp_path / "alpha_out.csv"
    header = ["sample", "app", "size"]
    rows = [
        ("S1", "bob",   "5"),
        ("S1", "alice", "8"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = read_df(out)
    assert df.loc[0, "sample"] == "S1"
    assert df.loc[0, "size"]   == "3"
    assert df.loc[0, "size_base"] == "8"


def test_same_app_diffing_itself(tmp_path):
    """
    If app1 == app2, we still pair each row with itself:
    numeric diff → 0, non‐numeric → same value.
    """
    inp = tmp_path / "self.csv"
    out = tmp_path / "self_out.csv"
    header = ["sample", "app", "size", "tag"]
    rows = [
        ("1", "bob", "10", "X"),
        ("1", "bob", "20", "Y"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="bob", app2="bob", output=out)
    df = read_df(out)
    # we should get one row with aggregated values and base columns
    assert list(df["size"]) == ["10"]
    assert list(df["tag"])  == ["X+Y"]
    assert list(df["app"])  == ["bob"]
    # base columns should have the same values since it's same app
    assert list(df["size_base"]) == ["15"]  # mean of 10,20
    assert list(df["tag_base"])  == ["X+Y"]


def test_same_app_diffing_itself2(tmp_path):
    """
    If app1 == app2, we still pair each row with itself:
    numeric diff → 0, non‐numeric → same value.
    """
    inp = tmp_path / "self.csv"
    out = tmp_path / "self_out.csv"
    header = ["sample", "app", "size", "tag", "tag2"]
    rows = [
        ("1", "bob", "10", "X", "a"),
        ("1", "bob", "20", "Y", "b"),
        ("1", "bob", "30", "X", "c"),
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="bob", app2="bob", output=out)
    df = read_df(out)
    # we should get one row with aggregated values and base columns
    assert list(df["size"]) == ["10"]  # variance: average diff between consecutive values (|20-10| + |30-20|)/(3-1) = 20/2 = 10
    assert list(df["tag"])  == ["X+Y+X"]
    assert list(df["tag2"])  == ["a+b+c"]
    assert list(df["app"])  == ["bob"]
    # base columns should have the same values since it's same app
    assert list(df["size_base"]) == ["20"]  # mean of 10,20,30
    assert list(df["tag_base"])  == ["X+Y+X"]
    assert list(df["tag2_base"])  == ["a+b+c"]


def test_empty_input_header_only_raises(tmp_path):
    """
    If the CSV has no data rows at all, both apps are missing
    and we should get a UserError.
    """
    inp = tmp_path / "empty.csv"
    out = tmp_path / "empty_out.csv"
    # write header only
    write_csv(inp, ["sample", "app", "size"], [])

    with pytest.raises(UserError):
        diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)


def test_missing_non_numeric_values(tmp_path):
    """
    If a non‐numeric column is empty on one side but present on the other,
    it will be slash‐joined (e.g. '/value' or 'value/').
    """
    inp = tmp_path / "miss.csv"
    out = tmp_path / "miss_out.csv"
    header = ["sample", "app", "tag"]
    rows = [
        ("1", "bob",   ""),        # empty → read as empty string
        ("1", "alice", "HELLO"),   # non‐empty
    ]
    write_csv(inp, header, rows)

    diff_samples_of_two_apps(inp, app1="alice", app2="bob", output=out)
    df = read_df(out)
    # tag = "HELLO/" because bob's side was empty, base = "HELLO"
    assert df.loc[0, "tag"] == "HELLO/"
    assert df.loc[0, "tag_base"] == "HELLO"
