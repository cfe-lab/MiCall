from pathlib import Path
from typing import Callable, Union
import pandas as pd

from micall.utils.user_error import UserError
from .logger import logger


def _can_parse_float(val: str) -> bool:
    try:
        float(val)
        return True
    except Exception:
        return False


def _parse_float(val: str) -> Union[float, None]:
    try:
        return float(val)
    except Exception:
        return None


def _is_numeric_by_values(col: pd.Series) -> bool:
    """Decide if a column is numeric based on its string values.

    Treat empty strings as missing; all other non-empty values must be
    convertible to float for the column to be considered numeric.
    """
    # Ensure we're working with Python strings; pandas shouldn't coerce types
    vals = ("" if v is None else str(v) for v in col.tolist())
    non_empty = [v for v in vals if v != ""]
    if not non_empty:
        # No data → not considered numeric to avoid producing numbers
        return False
    return all(_can_parse_float(v) for v in non_empty)


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1) Read everything as strings; disable NA/NaN inference entirely
    df = pd.read_csv(
        input,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    cols = tuple(df.columns)  # keep the original column order

    # 2) Fail fast if either app is not present
    apps_present = set(df['app'])
    missing = {app1, app2} - apps_present
    if missing:
        raise UserError("App(s) not found in input: %s.", ', '.join(map(repr, missing)))

    # 3) Split into the two apps
    df1 = df[df['app'] == app1]
    df2 = df[df['app'] == app2]

    def tostr(x) -> str:
        if x is None:
            return ''
        if isinstance(x, float):
            try:
                if x.is_integer():
                    return str(int(x))
            except Exception:
                pass
            # Use default string conversion for floats
            return str(x)
        return str(x)

    # 4) If there are multiple rows per (sample, app), collapse them:
    #    - numeric cols → mean
    #    - non-numeric cols → concat all values
    def concat_all(xs: pd.Series) -> str:
        return '+'.join(map(tostr, xs.tolist()))

    def _format_number(n: Union[float, int, None]) -> str:
        if n is None:
            return ''
        if isinstance(n, float):
            try:
                if n.is_integer():
                    return str(int(n))
            except Exception:
                pass
        return str(n)

    def mean_as_str(xs: pd.Series) -> str:
        vals = [v for v in (tostr(v) for v in xs.tolist()) if v != '']
        nums = [float(v) for v in vals if _can_parse_float(v)]
        if not nums:
            return ''
        mean_val = sum(nums) / len(nums)
        return _format_number(mean_val)

    def diff_all_as_str(xs: pd.Series) -> str:
        # Calculate average absolute difference across consecutive values (numeric only)
        vals = [v for v in (tostr(v) for v in xs.tolist()) if v != '']
        nums = [float(v) for v in vals if _can_parse_float(v)]
        if len(nums) <= 1:
            return '0'
        diffs = [abs(a - b) for a, b in zip(nums, nums[1:])]
        avg = sum(diffs) / len(diffs)
        return _format_number(avg)

    agg_map: dict[str, Union[str, Callable[[pd.Series], str]]] = {}
    # For same-app comparison, we need separate aggregation for base columns
    base_agg_map: dict[str, Union[str, Callable[[pd.Series], str]]] = {}

    # Pre-compute which columns appear numeric based on their string values
    numeric_cols: dict[str, bool] = {
        c: _is_numeric_by_values(df[c]) if c not in ('sample', 'app') else False
        for c in cols
    }

    for c in cols:
        if c == 'sample':
            continue
        elif c == 'app':
            agg_map[c] = 'first'
            base_agg_map[c] = 'first'
        elif numeric_cols.get(c, False):
            if app1 == app2:
                agg_map[c] = diff_all_as_str
                base_agg_map[c] = mean_as_str  # base columns show the actual mean values
            else:
                agg_map[c] = mean_as_str
                base_agg_map[c] = mean_as_str
        else:
            agg_map[c] = concat_all
            base_agg_map[c] = concat_all

    df1 = df1.groupby('sample', as_index=False).agg(agg_map)
    df2 = df2.groupby('sample', as_index=False).agg(agg_map)

    # For same-app comparison, create base DataFrames with mean aggregation
    if app1 == app2:
        df1_base = df[df['app'] == app1].groupby('sample', as_index=False).agg(base_agg_map)
    else:
        df1_base = df1

    # 4) Find the set of samples that appear in both
    common_samples = sorted(
        set(df1['sample']).intersection(df2['sample'])
    )

    # 5) Group by sample, reset index so we can pair‐up later
    grouped1 = { s: grp.reset_index(drop=True)
                 for s, grp in df1.groupby('sample') if s in common_samples }
    grouped2 = { s: grp.reset_index(drop=True)
                 for s, grp in df2.groupby('sample') if s in common_samples }

    # Also group base DataFrames for base column values
    grouped1_base = { s: grp.reset_index(drop=True)
                      for s, grp in df1_base.groupby('sample') if s in common_samples }

    # 6) Build output records by pairing rows of each sample
    out_recs = []
    for sample in common_samples:
        A = grouped1[sample]
        B = grouped2[sample]
        A_base = grouped1_base[sample]
        n = min(len(A), len(B))

        if len(A) != len(B):
            logger.warning("Sample %r has %s rows in %r "
                           "but %s rows in %r; using first %s pairs.",
                           str(sample), len(A), app1, len(B), app2, n,
                           )

        for i in range(n):
            rowA = A.iloc[i]
            rowB = B.iloc[i]
            rowA_base = A_base.iloc[min(i, len(A_base)-1)]  # Handle potential size difference
            rec = {'sample': sample}

            for col in cols:
                if col == 'sample':
                    continue

                a = rowA[col]
                b = rowB[col]
                a_base = rowA_base[col]  # Use base value for base columns
                column_is_numeric = numeric_cols.get(col, False)

                if app1 == app2:
                    rec[col] = tostr(a)
                    # For same app comparison, also add base column
                    if col != 'app':
                        rec[f"{col}_base"] = tostr(a_base)  # Use the mean/aggregated base value
                    continue

                # only treat as numeric when both sides are parseable numbers
                if column_is_numeric and _can_parse_float(str(a)) and _can_parse_float(str(b)):
                    aval = _parse_float(str(a))
                    bval = _parse_float(str(b))
                    if aval is not None and bval is not None:
                        diffval = aval - bval
                        rec[col] = _format_number(diffval)
                    else:
                        # Fallback to string diff if any side failed
                        sa, sb = tostr(a), tostr(b)
                        rec[col] = sa if sa == sb else f"{sa}/{sb}"
                else:
                    # non‐numeric: "L/R" or just "L" if equal
                    sa, sb = tostr(a), tostr(b)
                    rec[col] = sa if sa == sb else f"{sa}/{sb}"

                # Add base column (app1 value) for all columns except 'app'
                if col != 'app':
                    if column_is_numeric and _can_parse_float(str(a_base)):
                        base_val = _parse_float(str(a_base))
                        rec[f"{col}_base"] = _format_number(base_val)
                    else:
                        rec[f"{col}_base"] = tostr(a_base)

            out_recs.append(rec)

    # 7) Build the final DataFrame with base columns added
    # Create new column order: original columns + base columns (excluding 'sample' and 'app')
    base_cols = [f"{col}_base" for col in cols if col not in ('sample', 'app')]
    new_cols = list(cols) + base_cols

    out_df = pd.DataFrame(out_recs, columns=new_cols)

    # 8) Write to CSV
    # Ensure all values are strings, and replace any accidental NaN with empty string
    out_df = out_df.fillna('').astype(str)
    out_df.to_csv(output, index=False)
