from pathlib import Path
from typing import Callable, Union
import pandas as pd

from micall.utils.user_error import UserError
from .logger import logger


def is_numeric_dtype(col: pd.Series) -> bool:
    return col.dtype.kind in ('i','u','f')


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1) Read everything
    df = pd.read_csv(input)
    cols = tuple(df.columns)  # keep the original column order

    # 2) Fail fast if either app is not present
    apps_present = set(df['app'])
    missing = {app1, app2} - apps_present
    if missing:
        raise UserError("App(s) not found in input: %s.", ', '.join(map(repr, missing)))

    # 3) Split into the two apps
    df1 = df[df['app'] == app1]
    df2 = df[df['app'] == app2]

    def tostr(x: pd.DataFrame) -> str:
        if pd.isna(x):
            return ''
        else:
            return str(x)

    # 4) If there are multiple rows per (sample, app), collapse them:
    #    - numeric cols → mean
    #    - non-numeric cols → concat all values
    def concat_all(xs: pd.Series) -> str:
        return '+'.join(map(tostr, xs.tolist()))

    def diff_all(xs: pd.Series) -> float:
        if len(xs) == 1:
            return 0

        diffs = (abs(a - b) for a, b in zip(xs, xs[1:]))
        return sum(diffs) / (len(xs) - 1)

    agg_map: dict[str, Union[str, Callable[[pd.Series], Union[float, str]]]] = {}
    # For same-app comparison, we need separate aggregation for base columns
    base_agg_map: dict[str, Union[str, Callable[[pd.Series], Union[float, str]]]] = {}

    for c in cols:
        if c == 'sample':
            continue
        elif c == 'app':
            agg_map[c] = 'first'
            base_agg_map[c] = 'first'
        elif is_numeric_dtype(df[c]):
            if app1 == app2:
                agg_map[c] = diff_all
                base_agg_map[c] = 'mean'  # base columns show the actual mean values
            else:
                agg_map[c] = 'mean'
                base_agg_map[c] = 'mean'
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
                column = df[col]

                if app1 == app2:
                    rec[col] = a
                    # For same app comparison, also add base column
                    if col != 'app':
                        rec[f"{col}_base"] = a_base  # Use the mean/aggregated base value
                    continue

                # only treat ints and floats as numeric; exclude bools
                if is_numeric_dtype(column):
                    # numeric diff
                    rec[col] = pd.to_numeric(a, errors='coerce') - \
                               pd.to_numeric(b, errors='coerce')
                    if float.is_integer(rec[col]):
                        rec[col] = int(rec[col])
                else:
                    # non‐numeric: "L/R" or just "L" if equal
                    sa, sb = tostr(a), tostr(b)
                    rec[col] = sa if sa == sb else f"{sa}/{sb}"

                # Add base column (app1 value) for all columns except 'app'
                if col != 'app':
                    if is_numeric_dtype(column):
                        base_val = pd.to_numeric(a_base, errors='coerce')
                        if not pd.isna(base_val) and float.is_integer(base_val):
                            rec[f"{col}_base"] = int(base_val)
                        else:
                            rec[f"{col}_base"] = base_val
                    else:
                        rec[f"{col}_base"] = tostr(a_base)

            out_recs.append(rec)

    # 7) Build the final DataFrame with base columns added
    # Create new column order: original columns + base columns (excluding 'sample' and 'app')
    base_cols = [f"{col}_base" for col in cols if col not in ('sample', 'app')]
    new_cols = list(cols) + base_cols

    out_df = pd.DataFrame(out_recs, columns=new_cols)

    # 8) Write to CSV
    out_df.to_csv(output, index=False)
