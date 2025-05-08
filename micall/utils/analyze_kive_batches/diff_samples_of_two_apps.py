from pathlib import Path
import pandas as pd
from pandas.api.types import is_numeric_dtype

from micall.utils.user_error import UserError
from .logger import logger


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

    # 4) Find the set of samples that appear in both
    common_samples = sorted(
        set(df1['sample']).intersection(df2['sample'])
    )

    # 5) Group by sample, reset index so we can pair‐up later
    grouped1 = { s: grp.reset_index(drop=True)
                 for s, grp in df1.groupby('sample') if s in common_samples }
    grouped2 = { s: grp.reset_index(drop=True)
                 for s, grp in df2.groupby('sample') if s in common_samples }

    # 6) Build output records by pairing rows of each sample
    out_recs = []
    for sample in common_samples:
        A = grouped1[sample]
        B = grouped2[sample]
        n = min(len(A), len(B))

        if len(A) != len(B):
            logger.warning("Sample %r has %s rows in %r "
                           "but %s rows in %r; using first %s pairs.",
                           str(sample), len(A), app1, len(B), app2, n,
                           )

        for i in range(n):
            rowA = A.iloc[i]
            rowB = B.iloc[i]
            rec = {'sample': sample}

            for col in cols:
                if col == 'sample':
                    continue

                a = rowA[col]
                b = rowB[col]

                if is_numeric_dtype(df[col]):
                    # numeric diff
                    rec[col] = pd.to_numeric(a, errors='coerce') - \
                               pd.to_numeric(b, errors='coerce')
                else:
                    # non‐numeric: "L/R" or just "L" if equal
                    sa, sb = str(a), str(b)
                    rec[col] = sa if sa == sb else f"{sa}/{sb}"

            out_recs.append(rec)

    # 7) Build the final DataFrame, forcing the original columns even if out_recs=[]
    out_df = pd.DataFrame(out_recs, columns=cols)

    # 8) Write to CSV
    out_df.to_csv(output, index=False)
