from pathlib import Path
import pandas as pd
from pandas.api.types import is_numeric_dtype


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1) Read everything
    df = pd.read_csv(input)

    # 2) Split into the two apps
    df1 = df[df['app'] == app1]
    df2 = df[df['app'] == app2]

    # 3) Find the set of samples that appear in both
    common_samples = sorted(
        set(df1['sample']).intersection(df2['sample'])
    )

    # 4) Prepare buckets by sample
    grouped1 = { s: grp.reset_index(drop=True)
                 for s, grp in df1.groupby('sample') if s in common_samples }
    grouped2 = { s: grp.reset_index(drop=True)
                 for s, grp in df2.groupby('sample') if s in common_samples }

    # 5) Prepare output records
    out_recs = []
    cols = list(df.columns)          # preserve original column order
    for sample in common_samples:
        A = grouped1[sample]
        B = grouped2[sample]
        n = min(len(A), len(B))

        if len(A) != len(B):
            print(f"Warning: sample {sample} has {len(A)} rows in {app1} "
                  f"but {len(B)} rows in {app2}; using first {n} pairs.")

        for i in range(n):
            rowA = A.iloc[i]
            rowB = B.iloc[i]
            rec = {'sample': sample}

            for col in cols:
                if col == 'sample':
                    continue

                a = rowA[col]
                b = rowB[col]

                # Numeric? subtract
                if is_numeric_dtype(df[col]):
                    try:
                        rec[col] = a - b
                    except Exception:
                        # fallback if dtype inference was wonky
                        rec[col] = pd.to_numeric(a, errors='coerce') \
                                 - pd.to_numeric(b, errors='coerce')
                else:
                    # Non‐numeric: "L/R" or just "L" if equal
                    sa = str(a)
                    sb = str(b)
                    rec[col] = sa if sa == sb else f"{sa}/{sb}"

            out_recs.append(rec)

    # 6) Build DataFrame and write CSV
    out_df = pd.DataFrame(out_recs)

    # ensure same column order
    out_df = out_df[cols]
    out_df.to_csv(output, index=False)
