from pathlib import Path
import pandas as pd
from pandas.api.types import is_numeric_dtype


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1. Read
    df = pd.read_csv(input)

    # 2. (Optional) If you know certain columns must be numeric,
    #    you can coerce them here.  E.g.:
    # for c in ['size', 'depth', 'run_time']:
    #     df[c] = pd.to_numeric(df[c], errors='coerce')

    # 3. Split into the two apps, index by sample
    df1 = df[df['app'] == app1].set_index('sample')
    df2 = df[df['app'] == app2].set_index('sample')

    # 4. Keep only samples present in both
    common = df1.index.intersection(df2.index)
    df1 = df1.loc[common]
    df2 = df2.loc[common]

    # 5. Build the output
    out = pd.DataFrame(index=common)
    out['sample'] = out.index  # restore 'sample' as column

    # Loop over every column in the original
    for col in df.columns:
        if col == 'sample':
            continue

        left  = df1[col]
        right = df2[col]

        # if the original column was numeric → subtract
        if is_numeric_dtype(df[col]):
            out[col] = left - right

        # otherwise build "L/R" or just "L" if L==R
        else:
            ls = left.astype(str)
            rs = right.astype(str)
            out[col] = ls.where(ls == rs, ls + '/' + rs)

    # 6. Reorder to original column order and write out
    out = out[df.columns]
    out.to_csv(output, index=False)
