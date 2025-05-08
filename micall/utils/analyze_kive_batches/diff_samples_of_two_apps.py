from pathlib import Path
import pandas as pd
import numpy as np


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1. Read the CSV
    df = pd.read_csv(input)

    # 2. Split out the two apps
    df1 = df[df['app'] == app1]
    df2 = df[df['app'] == app2]

    # 3. Restrict to samples that appear in both
    common = set(df1['sample']).intersection(df2['sample'])
    df1 = df1[df1['sample'].isin(common)].set_index('sample')
    df2 = df2[df2['sample'].isin(common)].set_index('sample')

    # 4. Give them distinct suffixes and join on sample
    df1 = df1.add_suffix('_left')
    df2 = df2.add_suffix('_right')
    merged = df1.join(df2, how='inner')

    # 5. Build the result by iterating original columns
    result = pd.DataFrame()
    # bring sample back as a column
    result['sample'] = merged.index

    # figure out which columns were numeric in the original df
    numeric_cols = df.select_dtypes(include='number').columns

    for col in df.columns:
        if col == 'sample':
            continue
        left = merged[f"{col}_left"]
        right = merged[f"{col}_right"]

        if col in numeric_cols:
            # numeric → subtract
            result[col] = left - right
        else:
            # non‐numeric → "Left/Right", or just "Left" if equal
            lstr = left.astype(str)
            rstr = right.astype(str)
            result[col] = np.where(lstr == rstr, lstr, lstr + '/' + rstr)

    # 6. Reorder to original column order and write out
    result = result[df.columns]
    result.to_csv(output, index=False)
