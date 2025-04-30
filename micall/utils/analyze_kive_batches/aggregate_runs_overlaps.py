
from pathlib import Path
import pandas as pd


def aggregate_runs_overlaps(input: Path, output: Path) -> None:
    # 1. Read the CSV
    df = pd.read_csv(input)

    # 2. Group by "app" and compute the stats.
    grouped = (
        df
        .groupby('app')
        .agg(
            avg_overlap_size       = ('overlap_size',       'mean'),
            avg_overlap_mismatches = ('overlap_mismatches', 'mean'),
            avg_overlap_pvalue     = ('overlap_pvalue',     'mean'),
            overlap_count          = ('app',                'count'),
        )
        .reset_index()
    )

    # 3. Write out to CSV
    grouped.to_csv(output, index=False)
