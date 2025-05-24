
from pathlib import Path
import pandas as pd


def aggregate_runs_stats(input: Path, output: Path) -> None:
    # 1. Read the CSV
    df = pd.read_csv(input)

    # 2. Group by "app" and compute the stats.
    grouped = (
        df
        .groupby('app')
        .agg(
            avg_concordance       = ('concordance',       'mean'),
            avg_depth             = ('depth',             'mean'),
            avg_alignment_score   = ('alignment_score',   'mean'),
            avg_mlen              = ('mlen',              'mean'),
            avg_total_mlen        = ('total_mlen',        'mean'),
            avg_overlap_count     = ('overlap_count',     'mean'),
            avg_number_of_contigs = ('number_of_contigs', 'mean'),
            avg_contigs_size      = ('avg_contigs_size',  'mean'),
            avg_run_time          = ('run_time',          'mean'),
            run_count             = ('app',               'count'),
        )
        .reset_index()
    )

    # 3. Write out to CSV
    grouped.to_csv(output, index=False)
