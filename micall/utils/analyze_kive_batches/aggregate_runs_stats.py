
from pathlib import Path
import pandas as pd
import tomllib
from typing import Optional


def aggregate_runs_stats(input: Path, output: Path, properties: Optional[Path] = None) -> None:
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
            avg_s_alignment_score = ('stitched_alignment_score', 'mean'),
            avg_mlen              = ('mlen',              'mean'),
            avg_total_mlen        = ('total_mlen',        'mean'),
            avg_overlap_count     = ('overlap_count',     'mean'),
            avg_number_of_contigs = ('number_of_contigs', 'mean'),
            avg_contigs_size      = ('avg_contigs_size',  'mean'),
            avg_soft_clips_count  = ('soft_clips_count',  'mean'),
            avg_exact_uncovered   = ('exact_uncovered',   'mean'),
            avg_run_time          = ('run_time',          'mean'),
            run_count             = ('app',               'count'),
        )
        .reset_index()
    )

    # 3. Sort by properties.toml order if provided
    if properties is not None:
        with properties.open("rb") as reader:
            props_obj = tomllib.load(reader)
        app_order = list(props_obj.keys())
        # Create a categorical type with the order from properties.toml
        grouped['app'] = pd.Categorical(grouped['app'], categories=app_order, ordered=True)
        grouped = grouped.sort_values('app')

    # 4. Write out to CSV
    grouped.to_csv(output, index=False)
