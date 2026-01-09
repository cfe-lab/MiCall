
from pathlib import Path
import pandas as pd
import tomllib
from typing import Optional


def aggregate_runs_overlaps(input: Path, output: Path, properties: Optional[Path] = None) -> None:
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
