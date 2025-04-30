
from pathlib import Path
import pandas as pd


def aggregate_runs_stats(input: Path, output: Path) -> None:
    # 1. Read the CSV
    df = pd.read_csv(input)

    # 2. Group by "assembler" and compute the same stats you had in SQL:
    #    AVG(concordance), AVG(depth), AVG(mlen),
    #    AVG(total_mlen), AVG(overlap_count),
    #    AVG(number_of_contigs), AVG(avg_contigs_size),
    #    AVG(run_time), SUM(run_time), COUNT(assembler)
    #
    # We can use Pandas "named aggregation" (pandas>=0.25):
    grouped = (
        df
        .groupby('assembler')
        .agg(
            avg_concordance       = ('concordance',       'mean'),
            avg_depth             = ('depth',             'mean'),
            avg_mlen              = ('mlen',              'mean'),
            avg_total_mlen        = ('total_mlen',        'mean'),
            avg_overlap_count     = ('overlap_count',     'mean'),
            avg_number_of_contigs = ('number_of_contigs', 'mean'),
            avg_contigs_size      = ('avg_contigs_size',  'mean'),
            avg_run_time          = ('run_time',          'mean'),
            count                 = ('assembler',         'count'),
        )
        .reset_index()
    )

    # 3. Write out to CSV
    grouped.to_csv(output, index=False)
