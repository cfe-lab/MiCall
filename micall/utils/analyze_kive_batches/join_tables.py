
from typing import Iterable, Set
from pathlib import Path
import pandas as pd

from micall.utils.new_atomic_file import new_atomic_text_file


def join_tables(inputs: Iterable[Path], column: str, output: Path) -> None:
    """
    Reads all CSVs in `inputs`, full-outer-joins them on `column`, and writes the result to `output`.

    Requirements:
      - Each CSV must have a header row.
      - `column` must be present in every CSV.
      - No two CSVs may share a non-index column name; if they do, a ValueError is raised.
      - Missing values in the output are written as empty strings.
    """
    inputs = list(inputs)
    if not inputs:
        raise ValueError("At least one input file must be provided.")

    # We'll collect DataFrames here
    dfs = []
    # To preserve the final column order, track non-index columns in the order seen
    col_order = []
    seen_cols: Set[str] = set()

    for path in inputs:
        df = pd.read_csv(path)

        if column not in df.columns:
            raise ValueError(f"Index column {column!r} not found in file {path}")

        # find non-index columns in this df
        this_cols = [c for c in df.columns if c != column]
        overlap = seen_cols.intersection(this_cols)
        if overlap:
            raise ValueError(
                "Column name overlap across inputs (excluding index column): "
                + ", ".join(sorted(overlap))
            )

        # record them for final ordering
        for c in this_cols:
            if c not in seen_cols:
                seen_cols.add(c)
                col_order.append(c)

        dfs.append(df)

    # Perform successive full outer joins
    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on=column, how="outer")

    # Reorder columns: first the index column, then all other columns in the order we saw them
    final_cols = [column] + col_order
    merged = merged.loc[:, final_cols]

    # Write out, with missing values as empty strings
    with new_atomic_text_file(output) as out:
        merged.to_csv(out, index=False, na_rep="")
