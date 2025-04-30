import csv
from typing import Iterable, Set, Dict
from pathlib import Path
from micall.utils.new_atomic_file import new_atomic_text_file


def join_tables(inputs: Iterable[Path], column: str, output: Path) -> None:
    """
    Full‐outer‐join all CSVs in `inputs` on `column`, writing to `output`.

    - Each CSV must have a header row.
    - `column` must be present in every CSV.
    - No two CSVs may share a non‐index column name (raises ValueError).
    - Missing values become empty strings.
    """

    inputs = tuple(inputs)
    if not inputs:
        raise ValueError("At least one input file must be provided.")

    # For each input file we will build a map: index_value -> {non_index_col: cell, ...}
    per_file_maps: list[Dict[str, Dict[str, str]]] = []
    seen_cols: Set[str] = set()      # all non-index cols seen so far
    col_order: list[str] = []        # their order

    for path in inputs:
        with path.open(newline='') as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError(f"No header found in {path}")
            if column not in reader.fieldnames:
                raise ValueError(f"Index column {column!r} not in {path}")

            # figure out this file's non-index columns
            this_cols = [c for c in reader.fieldnames if c != column]
            overlap = seen_cols.intersection(this_cols)
            if overlap:
                raise ValueError(
                    "Column name overlap across inputs (excluding index column): "
                    + ", ".join(sorted(overlap))
                )

            # record new columns in order
            for c in this_cols:
                seen_cols.add(c)
                col_order.append(c)

            # build map index_value -> row-dict-of-non-index-cols
            m: Dict[str, Dict[str, str]] = {}
            for row in reader:
                key = row[column]
                # extract only the non-index columns
                m[key] = {c: row[c] for c in this_cols}
            per_file_maps.append(m)

    # collect all index-values in the order we see them
    all_index_values: list[str] = []
    seen_idx: Set[str] = set()
    for m in per_file_maps:
        for key in m:
            if key not in seen_idx:
                seen_idx.add(key)
                all_index_values.append(key)

    # the full output header
    final_cols = [column] + col_order

    # write out via DictWriter, filling missing cells with ""
    with new_atomic_text_file(output) as out_f:
        writer = csv.DictWriter(out_f, fieldnames=final_cols)
        writer.writeheader()

        for key in all_index_values:
            # start every row with empty strings
            this: Dict[str, str] = {col: "" for col in final_cols}
            this[column] = key
            # fill in any values from each file
            for m in per_file_maps:
                vals = m.get(key)
                if vals:
                    this.update(vals)
            writer.writerow(this)
