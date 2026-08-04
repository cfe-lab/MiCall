import csv
import itertools
from typing import Iterable, List, Dict, Set
from pathlib import Path
from micall.utils.new_atomic_file import new_atomic_text_file


def join_tables(inputs: Iterable[Path], column: str, output: Path) -> None:
    """
    Full‐outer‐join all CSVs in `inputs` on `column`, writing to `output`.
    Performs a Cartesian product for duplicate keys:
      if file A has N rows with key X and file B has M rows with key X,
      the output will have N*M rows for X.

    - Each CSV must have a header row.
    - `column` must be present in every CSV.
    - No two CSVs may share a non‐index column name (raises ValueError).
    - Missing values become empty strings.
    """

    inputs = tuple(inputs)
    if not inputs:
        raise ValueError("At least one input file must be provided.")

    # For each input file, build a map: index_value -> list of {non_index_col: cell}
    per_file_maps: List[Dict[str, List[Dict[str, str]]]] = []
    seen_cols: Set[str] = set()      # all non-index cols seen so far
    col_order: List[str] = []  # their global order

    for path in inputs:
        with path.open(newline='') as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError(f"No header found in {path}")
            if column not in reader.fieldnames:
                raise ValueError(f"Index column {column!r} not in {path}")

            # this file's non-index columns
            this_cols = [c for c in reader.fieldnames if c != column]
            overlap = seen_cols.intersection(this_cols)
            if overlap:
                raise ValueError(
                    "Column name overlap across inputs (excluding index column): "
                    + ", ".join(sorted(overlap))
                )
            for c in this_cols:
                seen_cols.add(c)
                col_order.append(c)

            # map key -> list of row‐dicts
            m: Dict[str, List[Dict[str, str]]] = {}
            for row in reader:
                key = row[column]
                entry = {c: row[c] for c in this_cols}
                m.setdefault(key, []).append(entry)
            per_file_maps.append(m)

    # collect all unique keys in first-seen order
    all_index_values = []
    seen_idx = set()
    for m in per_file_maps:
        for key in m:
            if key not in seen_idx:
                seen_idx.add(key)
                all_index_values.append(key)

    # final header
    final_cols = [column] + col_order

    with new_atomic_text_file(output) as out_f:
        writer = csv.DictWriter(out_f, fieldnames=final_cols)
        writer.writeheader()

        for key in all_index_values:
            # gather for each file: list of dicts (or [{}] if missing)
            lists_of_dicts = []
            for m in per_file_maps:
                rows = m.get(key)
                if rows:
                    lists_of_dicts.append(rows)
                else:
                    # no rows for this file → one "empty" row so product still yields one
                    lists_of_dicts.append([{}])

            # for every combination, merge into one output row
            for combo in itertools.product(*lists_of_dicts):
                row_out = {col: "" for col in final_cols}
                row_out[column] = key
                for part in combo:
                    row_out.update(part)
                writer.writerow(row_out)
