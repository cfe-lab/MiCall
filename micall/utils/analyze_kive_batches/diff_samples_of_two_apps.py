
from pathlib import Path
# import pandas as pd


def diff_samples_of_two_apps(input: Path, app1: str, app2: str, output: Path) -> None:
    # 1. Read the CSV
    # df = pd.read_csv(input)

    # TODO:
    #
    #   This function is based on sample names. Column "sample" is the "primary unique key".
    #
    #   Select samples that are in both `app1` and `app2`. I.e. rows where "app" == "app1" and "app" == "app2" for a single "sample".
    #   Let $LEFT be the value of a column for `app1`'s version of sample.
    #   Let $RIGHT be the value of a column for `app2`'s version of sample.
    #   For each numeric column, compute the difference, `$LEFT - $RIGHT`, and store in that column for the resulting table.
    #   For each non-numeric column, set the value to be `"$LEFT/$RIGHT"`. If `$LEFT == $RIGHT`, then just `$LEFT`, without slash comparison.
    #
    #   For example for `input` like
    #
    #      | sample | app   | size | type |
    #      | 1      | bob   | 50   | x    |
    #      | 1      | celia | 60   | y    |
    #      | 1      | alice | 70   | z    |
    #
    #   And `app1 == alice` and `app2 = bob`
    #
    #   The output should be
    #
    #      | sample | app         | size | type |
    #      | 1      | alice/bob   | 20   | x/z  |
    #

    # 3. Write out to CSV
    # result.to_csv(output, index=False)
    raise NotImplementedError()
