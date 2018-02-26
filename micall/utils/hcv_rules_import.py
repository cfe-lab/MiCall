from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple
from functools import partial

from openpyxl import load_workbook

READY_TABS = ('NS3_GT1a', 'NS3_GT1b')
PHENOTYPE_SCORES = {'likely susceptible': 0,
                    'resistance possible': 4,
                    'resistance likely': 8,
                    'effect unknown': 'effect unknown'}
RuleSet = namedtuple('RuleSet', 'drug_name phenotype_column mutations')


def parse_args():
    parser = ArgumentParser(description='Read HCV rules from a spreadsheet.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('spreadsheet',
                        type=partial(load_workbook, read_only=False))
    return parser.parse_args()


def main():
    args = parse_args()
    for ws in args.spreadsheet:
        if ws.title not in READY_TABS:
            continue
        print(ws.title)
        if not ws.title.startswith('NS'):
            print('Skipped.')
            continue
        header_rows = []
        rule_sets = None
        for row in ws.rows:
            if rule_sets is None:
                if row[0].value and row[0].value.startswith('WT'):
                    rule_sets = find_rule_sets(ws, header_rows)
                    find_phenotype_scores(row, rule_sets)
                else:
                    header_rows.append(row)
            elif row[0].value is not None:
                find_phenotype_scores(row, rule_sets)
        if not rule_sets:
            print('No rules found!')
            continue
        for rule_set in rule_sets:
            print(rule_set.drug_name, rule_set.mutations)


def find_phenotype_scores(row, rule_sets):
    mutation = row[0].value
    for rule_set in rule_sets:
        phenotype_cell = row[rule_set.phenotype_column]
        phenotype = phenotype_cell.value
        if phenotype is not None:
            try:
                score = PHENOTYPE_SCORES[phenotype.strip(' ?*').lower()]
            except Exception:
                raise ValueError('Unknown phenotype at {}: {!r}.'.format(
                    phenotype_cell.coordinate,
                    phenotype))
            if score:
                rule_set.mutations[mutation] = score


def find_rule_sets(ws, header_rows):
    drug_row_num = len(header_rows) - 1
    drug_ranges = []
    for cell_range in ws.merged_cells.ranges:
        if cell_range.max_row == drug_row_num:
            drug_ranges.append(cell_range)
    rule_sets = []
    drug_row = header_rows[-2]
    for cell in drug_row:
        if cell.value:
            # drug_rules = RuleSet(cell.value, {})
            for drug_range in drug_ranges:
                if cell.coordinate in drug_range:
                    for col in range(drug_range.min_col,
                                     drug_range.max_col + 1):
                        header = header_rows[-1][col].value
                        if header == 'Phenotype':
                            rule_sets.append(RuleSet(cell.value, col, {}))
                            break
                    else:
                        raise ValueError(
                            'No Phenotype column between columns {} and {}.'.format(
                                drug_range.min_col,
                                drug_range.max_col))
                    break
    return rule_sets


main()
