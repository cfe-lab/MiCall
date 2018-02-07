""" This script converts the HCV rules into a markdown file.

It is part of the project web site.
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

import os
from itertools import groupby
from operator import itemgetter

import sys
import yaml


def parse_args():
    parser = ArgumentParser(description="Convert HCV rules to a markdown page.",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    micall_path = os.path.dirname(os.path.dirname(__file__))
    default_rules_path = os.path.join(micall_path,
                                      'resistance',
                                      'hcv_rules.yaml')
    default_config_path = os.path.join(micall_path,
                                       'resistance',
                                       'genreport.yaml')
    default_display_path = os.path.join(os.path.dirname(micall_path),
                                        'docs',
                                        'hcv_rules.md')
    parser.add_argument('-r',
                        '--rules',
                        type=FileType(),
                        default=default_rules_path,
                        help='file with HCV rules in YAML format')
    parser.add_argument('-c',
                        '--config',
                        type=FileType(),
                        default=default_config_path,
                        help='report configuration in YAML format')
    parser.add_argument('-d',
                        '--display',
                        type=FileType('w'),
                        default=default_display_path,
                        help='file to write markdown display into')
    return parser.parse_args()


def load_drug_names(config_file, report_prefix='Hepatitis C'):
    report_config = yaml.safe_load(config_file)
    for report in report_config:
        if report['report_title'].startswith(report_prefix):
            known_drugs = report['known_drugs']
            return {drug[0]: drug[1]
                    for drug_class in known_drugs.values()
                    for drug in drug_class}
    raise ValueError('No report title starts with {!r}.'.format(report_prefix))


def write_report(report_file, rules):
    header = """\
---
title: CFE HCV Algorithm
description: version cfe-hcv 1.5 - 07-Dec-2016
---

Rules are applied for each drug and each genotype found in a sample. Each rule
contains a mutation description, the `=>` symbol, and an action. A mutation
description is the position number of an amino acid within the gene region,
plus a variant amino acid. An exclamation mark in the mutation means any
variant except those listed. A "TRUE" mutation is always applied for that
genotype. An action is either a number that adds to the score, or a name in
quotes. The program adds up all the scores and collects all the names, before
deciding on a resistance report.

* If the score is 4, then it reports "Resistance Possible".
* If the score is 8 or more, then it reports "Resistance Likely".
* If the score is 0 with the "Not available" name, then it reports "Resistance
    Interpretation Not Available".
* If the score is 0 with the "Not indicated" name, then it reports
    "Not Indicated".
* If the score is 0 with the "Effect unknown" name, then it reports "Mutations
    Detected; Effect Unknown".
* If there are gaps in coverage for some of the mutation locations, then it
    reports "Sequence does not meet quality-control standards".
* Otherwise, it reports "Likely Susceptible".

Here are the rules used for each of the drugs.

| Gene | Drug | Genotype | Rules |
|------|------|----------|-------|
"""
    report_file.write(header)
    name_width = max(len(rule['name']) for rule in rules)
    for region, drugs in groupby(rules, itemgetter('region')):
        for drug_num, drug in enumerate(drugs):
            drug_name = drug['name']
            for genotype_num, genotype in enumerate(drug['genotypes']):
                rule_text = genotype['rules'].replace(' => ', '=>')
                report_file.write('| {:4} | {:{}} | {:2} | {} |\n'.format(
                    region,
                    drug_name,
                    name_width,
                    genotype['genotype'],
                    rule_text))
                region = ''
                drug_name = ''


def main():
    args = parse_args()
    drug_names = load_drug_names(args.config)
    rules = yaml.safe_load(args.rules)
    for rule in rules:
        drug_code = rule['code']
        rule['name'] = drug_names[drug_code]
        rule['region'] = rule['genotypes'][0]['region']
    rules.sort(key=itemgetter('region', 'name'))
    write_report(args.display, rules)


main()
