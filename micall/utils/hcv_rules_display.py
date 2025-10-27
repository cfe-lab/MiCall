""" This script converts the HCV rules into a markdown file.

It is part of the project web site.
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

import os
from itertools import groupby
from operator import itemgetter

import yaml
from collections import defaultdict
from pyvdrm.hcvr import HCVR, AsiScoreCond, ScoreList, ScoreExpr, BoolTrue

from micall.resistance.asi_algorithm import HCV_RULES_VERSION, HCV_RULES_DATE


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
                                        f'hcv_rules{HCV_RULES_VERSION}.md')
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
    header = f"""\
---
title: CFE HCV Algorithm
description: version cfe-hcv {HCV_RULES_VERSION} - {HCV_RULES_DATE}
---

Rules are applied for each drug and each genotype found in a sample. Each rule
contains a score and a list of the mutation descriptions that receive that
score. A mutation is represented by a reference amino acid, the position of the
amino acid within the gene region, and the variant amino acid detected. An
exclamation mark in the mutation means any variant except those listed. "TRUE"
in the mutations column means that the value in the score column is always used
for that genotype, regardless of any mutations. A score is either a number, or
a name in quotes. The program adds up all the numerical scores and collects all
the names, before deciding on a resistance report.

* If the score is 4, then it reports "Resistance Possible".
* If the score is 8 or more, then it reports "Resistance Likely".
* If the score is 0 with the "Not available" name, then it reports "Resistance
    Interpretation Not Available", which means that resistance interpretation
    is not available for the specified genotype and drug.
* If the score is 0 with the "Not indicated" name, then it reports
    "Not Indicated", which means that the drug is not indicated for use in
    Canada.
* If the score is 0 with the "Effect unknown" name, then it reports "Mutations
    Detected; Effect Unknown".
* If there is insufficient coverage due to poor amplification, sequencing, or
    mapping, then it reports "Sequence does not meet quality-control standards".
* Otherwise, it reports "Likely Susceptible".

Here are the rules used for each of the drugs. Below that is a description of
how the scores are chosen.

| Gene | Drug | Genotype | Score | Mutations |
|------|------|----------|-------|-----------|
"""
    footer = """\

The overall resistance scores of a mutant strain are derived from drug
susceptibility and clinical observations.

The drug susceptibility of a mutant is compared with the drug susceptibility of
the wild type, and the fold change is reported. For example, a mutant with five
times lower drug susceptibility than the wild type reports a fold change of 5.

Some mutants have a single amino acid changed from the wild type, and some
mutants have more than one change.

Another type of study reports mutations that are observed in clinical patients
with virological failure.

Based on all these study results, we assign the scores shown in the table
above. For each drug and genotype combination, we choose a middle range of fold
change values for drug susceptibility. Any mutations that fall in that middle
range get a score of 4 (resistance possible), mutations above that range get a
score of 8 (resistance likely), and mutations below that range get no score
(susceptible). Mutations observed in clinical patients with virological failure
get a score of 4 (resistance possible).

If reports for mutants with more than one change disagree with the reports for
the individual changes, then the combination takes precedence.

As an example, here are the fold change ranges for NS3 drugs, genotype 1a:

| Drugs | Middle Range |
|-------|--------------|
|Grazoprevir, Paritaprevir, Glecaprevir, Simeprevir, Asunaprevir| 5 - 10 |
|Voxilaprevir| 2.6 - 10* |

`*` Some studies only reported ranges of fold changes, so we used the 2.5 - 20
range as a middle range for those studies.
"""
    report_file.write(header)
    name_width = max(len(rule['name']) for rule in rules)
    displays = ('8',
                '4',
                'Not available',
                'Not indicated',
                'Effect unknown')
    for region, drugs in groupby(rules, itemgetter('region')):
        for drug_num, drug in enumerate(drugs):
            drug_name = drug['name']
            for genotype_num, genotype in enumerate(drug['genotypes']):
                rule_text = genotype['rules']
                genotype_text = genotype['genotype']
                mutations = parse_rule(rule_text)
                for display in displays:
                    display_mutations = mutations[display]
                    if display_mutations:
                        report_file.write(
                            '| {:4} | {:{}} | {:2} | {} | {} |\n'.format(
                                region,
                                drug_name,
                                name_width,
                                genotype_text,
                                display,
                                ', '.join(display_mutations)))
                        region = ''
                        drug_name = ''
                        genotype_text = ''
    report_file.write(footer)


def parse_rule(rule_text):
    rule = HCVR(rule_text)
    assert isinstance(rule.dtree, AsiScoreCond), repr(rule.dtree)
    score_list = rule.dtree.children[0]
    assert isinstance(score_list, ScoreList), repr(score_list)
    mutations = defaultdict(list)
    for score_expr in score_list.children:
        assert isinstance(score_expr, ScoreExpr)
        children = score_expr.children
        if len(children) == 4:
            # Flag
            mutation, _, flag, _ = children
            display = flag
        elif len(children) == 2:
            # Numeric score
            mutation, score = children
            display = score
        else:
            raise ValueError('Unexpected score expression: {!r}'.format(
                children))
        if isinstance(mutation, BoolTrue):
            mutation = 'TRUE'
        else:
            mutation = str(mutation.mutations)
        mutations[display].append(mutation)
    return mutations


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

if __name__ == '__main__':
    main()
