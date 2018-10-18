#! /usr/bin/env python3.6
import os

# typing is only supported from python3.5 onwards
# from typing import List, Dict, Tuple

from argparse import ArgumentParser, FileType
import csv
from collections import defaultdict, Counter, namedtuple

import yaml

import micall.resistance.pdfreport as pdfreport

REPORT_CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'genreport.yaml')


def parse_args():
    parser = ArgumentParser(
        description="""Generate a resistance report in PDF form from a \
resistance call file and a mutations call file, both in csv format""")
    parser.add_argument(
        'resistance_csv', type=FileType('r'), help='resistance call file')
    parser.add_argument(
        'mutations_csv', type=FileType('r'), help='mutation call file')
    parser.add_argument(
        'res_report_pdf',
        type=FileType('wb'),
        help='A resistance report in PDF format.')
    parser.add_argument(
        "-s", "--sample", type=str, help='An optional sample name')
    return parser.parse_args()


def read_config(git_version):
    """Read in a configuration file for generating reports."""
    with open(REPORT_CONFIG_PATH) as fi:
        virus_configs = yaml.safe_load(fi)
    report_templates = []
    for virus_config in virus_configs:
        virus_config['generated_by_text'] = (
            virus_config['generated_by_text'].format(git_version or ''))
        report_templates.append(ReportTemplate(virus_config))
    return report_templates


ReportPage = namedtuple('ReportPage', 'resistance_calls mutations')


class ReportTemplate:
    def __init__(self, virus_config, raise_missing=False):
        self.virus_config = virus_config
        self.genotype_pages = defaultdict(lambda: ReportPage({}, {}))
        err_string = "Error in configuration file"
        if not isinstance(virus_config, dict):
            raise RuntimeError("""Configuration in {} must be a
    single dict class, but found a {}""".format(REPORT_CONFIG_PATH, type(virus_config)))
        # The following keys must be present in the dict
        known_keys = frozenset([
            'known_drugs', 'known_drug_classes', 'known_regions',
            'resistance_level_colours', 'disclaimer_text', "generated_by_text",
            "report_title", 'failure_message'
        ])
        present_keys = set(virus_config.keys())
        missing_keys = known_keys - present_keys
        unknown_keys = present_keys - known_keys
        if unknown_keys:
            message = "Unknown configuration: {}.".format(', '.join(
                sorted(unknown_keys)))
            raise ValueError(message)
        if missing_keys and raise_missing:
            message = "Missing configuration: {}.".format(', '.join(
                sorted(missing_keys)))
            raise ValueError(message)
        # --- do some sanity checks of the configuration file ---
        # check resistance_level_colours
        fld_name = 'resistance_level_colours'
        coltab = virus_config.get(fld_name, {})
        # convert this information into a dict for later.
        virus_config["res_level_dct"] = res_level_dct = {}
        for level, lev_tup in coltab.items():
            if len(lev_tup) != 3:
                raise RuntimeError("{}: {} must have 3 entries: {}".format(
                    err_string, fld_name, lev_tup))
            lev_str, bg_col, fg_col = lev_tup
            types_ok = isinstance(lev_str, str) and isinstance(
                bg_col, int) and isinstance(fg_col, int)
            if not types_ok:
                raise RuntimeError("{}: {} string, int, int expected '{}'".format(
                    err_string, fld_name, lev_tup))
            res_level_dct[level] = lev_str.upper()
        # check known_regions and turn it into a set of strings
        fld_name = 'known_regions'
        known_regions = virus_config.get(fld_name, [])
        if not isinstance(known_regions, list) or\
           sum([isinstance(s, str) for s in known_regions]) != len(known_regions):
            raise RuntimeError(
                "{}: {} must be a list of strings".format(err_string, fld_name))
        virus_config[fld_name] = frozenset(known_regions)
        # check known_drug_classes and turn it into a set of strings
        fld_name = 'known_drug_classes'
        known_drug_classes = virus_config.get(fld_name, [])
        dc_names = [drug_class[0] for drug_class in known_drug_classes]
        if not isinstance(known_drug_classes, list) or\
           sum([isinstance(s, str) for s in dc_names]) != len(known_drug_classes):
            raise RuntimeError(
                "{}: {} must be a list of strings".format(err_string, fld_name))
        virus_config['drug_class_tableheaders'] = dict(known_drug_classes)
        known_drug_classes = virus_config[fld_name] = frozenset(dc_names)
        virus_config["known_dclass_list"] = dc_names

        # check known_drugs
        fld_name = "known_drugs"
        drug_dct = virus_config.get(fld_name, {})
        if not isinstance(drug_dct, dict) or\
           sum([isinstance(s, str) for s in drug_dct.keys()]) != len(drug_dct):
            raise RuntimeError("{}: {} must be a dict with strings as keys".format(
                err_string, fld_name))
        # make sure each drug class is defined
        got_drug_classes = set(drug_dct.keys())
        if known_drug_classes != got_drug_classes:
            raise RuntimeError("{}: {} inconsistent drug_classes".format(
                err_string, fld_name))
        # Drug names and codes must be unique.
        drug_code_counts = Counter(drug[0]
                                   for drug_class in drug_dct.values()
                                   for drug in drug_class)
        drug_name_counts = Counter(drug[1]
                                   for drug_class in drug_dct.values()
                                   for drug in drug_class)
        duplicate_drug_codes = ', '.join(sorted(
            drug_code
            for drug_code, count in drug_code_counts.items()
            if count > 1))
        duplicate_drug_names = ', '.join(sorted(
            drug_name
            for drug_name, count in drug_name_counts.items()
            if count > 1))
        if duplicate_drug_codes:
            raise RuntimeError("{}: {} duplicate drug identifiers: {}.".format(
                err_string,
                fld_name,
                duplicate_drug_codes))
        if duplicate_drug_names:
            raise RuntimeError("{}: {} duplicate drug names: {}.".format(
                err_string,
                fld_name,
                duplicate_drug_names))
        # now generate a helper dict: drug : (drug_class, drugname)
        virus_config["drug_dct"] = dct = {}
        for d_class, ll in drug_dct.items():
            for d, d_name in ll:
                dct[d] = (d_class, d_name)

    def register_regions(self, regions):
        """ Register this report page for all of its known regions. """
        for region in self.virus_config.get('known_regions', []):
            regions[region] = self

    def register_drug_classes(self, drug_classes):
        for drug_class in self.virus_config.get('known_drug_classes', []):
            drug_classes[drug_class] = self

    def __repr__(self):
        report_title = self.virus_config.get('report_title')
        return "ReportPage({{'report_title': {!r}}})".format(report_title)

    def get_reported_genotypes(self):
        return sorted(self.genotype_pages.keys())

    def get_reported_drug_classes(self, genotype):
        reported_drug_codes = set(self.genotype_pages[genotype].resistance_calls.keys())
        return {class_code
                for class_code, drugs in self.virus_config['known_drugs'].items()
                if any(drug_code in reported_drug_codes
                       for drug_code, _ in drugs)}


def read_mutations(drug_classes, csv_file):
    """Read in a mutations file from CSV.

    Returns a list of dictionaries.
    """
    data_lst = list(csv.DictReader(csv_file, restkey="dummy"))
    tmp_dct = defaultdict(list)
    for od_num, od in enumerate(data_lst):
        d_class, mut_str = od['drug_class'], od["mutation"]
        genotype = od['genotype']
        tmp_dct[(genotype, d_class)].append('{}({:.0f}%)'.format(
            mut_str,
            100*float(od['prevalence'])))
    for (genotype, drug_class), mutations in tmp_dct.items():
        report_template = drug_classes[drug_class]
        report_page = report_template.genotype_pages[genotype]
        mutation_display = ', '.join(mutations)
        report_page.mutations[drug_class] = mutation_display


def read_resistance(regions, csv_file):
    """Read in a resistance call file from CSV.

    :param regions: {region: ReportTemplate} each ReportPage will receive the
        resistance calls from its regions
    :param csv_file: the resistance calls to load
    """
    data_lst = list(csv.DictReader(csv_file, restkey="dummy"))
    for od in data_lst:
        template = regions[od['region']]
        level = int(od['level'])
        drug_id = od['drug']
        genotype = od['genotype']
        report_page = template.genotype_pages[genotype]
        report_page.resistance_calls[drug_id] = (level, od["level_name"])


def gen_report(resistance_csv,
               mutations_csv,
               res_report_pdf,
               sample_name=None,
               git_version=None):
    """ Generate a PDF report file.

    :param resistance_csv: open CSV file with resistance calls
    :param mutations_csv: open CSV file with mutation calls
    :param res_report_pdf: PDF file name to write to
    :param sample_name: name to describe the sample on the report
    :param git_version: source code version to display
    """
    report_templates = read_config(git_version)
    regions = {}
    drug_classes = {}
    for report_template in report_templates:
        report_template.register_regions(regions)
        report_template.register_drug_classes(drug_classes)
    read_resistance(regions, resistance_csv)
    read_mutations(drug_classes, mutations_csv)

    pdfreport.write_report_one_column(report_templates, res_report_pdf, sample_name)


def main():
    args = parse_args()
    gen_report(args.resistance_csv, args.mutations_csv, args.res_report_pdf,
               args.sample)


if __name__ == '__main__':
    main()
