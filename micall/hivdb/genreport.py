#! /usr/bin/env python3.4
import os

# typing is only supported from python3.5 onwards
# from typing import List, Dict, Tuple

from argparse import ArgumentParser, FileType
import csv
import yaml

try:
    import pdfreport
except ImportError:
    import micall.hivdb.pdfreport as pdfreport

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


def read_config():
    """Read in a configuration file for generating reports."""
    cfg_name = REPORT_CONFIG_PATH
    with open(cfg_name, "r") as fi:
        try:
            cfg_dct = yaml.safe_load(fi)
        except yaml.scanner.ScannerError as e:
            raise RuntimeError(
                "Scanning Error reading FSM from '{}'\n{}".format(cfg_name, e))
        except yaml.parser.ParserError as e:
            raise RuntimeError(
                "Parse Error reading FSM from '{}'\n{}".format(cfg_name, e))
        except:
            print("failed to read config file {}".format(cfg_name))
            raise
    err_string = "Error in configuration file"
    if not isinstance(cfg_dct, dict):
        raise RuntimeError("""Configuration in {} must be a
single dict class, but found a {}""".format(REPORT_CONFIG_PATH, type(cfg_dct)))
    # The following keys must be present in the dict
    KNOWN_KEYS = frozenset([
        'known_drugs', 'known_drug_classes', 'known_regions',
        'resistance_level_colours', 'disclaimer_text', "generated_by_text",
        "report_title"
    ])
    PRESENT_KEYS = set(cfg_dct.keys())
    missing_keys = KNOWN_KEYS - PRESENT_KEYS
    unknown_keys = PRESENT_KEYS - KNOWN_KEYS
    if missing_keys != set() or unknown_keys != set():
        print("Missing keys: {}".format(", ".join(missing_keys)))
        print("Unknown keys: {}".format(", ".join(unknown_keys)))
        raise RuntimeError(err_string)
    # --- do some sanity checks of the configuration file ---
    # check resistance_level_colours
    # there must be exactly six resistance levels colours [0--5]
    fld_name = 'resistance_level_colours'
    coltab = cfg_dct[fld_name]
    if list(range(0, 6)) != list(coltab.keys()):
        raise RuntimeError(
            "{}: {} must have levels [0..5]".format(err_string, fld_name))
    # convert this information into a dict for later.
    cfg_dct["res_level_dct"] = res_level_dct = {}
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
    known_regions = cfg_dct[fld_name]
    if not isinstance(known_regions, list) or\
       sum([isinstance(s, str) for s in known_regions]) != len(known_regions):
        raise RuntimeError(
            "{}: {} must be a list of strings".format(err_string, fld_name))
    cfg_dct[fld_name] = frozenset(known_regions)
    # check known_drug_classes and turn it into a set of strings
    fld_name = 'known_drug_classes'
    known_drug_classes = cfg_dct[fld_name]
    dc_names, dc_tableh = zip(*known_drug_classes)
    if not isinstance(known_drug_classes, list) or\
       sum([isinstance(s, str) for s in dc_names]) != len(known_drug_classes):
        raise RuntimeError(
            "{}: {} must be a list of strings".format(err_string, fld_name))
    cfg_dct['drug_class_tableheaders'] = dict(known_drug_classes)
    known_drug_classes = cfg_dct[fld_name] = frozenset(dc_names)
    cfg_dct["known_dclass_list"] = dc_names

    # check known_drugs
    fld_name = "known_drugs"
    drug_dct = cfg_dct[fld_name]
    if not isinstance(drug_dct, dict) or\
       sum([isinstance(s, str) for s in drug_dct.keys()]) != len(drug_dct):
        raise RuntimeError("{}: {} must be a dict with strings as keys".format(
            err_string, fld_name))
    # make sure each drug class is defined
    got_drug_classes = set(drug_dct.keys())
    if known_drug_classes != got_drug_classes:
        raise RuntimeError("{}: {} inconsistent drug_classes".format(
            err_string, fld_name))
    # drug_names must be unique
    numdrugs, drug_set, drugname_set = 0, set(), set()
    for ll in drug_dct.values():
        numdrugs += len(ll)
        dd_set, dn_set = [set(l) for l in zip(*ll)]
        drug_set |= dd_set
        drugname_set |= dn_set
    if len(drug_set) != numdrugs:
        raise RuntimeError("{}: {} drug identifiers are not unique!".format(err_string, fld_name))
    if len(drugname_set) != numdrugs:
        raise RuntimeError("{}: {} drug_names are not unique!".format(err_string, fld_name))
    # now generate a helper dict: drug : (drug_class, drugname)
    cfg_dct["drug_dct"] = dct = {}
    for d_class, ll in drug_dct.items():
        for d, d_name in ll:
            dct[d] = (d_class, d_name)
    return cfg_dct


def read_mutations(cfg_dct, csv_file):
    """Read in a resistence call file from CSV.
    Returns a list of dictionaries.
    """
    err_string = "Error in mutations file '{}'".format(csv_file.name)
    exp_set = frozenset("drug_class,mutation".split(","))
    try:
        data_lst = list(csv.DictReader(csv_file, restkey="dummy"))
    except csv.Error:
        print("{}: Failed to read csv file".format(err_string))
        raise
    # make sure that all lines have exactly the required fields
    if sum([set(od.keys()) == exp_set for od in data_lst]) != len(data_lst):
        raise RuntimeError("{}: unexpected data found.".format(err_string))
    # simple sanity check of the fields
    # generate a dict of drug_class -> string of mutations
    # also set the strings to 'NONE' if appropriate
    tmp_dct = {}
    known_drug_class_set = cfg_dct['known_drug_classes']
    glob_err = False
    for od_num, od in enumerate(data_lst):
        line_err = False
        d_class, mut_str = od['drug_class'], od["mutation"]
        if d_class not in known_drug_class_set:
            print("unknown drug_class {}".format(d_class))
            line_err = True
        if line_err:
            print("{}: dataset {}: error with {}\n\n".format(
                err_string, od_num + 1, od))
            glob_err = True
        tmp_dct.setdefault(d_class, []).append(mut_str)
    if glob_err:
        raise RuntimeError(
            "{}: fatal error(s) detected: giving up...".format(err_string))
    cfg_dct["mutation_strings"] = mut_dct = {}
    for d_class in known_drug_class_set:
        mut_str = ", ".join(tmp_dct[d_class]) if d_class in tmp_dct else "None"
        mut_dct[d_class] = "Relevant {} Mutations: ".format(d_class) + mut_str
    return data_lst


def read_resistance(cfg_dct, csv_file):
    """Read in a resistence call file from CSV.
    Returns a list of dictionaries.
    """
    err_string = "Error in resistance file '{}'".format(csv_file.name)
    exp_set = frozenset(
        "region,drug_class,drug,drug_name,level,level_name,score".split(","))
    try:
        data_lst = list(csv.DictReader(csv_file, restkey="dummy"))
    except csv.Error:
        print("{}: Failed to read csv file".format(err_string))
        raise
    # print("GOO", res_data_lst)
    # make sure that all lines have exactly the required fields
    if sum([set(od.keys()) == exp_set for od in data_lst]) != len(data_lst):
        raise RuntimeError("{}: unexpected data found.".format(err_string))
    # simple sanity check of the fields
    # also create a dictionary key: drug_id (each drug must only be reported once)
    cfg_dct["res_results"] = resdct = {}
    known_regions_set = cfg_dct['known_regions']
    known_drug_class_set = cfg_dct['known_drug_classes']
    known_drugs_dct = cfg_dct['drug_dct']
    known_drugs_nameset = frozenset(known_drugs_dct.keys())
    res_level_dct = cfg_dct["res_level_dct"]
    known_res_levels = frozenset(res_level_dct.keys())
    known_res_level_names = frozenset(res_level_dct.values())
    glob_err = False
    for od_num, od in enumerate(data_lst):
        line_err = False
        if od['region'] not in known_regions_set:
            print("unknown region {}".format(od['region']))
            line_err = True
        if od['drug_class'] not in known_drug_class_set:
            print("unknown drug_class {}".format(od['drug_class']))
            line_err = True
        drug_id = od['drug']
        if drug_id not in known_drugs_nameset:
            print("unknown drug {}".format(drug_id))
            line_err = True
        # NOTE Until now, all elements are strings; convert level into int and score into float
        level_str, level_name, score_str = od['level'], od['level_name'], od['score']
        try:
            level = int(level_str)
        except ValueError:
            print("level is not an int {}".format(level_str))
            line_err = True
        if drug_id in resdct:
            print("multiple resistance results for drug {}".format(drug_id))
            line_err = True
        else:
            resdct[drug_id] = (level, od["level_name"])
        try:
            float(score_str)
        except ValueError:
            print("score is not a float {}".format(score_str))
            line_err = True
        if level not in known_res_levels:
            print("unknown resistance level '{}'".format(level))
            line_err = True
        if level_name.upper() not in known_res_level_names:
            print("unknown resistance level_name '{}'".format(level_name))
            line_err = True
        # TODO: check for score values here --
        if line_err:
            print("{}: dataset {}: error with {}\n\n".format(
                err_string, od_num + 1, od))
            glob_err = True
    if glob_err:
        raise RuntimeError(
            "{}: fatal error(s) detected: giving up...".format(err_string))
    return data_lst


def gen_report(resistance_csv, mutations_csv, res_report_pdf,
               sample_name=None):
    """Read in resistance_csv, validate its input and generate
    a PDF report to res_report+pdf"""
    cfg_dct = read_config()
    res_lst = read_resistance(cfg_dct, resistance_csv)
    mut_lst = read_mutations(cfg_dct, mutations_csv)

    pdfreport.write_report_one_column(cfg_dct, res_lst, mut_lst, res_report_pdf, sample_name)


def main():
    args = parse_args()
    gen_report(args.resistance_csv, args.mutations_csv, args.res_report_pdf,
               args.sample)


if __name__ == '__main__':
    main()
