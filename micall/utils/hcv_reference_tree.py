#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import os
import re
from collections import Counter
from subprocess import run

import ete3
import requests

""" Retrieve HCV reference sequences for building a tree.

Download FastTree from http://www.microbesonline.org/fasttree/

Download FigTree for viewing the trees from http://tree.bio.ed.ac.uk/software/figtree/

If we need more sequences, this query looks useful:
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore\
&term=hcv[Organism]+AND+NS5A[Gene+Name]+AND+1200:10000[Sequence+Length]&usehistory=y

Then get the details with this:
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\
&rettype=gb&retmode=xml&retmax=10&query_key=1&WebEnv=???
"""

logger = logging.getLogger(__name__)
EXCLUDED_ACCESSIONS = ('HQ537007',  # 1a that doesn't tree with the other 7.
                       'EU246940')  # 6u that doesn't tree with the other 2.


def parse_args():
    parser = ArgumentParser(
        description='Fetch and filter HCV references for all subtypes.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder',
                        default='hcv_reference_tree',
                        help='folder to download references to')

    return parser.parse_args()


def fetch_raw_hcv(folder):
    raw_hcv_path = os.path.join(folder, 'raw_hcv.fasta')
    if os.path.exists(raw_hcv_path):
        return raw_hcv_path
    logger.info('Fetching raw HCV sequences.')
    response = requests.post(
        'https://hcv.lanl.gov/cgi-bin/NEWALIGN/align.cgi',
        data={'ORGANISM': 'HCV',
              'ALIGN_TYPE': 'REF',
              'SUBORGANISM': 'HCV',
              'PRE_USER': 'predefined',
              'REGION': 'GENOME',
              'START': '',
              'END': '',
              'GENO_SUB': 'ALL',
              'BASETYPE': 'DNA',
              'YEAR': '2014',
              'FORMAT': 'fasta',
              'submit': 'Get Alignment'})
    match = re.search(r'<a href="([^"]*)">Download', response.text)
    assert match
    response = requests.get('https://www.hcv.lanl.gov' + match.group(1))
    with open(raw_hcv_path, 'w') as raw_hcv:
        raw_hcv.write(response.text)
    return raw_hcv_path


def filter_hcv(raw_hcv_path, excluded=tuple()):
    filtered_hcv_path = os.path.join(os.path.dirname(raw_hcv_path),
                                     'filtered_hcv.fasta')
    if os.path.exists(filtered_hcv_path):
        return filtered_hcv_path

    logger.info('Filtering HCV.')
    with open(raw_hcv_path) as raw_hcv, open(filtered_hcv_path, 'w') as filtered_hcv:
        invalid_headers = filter_hcv_fasta(raw_hcv, filtered_hcv, excluded)

    if invalid_headers:
        message = ', '.join(sorted(invalid_headers))
        logger.warning('Invalid headers (%d) in FASTA: %s.',
                       len(invalid_headers),
                       message)
    return filtered_hcv_path


def get_subtype(header):
    parts = header.lower().split('.')
    assert parts[0] == '>ref'
    match = re.match(r'(con_)?(\d+[a-z])(\(\d+\))?$', parts[1])
    return match and match.group(2)


def filter_hcv_fasta(raw_hcv, filtered_hcv, excluded=tuple()):
    seen_headers = set()
    invalid_headers = []
    is_skipping = False
    for line in raw_hcv:
        if line.startswith('>'):
            line = line.strip()
            is_skipping = line in seen_headers
            if is_skipping:
                invalid_headers.append(line + ' (duplicate)')
                continue
            for exclusion in excluded:
                if exclusion in line:
                    is_skipping = True
            if is_skipping:
                invalid_headers.append(line + ' (excluded)')
                continue
            seen_headers.add(line)
            subtype = get_subtype(line)
            is_skipping = subtype is None
            if is_skipping:
                invalid_headers.append(line)
            else:
                line = line.replace('(', '[').replace(')', ']')
                line = f'{line}-{subtype}\n'
        if not is_skipping:
            filtered_hcv.write(line)
    return invalid_headers


def build_tree(filtered_hcv_path):
    tree_path = os.path.join(os.path.dirname(filtered_hcv_path),
                             'hcv_references.tree')
    if os.path.exists(tree_path):
        return tree_path

    logger.info('Building tree.')
    with open(tree_path, 'wb') as tree_file:
        run(['FastTree', '-nt', filtered_hcv_path],
            check=True,
            stdout=tree_file)
    return tree_path


def check_tree(tree_source, report_file=None):
    tree = ete3.Tree(tree_source)
    subtypes = Counter()
    genotypes = Counter()
    for leaf in tree:
        subtype = leaf.name.split('-')[-1]
        subtypes[subtype] += 1
        genotype = subtype[:-1]
        genotypes[genotype] += 1
        leaf.add_features(subtype=subtype, genotype=genotype)
    should_print = False
    for subtype, count in sorted(subtypes.items()):
        if count < 2:
            continue
        is_monophyletic, status, nodes = tree.check_monophyly(
            values=[subtype],
            target_attr='subtype')
        if not is_monophyletic:
            node_names = sorted(node.name for node in nodes)
            print(subtype, status, ', '.join(node_names), file=report_file)
            should_print = True
    for genotype, count in sorted(genotypes.items()):
        if count < 2:
            continue
        is_monophyletic, status, nodes = tree.check_monophyly(
            values=[genotype],
            target_attr='genotype')
        if status == 'paraphyletic' and genotype == '6':
            # We let genotype 6 get away with paraphyletic.
            continue
        if not is_monophyletic:
            node_names = sorted(node.name for node in nodes)
            print(genotype, status, ', '.join(node_names), file=report_file)
            should_print = True

    if should_print:
        print(tree, file=report_file)


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s[%(levelname)s]%(module)s:%(lineno)d - %(message)s')
    logger.info('Starting.')

    try:
        os.makedirs(args.folder, exist_ok=True)
        raw_hcv_path = fetch_raw_hcv(args.folder)
        filtered_hcv_path = filter_hcv(raw_hcv_path, EXCLUDED_ACCESSIONS)
        tree_path = build_tree(filtered_hcv_path)
        check_tree(tree_path)
        logger.info('Done.')
    except Exception:
        logger.error('Failed.')
        raise


if __name__ == '__main__':
    main()
