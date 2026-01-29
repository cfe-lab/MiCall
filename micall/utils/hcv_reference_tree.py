#!/usr/bin/env python3
import shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import os
import re
from collections import Counter, defaultdict
from csv import DictReader
from itertools import chain
from subprocess import run

import ete3
from Levenshtein import distance
import requests
import sys

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
    parser.add_argument('consensus_files',
                        nargs='*',
                        help='files to read MAX consensus sequences from for testing')

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
              'submit': 'Get Alignment'},
        verify=False)
    match = re.search(r'<a href="([^"]*)">Download', response.text)
    assert match
    response = requests.get('https://www.hcv.lanl.gov' + match.group(1),
                            verify=False)
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


def build_tree(fasta_path, check_cache=False):
    root, ext = os.path.splitext(fasta_path)
    tree_path = root + '.tree'
    if check_cache and os.path.exists(tree_path):
        return tree_path

    logger.info('Building tree.')
    with open(tree_path, 'wb') as tree_file:
        run(['FastTree', '-quiet', '-nt', fasta_path],
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
        previous_ref = None
        sample_nodes = []
        for node in tree:
            is_reference = node.name.startswith('Ref.')
            if is_reference:
                if sample_nodes:
                    next_ref = node
                    for sample_node in sample_nodes:
                        next_dist = tree.get_distance(sample_node, next_ref)
                        neighbouring_subtypes = set()
                        if previous_ref is None:
                            prev_dist = None
                        else:
                            prev_dist = tree.get_distance(previous_ref, sample_node)
                            if prev_dist <= next_dist:
                                if sample_node.subtype == previous_ref.subtype:
                                    continue
                                neighbouring_subtypes.add(previous_ref.subtype)
                        if prev_dist is None or next_dist <= prev_dist:
                            if sample_node.subtype == next_ref.subtype:
                                continue
                            neighbouring_subtypes.add(next_ref.subtype)
                        neighbours = ', '.join(sorted(neighbouring_subtypes))
                        print(f'{sample_node.name} treed with {neighbours}', file=report_file)
                    sample_nodes.clear()
                previous_ref = node
            else:
                sample_nodes.append(node)
        print(tree, file=report_file)


def iterate_fasta(fasta):
    header = sequence = ''
    for line in chain(fasta, ['>']):
        if line.startswith('>'):
            if header:
                yield header, sequence
            header = line[1:].strip()
            sequence = ''
        else:
            sequence += line.strip()


def check_distances(combined_hcv, report_file=None):
    references = [(header, sequence.replace('-', ''))
                  for header, sequence in iterate_fasta(combined_hcv)
                  if header.startswith('Ref.')]

    combined_hcv.seek(0)
    samples = ((header, sequence.replace('-', ''))
               for header, sequence in iterate_fasta(combined_hcv)
               if header.startswith('Sample.'))
    for header, sequence in samples:
        reported_genotype = header.split('-')[-1]
        reported_ref = best_ref = min_distance = reported_distance = None
        reported_size = best_size = 0
        for ref_header, ref_seq in references:
            ref_genotype = ref_header.split('-')[-1]
            d = distance(sequence, ref_seq)
            if (ref_genotype == reported_genotype and
                    (reported_distance is None or d < reported_distance)):
                reported_distance = d
                reported_ref = ref_header
                reported_size = len(ref_seq)
            if min_distance is None or d < min_distance:
                min_distance = d
                best_ref = ref_header
                best_size = len(ref_seq)
        if min_distance != reported_distance:
            best_genotype = best_ref.split('-')[-1]
            print(f'Reported {reported_genotype}, but {best_genotype} is '
                  f'closer: {header}(0/{len(sequence)}), '
                  f'{reported_ref}({reported_distance}/{reported_size}), '
                  f'{best_ref}({min_distance}/{best_size}).',
                  file=report_file)


def combine_samples(filtered_hcv, consensus_file, coverage_scores, combined_file):
    """ Combine MAX consensus for each sample with HCV reference sequences. """
    reader = DictReader(coverage_scores)
    resistance_regions = ('NS3', 'NS5a', 'NS5b')
    covered_seeds = defaultdict(set)  # {(sample, seed): {region}}
    for row in reader:
        if row['on.score'] == '4':
            region = row['region'].split('-')[-1]
            if region in resistance_regions:
                region += '*'
            key = (row['sample'], row['seed'])
            covered_seeds[key].add(region)
    shutil.copyfileobj(filtered_hcv, combined_file)
    reader = DictReader(consensus_file)
    for row in reader:
        cutoff = row['consensus-percent-cutoff']
        if cutoff != 'MAX':
            continue
        sample_name = row['sample']
        seed = row['region']
        covered_regions = covered_seeds[(sample_name, seed)]
        if not covered_regions:
            continue
        region_text = '~'.join(sorted(covered_regions))
        seed_parts = seed.split('-')
        if seed_parts[0] != 'HCV':
            continue
        subtype = seed_parts[-1]
        combined_file.write(f'>Sample.{subtype}.{sample_name}~{region_text}-{subtype}\n')
        combined_file.write(row['sequence'])
        combined_file.write('\n')


def align_samples(combined_hcv_path, aligned_hcv_path):
    with open(aligned_hcv_path, 'wb') as aligned_hcv:
        run(['mafft', '--quiet', combined_hcv_path], stdout=aligned_hcv)


def check_sample_trees(filtered_hcv_path, consensus_files):
    working_path = os.path.dirname(filtered_hcv_path)
    combined_path = os.path.join(working_path, 'combined.fasta')
    aligned_path = os.path.join(working_path, 'aligned.fasta')
    with open(filtered_hcv_path) as filtered_hcv:
        for consensus_file_name in consensus_files:
            filtered_hcv.seek(0)
            coverage_path = os.path.join(os.path.dirname(consensus_file_name),
                                         'coverage_scores.csv')
            with open(consensus_file_name) as consensus_file, \
                    open(coverage_path) as coverage_file, \
                    open(combined_path, 'w') as combined_file:
                combine_samples(filtered_hcv, consensus_file, coverage_file, combined_file)
            logger.info('Checking %s', consensus_file_name)
            with open(combined_path) as combined_file:
                check_distances(combined_file)
            align_samples(combined_path, aligned_path)
            tree_path = build_tree(aligned_path)
            check_tree(tree_path)


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format='%(asctime)s[%(levelname)s]%(module)s:%(lineno)d - %(message)s')
    logger.info('Starting.')

    try:
        os.makedirs(args.folder, exist_ok=True)
        raw_hcv_path = fetch_raw_hcv(args.folder)
        filtered_hcv_path = filter_hcv(raw_hcv_path, EXCLUDED_ACCESSIONS)
        tree_path = build_tree(filtered_hcv_path, check_cache=True)
        check_tree(tree_path)
        check_sample_trees(filtered_hcv_path, args.consensus_files)
        logger.info('Done.')
    except Exception:
        logger.error('Failed.')
        raise


if __name__ == '__main__':
    main()
