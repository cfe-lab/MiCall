import csv
import os
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

import requests

from micall.monitor import qai_helper


def parse_args():
    parser = ArgumentParser(description='Upload project definitions to QAI.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--qai_server',
        default=os.environ.get('MICALL_QAI_SERVER', 'http://localhost:3000'),
        help='server to post reviews on')
    parser.add_argument(
        '--qai_user',
        default=os.environ.get('MICALL_QAI_USER', 'bob'),
        help='user name for QAI server')
    parser.add_argument(
        '--qai_password',
        default=SUPPRESS,
        help='password for QAI server (default not shown)')
    parser.add_argument(
        '--pipeline_version',
        default='0-dev',
        help='version number')

    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def main():
    filename = 'HIV1_COM_2015_genome_DNA.csv'
    args = parse_args()

    if not os.path.exists(filename):
        form = {'ORGANISM': 'HIV',
                'ALIGN_TYPE': 'COM',
                'SUBORGANISM': 'HIV1',
                'PRE_USER': 'predefined',
                'REGION': 'GENOME',
                'START': '',
                'END': '',
                'GENO_SUB': 'All',
                'BASETYPE': 'DNA',
                'YEAR': '2015',
                'FORMAT': 'csv',
                'submit': 'Get Alignment'}
        response = requests.post("https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi",
                                 data=form)

        response.raise_for_status()
        # print(response.text)

        match = re.search(r'<pre>(.*)</pre>', response.text, re.DOTALL)
        with open(filename, 'w') as f:
            f.write(match.group(1))

    with qai_helper.Session() as session:
        session.login(args.qai_server,
                      args.qai_user,
                      args.qai_password)

        seed_groups = session.get_json("/lab_miseq_seed_groups")
        seed_group_name = 'HIV1-seed'
        for seed_group in seed_groups:
            if seed_group['name'] == seed_group_name:
                break
        else:
            raise RuntimeError('Seed group {} not found.'.format(seed_group_name))
        old_regions = session.get_json("/lab_miseq_regions", retries=0)
        hiv_seeds = {region['name']: region
                     for region in old_regions
                     if region['seed_group_id'] == seed_group['id']}
        del old_regions

        clean_count = 0
        dirty_count = 0
        recombinant_names = []
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            for description, seed_seq in reader:
                seed_seq = seed_seq.replace('-', '')
                name_fields = description.split('.')
                subtype, country = name_fields[:2]
                accession = name_fields[-1]
                if subtype[0].isdigit():
                    recombinant_names.append(description)
                    continue
                seed_name = '-'.join(('HIV1', subtype, country, accession, 'seed'))

                groups = re.findall(r'([^ACGT]+)', seed_seq)
                if groups:
                    dirty_count += 1
                    print('Unexpected bases found in {}: {}'.format(
                        seed_name,
                        ', '.join(groups)))
                else:
                    clean_count += 1
                old_region = hiv_seeds.pop(seed_name, None)
                if old_region:
                    old_seq = ''.join(old_region['reference'])
                    if old_seq != seed_seq:
                        print('expected: ' + seed_seq)
                        print('found:    ' + old_seq)
                        raise RuntimeError('Seed sequence {} does not match.'.format(
                            seed_name))
                elif len(seed_name) > 30:
                    print('Name too long: {!r}.'.format(seed_name))
                else:
                    session.post_json(
                        "/lab_miseq_regions",
                        {'name': seed_name,
                         'description': description,
                         'is_nucleotide': True,
                         'reference': seed_seq,
                         'seed_group_id': seed_group['id']})
        if recombinant_names:
            print(f'Skipped {len(recombinant_names)} recombinants:',
                  ', '.join(sorted(recombinant_names)))
        if hiv_seeds:
            seed_names = sorted(hiv_seeds.keys())
            should_delete = True
            print('Left over seeds:')
            if not should_delete:
                print(', '.join(seed_names))
            else:
                for seed_name in seed_names:
                    print(seed_name)
                    seed_id = hiv_seeds[seed_name]['id']
                    session.delete('{}/lab_miseq_regions/{}'.format(
                        args.qai_server,
                        seed_id))

        print('Done with {} clean and {} dirty.'.format(clean_count, dirty_count))

if __name__ == '__main__':
    main()
