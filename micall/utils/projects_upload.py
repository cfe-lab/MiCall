import json
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from operator import itemgetter

from micall.core.project_config import ProjectConfig
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
    parser.add_argument(
        '--update_sequences',
        action='store_true',
        help="Update region sequences if they don't match.")

    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def main():
    args = parse_args()
    project_config = ProjectConfig.loadDefault()
    with open('../project_scoring.json', 'rU') as scoring_file:
        scoring_config = json.load(scoring_file)
    with qai_helper.Session() as session:
        session.login(args.qai_server,
                      args.qai_user,
                      args.qai_password)

        pipelines = session.get_json(
            "/lab_miseq_pipelines?version=" + args.pipeline_version,
            retries=0)
        if pipelines:
            raise RuntimeError('Pipeline {} already exists.'.format(
                args.pipeline_version))

        seed_groups = session.get_json("/lab_miseq_seed_groups")
        seed_group_ids = dict(map(itemgetter('name', 'id'), seed_groups))
        old_regions = session.get_json("/lab_miseq_regions", retries=0)
        regions = dict(((region['name'], region) for region in old_regions))
        for region_name, region_data in project_config.config['regions'].items():
            ref_seq = ''.join(region_data['reference'])
            region = regions.get(region_name)
            if region is None:
                seed_group_name = region_data['seed_group']
                seed_group_id = seed_group_ids.get(seed_group_name)
                if seed_group_id is None and seed_group_name:
                    seed_group = session.post_json("/lab_miseq_seed_groups",
                                                   {'name': seed_group_name})
                    seed_group_id = seed_group['id']
                    seed_group_ids[seed_group_name] = seed_group_id
                region = session.post_json(
                    "/lab_miseq_regions",
                    {'name': region_name,
                     'is_nucleotide': region_data['is_nucleotide'],
                     'reference': ref_seq,
                     'seed_group_id': seed_group_id})
                regions[region_name] = region
            elif region['reference'] != ref_seq:
                print("Reference doesn't match:", region_name)
                if args.update_sequences:
                    region['reference'] = ref_seq
                    session.post_json(f"/lab_miseq_regions/{region['id']}",
                                      region)

        pipeline = session.post_json("/lab_miseq_pipelines",
                                     {'version': args.pipeline_version})
        pipeline_id = pipeline['id']

        old_projects = session.get_json("/lab_miseq_projects", retries=0)
        projects = dict(((project['name'], project) for project in old_projects))
        for project_name, project_data in project_config.config['projects'].items():
            project = projects.get(project_name)
            if project is None:
                project = session.post_json(
                    "/lab_miseq_projects",
                    {'name': project_name,
                     'max_variants': project_data['max_variants']})
            project_version = session.post_json("/lab_miseq_project_versions",
                                                {'pipeline_id': pipeline_id,
                                                 'project_id': project['id']})
            for i, region_data in enumerate(project_data['regions']):
                scoring_data = scoring_config['projects'][project_name]['regions'][i]
                coordinate_region = regions[region_data['coordinate_region']]
                seed_region = regions[region_data['seed_region_names'][0]]
                seed_group_id = seed_region['seed_group_id']
                project_region = session.post_json(
                    "/lab_miseq_project_regions",
                    {'project_version_id': project_version['id'],
                     'coordinate_region_id': coordinate_region['id'],
                     'min_coverage1': scoring_data['min_coverage1'],
                     'min_coverage2': scoring_data['min_coverage2'],
                     'min_coverage3': scoring_data['min_coverage3'],
                     'seed_group_id': seed_group_id})

                for key_position in scoring_data['key_positions']:
                    session.post_json("/lab_miseq_key_positions",
                                      {'project_region_id': project_region['id'],
                                       'start_pos': key_position['start_pos'],
                                       'end_pos': key_position['end_pos']})

    print("Done.")


if __name__ == "__main__":
    main()
