import json
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from operator import itemgetter
from pathlib import Path
import logging
import sys
import time
import urllib3

from micall.core.project_config import ProjectConfig
from micall.utils.externals import ProjectsScoringFile
from micall.monitor import qai_helper

logger = logging.getLogger(__name__)

# Suppress urllib3 connection retry warnings
urllib3.connectionpool.log.setLevel(logging.ERROR)


def configure_logging(args):
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(
        level=logger.level,
        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s')


def find_missing_scoring_regions(project_config, scoring_config):
    logger.debug('Checking scoring coverage for amino-acid project regions.')
    scoring_projects = scoring_config.get('projects', {})
    missing = []
    for project_name, project_data in project_config.config['projects'].items():
        scoring_project = scoring_projects.get(project_name, {})
        scoring_regions = {
            region_data['coordinate_region']
            for region_data in scoring_project.get('regions', [])
        }
        logger.debug(
            'Project %s has %d scored regions in project_scoring.json.',
            project_name,
            len(scoring_regions))
        for region_data in project_data['regions']:
            region_name = region_data['coordinate_region']
            if project_config.config['regions'][region_name]['is_nucleotide']:
                continue
            if region_name not in scoring_regions:
                missing.append((project_name, region_name))
                logger.debug(
                    'Missing scoring mapping found: %s/%s',
                    project_name,
                    region_name)
    logger.debug('Scoring coverage check complete. Missing mappings: %d', len(missing))
    return missing


def fetch_region_by_name(session, region_name, max_attempts=5, wait_seconds=0.2):
    logger.debug('Fetching region by name from QAI: %s', region_name)
    for attempt in range(1, max_attempts + 1):
        for region in session.get_json('/lab_miseq_regions'):
            if region['name'] == region_name:
                logger.debug('Fetched region %s with id %s', region_name, region.get('id'))
                return region

        if attempt < max_attempts:
            delay = wait_seconds * (2 ** (attempt - 1))
            logger.debug(
                'Region %s not visible yet (attempt %d/%d). Retrying in %.2fs.',
                region_name,
                attempt,
                max_attempts,
                delay)
            time.sleep(delay)

    raise RuntimeError(f"Region {region_name!r} was not returned by QAI after create/update.")


def fetch_pipeline_by_version(session, pipeline_version, max_attempts=5, wait_seconds=0.2):
    logger.debug('Fetching pipeline by version from QAI: %s', pipeline_version)
    pipeline_path = '/lab_miseq_pipelines?version=' + pipeline_version
    for attempt in range(1, max_attempts + 1):
        pipelines = session.get_json(pipeline_path)
        if pipelines:
            if len(pipelines) > 1:
                logger.warning('Pipeline lookup for %s returned %d rows; using the first.',
                               pipeline_version,
                               len(pipelines))
            return pipelines[0]

        if attempt < max_attempts:
            delay = wait_seconds * (2 ** (attempt - 1))
            logger.debug(
                'Pipeline %s not visible yet (attempt %d/%d). Retrying in %.2fs.',
                pipeline_version,
                attempt,
                max_attempts,
                delay)
            time.sleep(delay)

    raise RuntimeError(
        f"Pipeline {pipeline_version!r} was not returned by QAI after create.")


def parse_args():
    # noinspection PyTypeChecker
    parser = ArgumentParser(description='Upload project definitions to QAI.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--qai_server',
        default=os.environ.get('MICALL_QAI_SERVER', 'http://localhost:4567'),
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

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true',
                                 help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def main() -> int:
    args = parse_args()
    configure_logging(args)
    logger.info('Starting project upload.')
    logger.debug(
        'Command options: qai_server=%s qai_user=%s pipeline_version=%s update_sequences=%s',
        args.qai_server,
        args.qai_user,
        args.pipeline_version,
        args.update_sequences)

    project_config = ProjectConfig.loadDefault()
    scoring_path = ProjectsScoringFile().path()
    logger.debug('Loading scoring config from %s', scoring_path)
    with scoring_path.open() as scoring_file:
        scoring_config = json.load(scoring_file)
    logger.debug(
        'Loaded %d projects and %d regions from projects.json, %d scoring projects from project_scoring.json.',
        len(project_config.config['projects']),
        len(project_config.config['regions']),
        len(scoring_config.get('projects', {})))

    missing_scoring_regions = find_missing_scoring_regions(project_config, scoring_config)
    if missing_scoring_regions:
        missing_summary = ', '.join(
            f'{project_name}/{region_name}'
            for project_name, region_name in sorted(missing_scoring_regions)
        )
        logger.error('Scoring config is missing amino-acid regions: %s', missing_summary)
        return 1
    logger.debug('Scoring validation passed for all amino-acid coordinate regions.')

    with qai_helper.Session() as session:
        logger.debug('Logging into QAI server %s as user %s.', args.qai_server, args.qai_user)
        session.login(args.qai_server,
                      args.qai_user,
                      args.qai_password)
        logger.debug('Login successful.')

        pipelines = session.get_json(
            "/lab_miseq_pipelines?version=" + args.pipeline_version)
        logger.debug('Pipeline lookup for %s returned %d entries.',
                     args.pipeline_version,
                     len(pipelines))
        if pipelines:
            logger.error('Pipeline %s already exists.', args.pipeline_version)
            return 1

        seed_groups = session.get_json("/lab_miseq_seed_groups")
        logger.debug('Fetched %d seed groups.', len(seed_groups))
        # noinspection PyTypeChecker
        seed_group_ids = dict(map(itemgetter('name', 'id'), seed_groups))
        old_regions = session.get_json("/lab_miseq_regions")
        logger.debug('Fetched %d existing regions from QAI.', len(old_regions))
        regions = dict(((region['name'], region) for region in old_regions))
        for region_name, region_data in project_config.config['regions'].items():
            ref_seq = ''.join(region_data['reference'])
            region = regions.get(region_name)
            logger.debug('Processing region %s (is_nucleotide=%s).',
                         region_name,
                         region_data['is_nucleotide'])
            if region is None:
                seed_group_name = region_data['seed_group']
                seed_group_id = seed_group_ids.get(seed_group_name)
                if seed_group_id is None and seed_group_name:
                    logger.debug('Creating new seed group %s for region %s.',
                                 seed_group_name,
                                 region_name)
                    seed_group = session.post_json("/lab_miseq_seed_groups",
                                                   {'name': seed_group_name})
                    seed_group_id = seed_group['id']
                    seed_group_ids[seed_group_name] = seed_group_id
                logger.info("Creating region %s with seed group %s", region_name, seed_group_name)
                region = session.post_json(
                    "/lab_miseq_regions",
                    {'name': region_name,
                     'is_nucleotide': region_data['is_nucleotide'],
                     'reference': ref_seq,
                     'seed_group_id': seed_group_id})
                if not region:
                    logger.debug('Region create for %s returned empty payload; fetching by name.',
                                 region_name)
                    region = fetch_region_by_name(session, region_name)
                regions[region_name] = region
                logger.debug('Region %s now tracked with id %s.', region_name, region.get('id'))
            elif region['reference'] != ref_seq:
                logger.warning("Reference doesn't match: %s", region_name)
                if args.update_sequences:
                    logger.debug('Updating reference sequence for region %s (id=%s).',
                                 region_name,
                                 region.get('id'))
                    region['reference'] = ref_seq
                    session.post_json(f"/lab_miseq_regions/{region['id']}",
                                      region)
                    refreshed_region = fetch_region_by_name(session, region_name)
                    regions[region_name] = refreshed_region
                else:
                    logger.debug('Skipping reference update for %s because --update_sequences is not set.',
                                 region_name)

        logger.info("Uploading project definitions to QAI server %s", args.qai_server)
        pipeline = session.post_json("/lab_miseq_pipelines",
                                     {'version': args.pipeline_version})
        if not pipeline:
            logger.debug('Pipeline create for %s returned empty payload; fetching by version.',
                         args.pipeline_version)
            pipeline = fetch_pipeline_by_version(session, args.pipeline_version)
        logger.info("Created pipeline version %s with id %d", args.pipeline_version, pipeline['id'])
        pipeline_id = pipeline['id']

        old_projects = session.get_json("/lab_miseq_projects")
        logger.debug('Fetched %d existing projects from QAI.', len(old_projects))
        projects = dict(((project['name'], project) for project in old_projects))
        for project_name, project_data in project_config.config['projects'].items():
            logger.debug('Uploading project %s with %d region mappings.',
                         project_name,
                         len(project_data['regions']))
            project = projects.get(project_name)
            if project is None:
                logger.debug('Project %s not found on QAI; creating it.', project_name)
                project = session.post_json(
                    "/lab_miseq_projects",
                    {'name': project_name,
                     'max_variants': project_data['max_variants']})
            else:
                logger.debug('Project %s already exists with id %s.', project_name, project.get('id'))
            project_version = session.post_json("/lab_miseq_project_versions",
                                                {'pipeline_id': pipeline_id,
                                                 'project_id': project['id']})
            logger.debug('Created project version id %s for project %s.',
                         project_version.get('id'),
                         project_name)
            for region_data in project_data['regions']:
                region_name = region_data['coordinate_region']
                if regions[region_name]["is_nucleotide"]:
                    logger.debug('Skipping nucleotide coordinate region %s for project %s.',
                                 region_name,
                                 project_name)
                    continue
                for region in scoring_config['projects'][project_name]['regions']:
                    if region['coordinate_region'] == region_name:
                        scoring_data = region
                        break
                else:
                    raise ValueError(f"Region {region_name!r} not present in scoring config file")
                coordinate_region = regions[region_name]
                seed_region = regions[region_data['seed_region_names'][0]]
                seed_group_id = seed_region['seed_group_id']
                logger.debug(
                    'Uploading project region for project=%s coordinate=%s seed=%s min_coverage=%s/%s/%s',
                    project_name,
                    region_name,
                    region_data['seed_region_names'][0],
                    scoring_data['min_coverage1'],
                    scoring_data['min_coverage2'],
                    scoring_data['min_coverage3'])
                project_region = session.post_json(
                    "/lab_miseq_project_regions",
                    {'project_version_id': project_version['id'],
                     'coordinate_region_id': coordinate_region['id'],
                     'min_coverage1': scoring_data['min_coverage1'],
                     'min_coverage2': scoring_data['min_coverage2'],
                     'min_coverage3': scoring_data['min_coverage3'],
                     'seed_group_id': seed_group_id})

                logger.debug('Project region id %s created with %d key positions.',
                             project_region.get('id'),
                             len(scoring_data['key_positions']))
                for key_position in scoring_data['key_positions']:
                    session.post_json("/lab_miseq_key_positions",
                                      {'project_region_id': project_region['id'],
                                       'start_pos': key_position['start_pos'],
                                       'end_pos': key_position['end_pos']})
                logger.debug('Finished key positions for project=%s coordinate=%s.',
                             project_name,
                             region_name)

    logger.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
