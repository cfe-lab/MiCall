import json
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from operator import itemgetter
import logging
from micall.utils.externals import ProjectsFile, ProjectsScoringFile

from collections import Counter
from copy import deepcopy
try:
    from micall.monitor import qai_helper
except ImportError:
    # Ignore import errors to allow tests without request module
    qai_helper = None


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def parse_args():
    parser = ArgumentParser(description='Dump project definitions from QAI.',
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

    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def check_key_positions(projects, warning_file):
    """ Complain if a coordinate region has two sets of key positions in the
    same project.

    @param projects: a dictionary of project definitions
    @param warning_file: an open file to write any warnings in
    """
    key_warning = ("WARNING: project {} has multiple sets of key positions for "
                   "coordinate region {}.\n")
    seed_warning = "WARNING: project {} has duplicate seed and coordinate: {}, {}\n"
    warnings = []
    for project_name, project in projects.items():
        key_position_coordinates = Counter()
        seed_coordinate_pairs = Counter()
        for region in project['regions']:
            coordinate_name = region['coordinate_region']
            if region['key_positions']:
                key_position_coordinates[coordinate_name] += 1
            for seed_name in region['seed_region_names']:
                seed_coordinate_pairs[(seed_name, coordinate_name)] += 1
        sorted_counts = key_position_coordinates.most_common()
        for coordinate_name, count in sorted_counts:
            if count > 1:
                warnings.append(key_warning.format(project_name,
                                                   coordinate_name))
            else:
                break
        sorted_counts = seed_coordinate_pairs.most_common()
        for (seed_name, coordinate_name), count in sorted_counts:
            if count > 1:
                warnings.append(seed_warning.format(project_name,
                                                    seed_name,
                                                    coordinate_name))
            else:
                break
    for warning in sorted(warnings):
        warning_file.write(warning)


def dump_json(json_object, filename):
    with open(filename, "w") as f:
        json.dump(json_object,
                  f,
                  sort_keys=True,
                  indent=2,
                  separators=(',', ': '))
        f.write('\n')


def main():
    args = parse_args()
    dump = {}
    used_regions = set()
    with qai_helper.Session() as session:
        session.login(args.qai_server,
                      args.qai_user,
                      args.qai_password)

        logger.info("Dumping project definitions from QAI server %s", args.qai_server)

        dump['regions'] = session.get_json("/lab_miseq_regions?mode=dump")
        dump['projects'] = session.get_json(
            "/lab_miseq_projects?mode=dump&pipeline=" +
            args.pipeline_version)

        regions_count = len(dump['regions'])
        projects_count = len(dump['projects'])
        logger.info("Dumped %d regions and %d projects", regions_count, projects_count)

        empty_projects = []
        for name, project in dump['projects'].items():
            project['regions'].sort(key=itemgetter('coordinate_region'))
            for region in project['regions']:
                used_regions.add(region['coordinate_region'])
                used_regions.update(region['seed_region_names'])
            if not project['regions']:
                empty_projects.append(name)
        for name in empty_projects:
            del dump['projects'][name]
        errors = dump['projects'].get('errors')
        if errors:
            raise RuntimeError('\n'.join(errors))
        check_key_positions(dump['projects'], sys.stdout)
    dump['regions'] = {key: value
                       for key, value in dump['regions'].items()
                       if key in used_regions}

    dump_scoring = deepcopy(dump)
    for project in dump['projects'].values():
        for region in project['regions']:
            del region['key_positions']
            del region['min_coverage1']
            del region['min_coverage2']
            del region['min_coverage3']

    with ProjectsFile().path() as projects_file_path:
        dump_json(dump, projects_file_path)
        logger.info("Wrote %s", projects_file_path)

    for project in dump_scoring['projects'].values():
        for region in project['regions']:
            name = region['coordinate_region']
            seq = ''.join(dump_scoring['regions'][name]['reference'])
            region['coordinate_region_length'] = len(seq)
    del dump_scoring['regions']

    with ProjectsScoringFile().path() as projects_scoring_file_path:
        dump_json(dump_scoring, projects_scoring_file_path)
        logger.info("Wrote %s", projects_scoring_file_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
