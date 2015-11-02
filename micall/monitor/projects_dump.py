import json
import sys

from micall.monitor import qai_helper
from micall import settings
from collections import Counter


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
    for project_name, project in projects.iteritems():
        regions_with_key_positions = set()
        seed_coordinate_pairs = Counter()
        for region in project['regions']:
            coordinate_name = region['coordinate_region']
            if region['key_positions']:
                if coordinate_name in regions_with_key_positions:
                    warnings.append(key_warning.format(project_name,
                                                       coordinate_name))
                regions_with_key_positions.add(coordinate_name)
            for seed_name in region['seed_region_names']:
                seed_coordinate_pairs[(seed_name, coordinate_name)] += 1
        sorted_counts = seed_coordinate_pairs.most_common()
        for (seed_name, coordinate_name), count in sorted_counts:
            if count > 1:
                warnings.append(seed_warning.format(seed_name, coordinate_name))
            else:
                break
    for warning in sorted(warnings):
        warning_file.write(warning)


def main():
    dump = {}
    with qai_helper.Session() as session:
        session.login(settings.qai_project_path,
                      settings.qai_project_user,
                      settings.qai_project_password)

        dump['regions'] = session.get_json("/lab_miseq_regions.json?mode=dump",
                                           retries=0)
        dump['projects'] = session.get_json(
            "/lab_miseq_projects.json?mode=dump&pipeline=" +
            settings.pipeline_version,
            retries=0)
        for project in dump['projects'].itervalues():
            project['regions'].sort()
        errors = dump['projects'].get('errors')
        if errors:
            raise StandardError('\n'.join(errors))
        check_key_positions(dump['projects'], sys.stdout)

    with open("../projects.json", "w") as f:
        json.dump(dump, f, sort_keys=True, indent=2, separators=(',', ': '))
        f.write('\n')

    print "Done."

if __name__ == "__main__":
    main()
