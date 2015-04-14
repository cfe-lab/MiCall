import json
import sys

import qai_helper
import settings

def check_key_positions(projects, warning_file):
    """ Complain if a coordinate region has two sets of key positions in the
    same project.
    
    @param projects: a dictionary of project definitions
    @param warning_file: an open file to write any warnings in
    """
    warnings = set() # set([(project, region)])
    for project_name, project in projects.iteritems():
        regions_with_key_positions = set()
        for region in project['regions']:
            if region['key_positions']:
                coordinate_name = region['coordinate_region']
                if coordinate_name in regions_with_key_positions:
                    warnings.add((project_name, coordinate_name))
                regions_with_key_positions.add(coordinate_name)
    for project_name, coordinate_name in sorted(warnings):
        warning_file.write(
            ("WARNING: project {} has multiple sets of key " +
            "positions for coordinate region {}.\n").format(
                project_name,
                coordinate_name))
    

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
        errors = dump['projects'].get('errors')
        if errors:
            raise StandardError('\n'.join(errors))
        check_key_positions(dump['projects'], sys.stdout)
        
    with open("projects.json", "w") as f:
        json.dump(dump, f, sort_keys=True, indent=2, separators=(',', ': '))
        f.write('\n')
    
    print "Done."

if __name__ == "__main__":
    main()
