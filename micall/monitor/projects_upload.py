from operator import itemgetter

from micall import settings
from micall.core.project_config import ProjectConfig
from micall.monitor import qai_helper


def main():
    project_config = ProjectConfig.loadDefault()
    with qai_helper.Session() as session:
        session.login(settings.qai_path,
                      settings.qai_user,
                      settings.qai_password)

        pipelines = session.get_json(
            "/lab_miseq_pipelines?version=" + settings.pipeline_version,
            retries=0)
        if pipelines:
            raise RuntimeError('Pipeline {} already exists.'.format(
                settings.pipeline_version))

        seed_groups = session.get_json("/lab_miseq_seed_groups")
        seed_group_ids = dict(map(itemgetter('name', 'id'), seed_groups))
        old_regions = session.get_json("/lab_miseq_regions", retries=0)
        regions = dict(((region['name'], region) for region in old_regions))
        for region_name, region_data in project_config.config['regions'].iteritems():
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
                     'reference': ''.join(region_data['reference']),
                     'seed_group_id': seed_group_id})
                regions[region_name] = region

        pipeline = session.post_json("/lab_miseq_pipelines",
                                     {'version': settings.pipeline_version})
        pipeline_id = pipeline['id']

        old_projects = session.get_json("/lab_miseq_projects", retries=0)
        projects = dict(((project['name'], project) for project in old_projects))
        for project_name, project_data in project_config.config['projects'].iteritems():
            project = projects.get(project_name)
            if project is None:
                project = session.post_json(
                    "/lab_miseq_projects",
                    {'name': project_name,
                     'max_variants': project_data['max_variants']})
            project_version = session.post_json("/lab_miseq_project_versions",
                                                {'pipeline_id': pipeline_id,
                                                 'project_id': project['id']})
            for region_data in project_data['regions']:
                coordinate_region = regions[region_data['coordinate_region']]
                seed_region = regions[region_data['seed_region_names'][0]]
                seed_group_id = seed_region['seed_group_id']
                project_region = session.post_json(
                    "/lab_miseq_project_regions",
                    {'project_version_id': project_version['id'],
                     'coordinate_region_id': coordinate_region['id'],
                     'min_coverage1': region_data['min_coverage1'],
                     'min_coverage2': region_data['min_coverage2'],
                     'min_coverage3': region_data['min_coverage3'],
                     'seed_group_id': seed_group_id})

                for key_position in region_data['key_positions']:
                    session.post_json("/lab_miseq_key_positions",
                                      {'project_region_id': project_region['id'],
                                       'start_pos': key_position['start_pos'],
                                       'end_pos': key_position['end_pos']})

    print "Done."

if __name__ == "__main__":
    main()
