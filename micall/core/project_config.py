"""
Defines an object class to parse a JSON object containing project
specifications -- e.g., seed and coordinate reference sequences.
"""

import json
import os


class ProjectConfig(object):
    @classmethod
    def search(cls, project_paths):
        projects = None
        for project_path in project_paths:
            try:
                with open(project_path, 'rU') as projects_file:
                    projects = cls()
                    projects.load(projects_file)
                    break
            except:
                projects = None

        if not projects:
            raise RuntimeError('No project definitions found in {!r}'.format(
                project_paths))
        return projects

    @classmethod
    def loadDefault(cls):
        file_path = os.path.dirname(__file__)
        project_paths = [os.path.join(file_path, 'projects.json'),
                         os.path.join(os.path.dirname(file_path), 'projects.json')]
        return cls.search(project_paths)

    def load(self, json_file):
        self.json_file = json_file.name
        self.config = json.load(json_file)

    def writeSeedFasta(self, fasta_file):
        """ Write seed references to a FASTA file.

        @param fasta_file: an open file
        """
        seed_region_set = set()
        for project in self.config['projects'].values():
            for region in project['regions']:
                seed_region_set.update(region['seed_region_names'])

        seed_region_list = list(seed_region_set)
        seed_name_map = {}  # {sequence: name}
        seed_region_list.sort()
        for name in seed_region_list:
            region = self.config['regions'][name]
            sequence = ''.join(region['reference'])
            duplicate_name = seed_name_map.get(sequence)
            if duplicate_name is not None:
                raise RuntimeError("Duplicate references: {} and {}.".format(
                    duplicate_name,
                    name))
            seed_name_map[sequence] = name
            fasta_file.write('>{name}\n{ref}\n'.format(name=name,
                                                       ref=sequence))

    def getReference(self, region_name):
        reference = self.config['regions'][region_name]['reference']
        return ''.join(reference).encode('utf-8')

    def getCoordinateReferences(self, seed_region):
        """ Find any coordinate references that are linked to a seed reference.

        @param seed_region: the name of a seed region
        @return {name: sequence} for any coordinate references, or {} if there
            are no linked references.
        """
        coord_refs = {}
        for project in self.config['projects'].values():
            for region in project['regions']:
                coord_region = region['coordinate_region']
                if seed_region in region['seed_region_names'] and coord_region:
                    coord_refs[coord_region] = self.getReference(coord_region)
        return coord_refs

    def getMaxVariants(self, coordinate_region):
        """ Find the maximum number of variants to report for a coordinate
        region.

        @param coordinate_region: The name of a coordinate region
        @return an integer with the maximum variants requested for any project
            that uses the coordinate region.
        """
        max_variants = 0
        for project in self.config['projects'].values():
            for region in project['regions']:
                if region['coordinate_region'] == coordinate_region:
                    max_variants = max(project['max_variants'], max_variants)
        return max_variants

    def getProjectSeeds(self, project_name):
        """ Return all the seed regions used by a project.

        @return a set of seed region names
        """

        seeds = set()
        for region in self.config['projects'][project_name]['regions']:
            seeds.update(region['seed_region_names'])

        return seeds

    def getSeedGroup(self, seed_region):
        """ Find the seed group a seed region belongs to.

        @param seed_region: the name of a seed region
        @return the name of a seed group
        """

        return self.config['regions'][seed_region]['seed_group']

    def getProjectRegions(self, seed_name, coordinate_name):
        """ Find all projects that use a specific combination of seed and
        coordinate references.

        @param seed_name: the name of a seed region
        @param coordinate_name: the name of a coordinate region
        @yield project names
        """
        project_names = self.config['projects'].keys()
        project_names.sort()
        for project_name in project_names:
            project = self.config['projects'][project_name]
            for region in project['regions']:
                if region['coordinate_region'] == coordinate_name and seed_name in region['seed_region_names']:
                    project_region = dict(region)
                    del project_region['coordinate_region']
                    del project_region['seed_region_names']
                    project_region['project_name'] = project_name
                    yield project_region

if __name__ == '__live_coding__':
    import unittest
    from micall.tests.project_config_test import ProjectConfigurationTest

    suite = unittest.TestSuite()
    suite.addTest(ProjectConfigurationTest("testProjectRegions"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
