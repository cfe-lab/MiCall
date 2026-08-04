import json
from typing import Dict, List

from micall.utils.externals import ProjectsFile, ProjectsScoringFile

G2P_SEED_NAME = "HIV1-CON-XX-Consensus-seed"


class ProjectConfig(object):
    @classmethod
    def search(cls, project_paths):
        projects = None
        for project_path in project_paths:
            try:
                with open(project_path, 'r') as projects_file:
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
        with ProjectsFile().path() as file_path:
            project_paths = [file_path]
            return cls.search(project_paths)

    @classmethod
    def loadScoring(cls):
        with ProjectsScoringFile().path() as file_path:
            project_paths = [file_path]
            return cls.search(project_paths)

    def load(self, json_file):
        self.config = json.load(json_file)

    def writeSeedFasta(self, fasta_file, excluded_seeds=None):
        """ Write seed references to a FASTA file.

        @param fasta_file: an open file
        @param excluded_seeds: a list of seed names to exclude from the file
        """
        seed_region_set = set()
        for project in self.config['projects'].values():
            for region in project['regions']:
                seed_region_set.update(region['seed_region_names'])

        if excluded_seeds:
            seed_region_set.difference_update(excluded_seeds)
        seed_region_list = list(seed_region_set)
        seed_name_map: Dict[str, str] = {}  # {sequence: name}
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
        return ''.join(reference)

    def getGenotypeReference(self, region_name):
        reference = self.config['genotype_references'][region_name]['reference']
        return ''.join(reference)

    def getNucReference(self, reference_name, region_start, region_end):
        try:
            reference = self.config['genotype_references'][reference_name]['reference']
        except KeyError:
            reference = self.config['regions'][reference_name]['reference']
        reference_full = ''.join(reference)
        region_sequence = reference_full[region_start-1:region_end]
        return region_sequence

    def isAmino(self, region_name):
        return not self.config['regions'][region_name]['is_nucleotide']

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

    def getAllReferences(self):
        return {name: self.getReference(name)
                for name in self.config['regions']}

    def getAllGenotypeReferences(self):
        return {name: self.getGenotypeReference(name)
                for name in self.config['genotype_references']}

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

    def getProjectRegions(self, seed_name, coordinate_name, excluded_projects=None):
        """ Find all project regions that use the seed and coordinate region.

        @param seed_name: the needed seed reference
        @param coordinate_name: the needed coordinate reference
        @param excluded_projects: a list of project names to exclude
        @return: a generator of project configuration dictionaries
        """
        project_names = set(self.config['projects'])
        if excluded_projects is not None:
            project_names.difference_update(excluded_projects)
        project_names: List[str] = sorted(project_names) # type: ignore[no-redef]
        for project_name in project_names:
            project = self.config['projects'][project_name]
            for region in project['regions']:
                if region['coordinate_region'] == coordinate_name and seed_name in region['seed_region_names']:
                    project_region = dict(region)
                    del project_region['coordinate_region']
                    del project_region['seed_region_names']
                    project_region['project_name'] = project_name
                    yield project_region
