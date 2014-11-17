import json
import os

import settings


class ProjectConfig(object):
    @classmethod
    def loadDefault(cls):
        projects = None
        project_paths = [settings.projects_json,
                         os.path.basename(settings.projects_json)]
        for project_path in project_paths:
            try:
                with open(project_path, 'rU') as projects_file:
                    projects = ProjectConfig()
                    projects.load(projects_file)
                    break
            except:
                pass
        
        if not projects:
            raise RuntimeError('No project definitions found in {!r}'.format(
                project_paths))
        return projects
    
    def load(self, json_file):
        self.config = json.load(json_file)
    
    def writeSeedFasta(self, fasta_file):
        seed_region_set = set()
        for project in self.config['projects'].itervalues():
            for region in project['regions']:
                seed_region_set.add(region['seed_region'])
        
        seed_region_list = list(seed_region_set)
        seed_region_list.sort()
        for name in seed_region_list:
            region = self.config['regions'][name]
            sequence = ''.join(region['reference'])
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
        for project in self.config['projects'].itervalues():
            for region in project['regions']:
                coord_region = region['coordinate_region']
                if region['seed_region'] == seed_region and coord_region:
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
        for project in self.config['projects'].itervalues():
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
            seeds.add(region['seed_region'])
        
        return seeds
