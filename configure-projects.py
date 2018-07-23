"""
Provides methods for inspecting and editing a MiCall project configuration
JSON file.

A project is a set of seed and coordinate reference sequences.
A seed is a sequenced that is used to initialize the iterative mapping
  of short reads.
A seed group is a set of seed references.
A coordinate reference is a sequence that defines a numbering system for
  the positions of nucleotides and amino acids.

"""

from micall.core.project_config import ProjectConfig
import argparse
import sys
import os

class ProjectEditor(ProjectConfig):
    def list_projects(self):
        """ Summary of all projects """
        #sys.stdout.write('\n\033[1mPROJECTS\033[0m\n')
        sys.stdout.write('{:>10}{:>9}{:>10}  {}\n'.format(
            'Name', 'Regions', 'Variants', 'Description')
        )
        project_names = list(self.config['projects'].keys())
        project_names.sort()

        for project_name in project_names:
            info = self.config['projects'][project_name]
            sys.stdout.write('{:>10}{:>9}{:>10}  {}{}\n'.format(
                project_name, len(info['regions']), info['max_variants'],
                info['description'][:40],
                '...' if len(info['description']) > 40 else ''
            ))
        sys.stdout.write('\n')



def parseArgs():
    parser = argparse.ArgumentParser(
        description='Inspect and edit a MiCall project configuration JSON file.'
    )
    parser.add_argument('-json',
                        help='Path to JSON file (default: micall/projects.json)')
    return parser.parse_args()


def list_project(pcfg, project, limit=3):
    for region in regions[:limit]:
        sys.stdout.write('    {}\n'.format(region['coordinate_region']))
    if len(regions) > limit:
        sys.stdout.write('    ({} others)\n'.format(len(regions) - limit))








def list_seeds(pcfg):
    seeds = pcfg.getProjectSeeds()
    return seeds



def list_cmds():
    sys.stdout.write('q, list')


def main():
    args = parseArgs()

    pcfg = ProjectEditor.loadDefault() if args.json is None else \
        ProjectEditor.load(args.json)

    sys.stdout.write('\033[1mMiCall project editor\033[0m\n')
    sys.stdout.write('loaded projects from file {}\n'.format(
        os.path.relpath(pcfg.json_file))
    )

    while True:
        sys.stdout.write('> ')
        cmd = input()
        if cmd == 'q' or cmd == 'quit':
            break

        action = cmd
        target = None
        if ' ' in cmd:
            action, target = cmd.split()

        if action == 'list':
            if target is None:
                pcfg.list_projects()
        elif action == 'add':
            pass
        elif action == 'edit':
            pass
        elif action == 'delete':
            pass


if __name__ == '__main__':
    main()
