"""
Provides methods for inspecting and editing a MiCall project configuration
JSON file.
"""

from micall.core.project_config import ProjectConfig
import argparse
import sys

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Inspect and edit a MiCall project configuration JSON file.'
    )
    parser.add_argument('-json',
                        help='Path to JSON file (default: micall/projects.json)')

    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--list', help='List all ',
                         choices=['projects', 'seeds', 'coords', 'all'])

    return parser.parse_args()

def list_projects(pcfg):
    sys.stdout.write('Projects:\n------------------\n')
    for project_name, info in pcfg.config['projects'].items():
        sys.stdout.write('%s : \n'.format(project_name))


def list_seeds(pcfg):
    seeds = pcfg.getProjectSeeds()
    return seeds

def main():
    args = parseArgs()
    pcfg = ProjectConfig.loadDefault() if args.json is None else \
        ProjectConfig.load(args.json)
    if args.list:
        if args.list == 'projects':
            list_projects(pcfg)


if __name__ == '__main__':
    main()
