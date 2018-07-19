"""
Provides methods for inspecting and editing a MiCall project configuration
JSON file.
"""

from micall.core.project_config import ProjectConfig
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Inspect and edit a MiCall project configuration JSON file.'
    )
    parser.add_argument('json',
                        help='Path to JSON file (default: ./projects.json)',
                        default='projects.json')
    parser.add_argument('target', choices=['seeds', 'coords', 'all'],
                        help='Target seed or coordinate references',
                        default='all')

    actions = parser.add_mutually_exclusive_group()
    actions.add_argument('--list', help='List all ')

    parser.parse_args()


def main():
    args = parseArgs()
    print(args)

if __name__ == '__main__':
    main()
