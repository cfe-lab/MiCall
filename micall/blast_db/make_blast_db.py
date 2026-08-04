import json
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from pathlib import Path
from subprocess import run

from micall.utils.externals import ProjectsFile

def parse_args():
    DEFAULT_DATABASE = str(Path(__file__).parent / 'refs.fasta')

    with ProjectsFile().path() as projects_file_path:
        DEFAULT_PROJECTS = str(projects_file_path)

    parser = ArgumentParser(description='Build Blast database from FASTA.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('json',
                        help='JSON file with projects and reference sequences',
                        nargs='?',
                        default=DEFAULT_PROJECTS,
                        type=FileType())
    parser.add_argument('fasta',
                        help='FASTA file with reference sequences',
                        nargs='?',
                        default=DEFAULT_DATABASE,
                        type=FileType('w'))
    return parser.parse_args()


def make_blast_db(projects_json, refs_fasta):
    with projects_json, refs_fasta:
        projects = json.load(projects_json)
        for name, region in projects['regions'].items():
            if region['seed_group'] is None:
                continue
            if name.startswith('HIV1-CRF'):
                # Exclude circulating recombinant forms (CRF).
                continue
            if name == 'HIV1-CON-XX-Consensus-seed':
                # Only used by G2P alignment.
                continue
            refs_fasta.write('>' + name + '\n')
            for line in region['reference']:
                refs_fasta.write(line)
                refs_fasta.write('\n')

    run(['makeblastdb',
         '-in', refs_fasta.name,
         '-parse_seqids',
         '-dbtype', 'nucl'],
        check=True)


def main():
    args = parse_args()
    make_blast_db(args.json, args.fasta)


if __name__ == '__main__':
    main()
