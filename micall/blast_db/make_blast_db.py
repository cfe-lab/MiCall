import json
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from pathlib import Path
from subprocess import run

DEFAULT_DATABASE = str(Path(__file__).parent / 'refs.fasta')
DEFAULT_PROJECTS = str(Path(__file__).parent.parent / 'projects.json')
COMPLETE_HIV_REFS = """\
HIV1-A1-CD-AM000053-seed
HIV1-B-FR-K03455-seed
HIV1-C-IN-KC156210-seed
HIV1-C-MW-KC156214-seed
HIV1-CPZ-TZ-JN091691-seed
HIV1-CPZ-US-AF103818-seed
HIV1-P-FR-GU111555-seed
HIV1-C-TZ-KC156220-seed
HIV1-CPZ-CM-FR686511-seed
HIV1-O-SN-AJ302646-seed
HIV1-O-BE-L20587-seed
HIV1-D-SN-AB485648-seed
HIV1-F1-RO-AB485658-seed
HIV1-G-GH-AB287004-seed
HIV1-D-KR-DQ054367-seed
HIV1-G-PT-FR846409-seed
HIV1-B-US-KT284371-seed
HIV1-A1-IN-KT152839-seed""".splitlines()


def parse_args():
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
            if (region['seed_group'] == 'HIV1-seed' and
                    name not in COMPLETE_HIV_REFS):
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
