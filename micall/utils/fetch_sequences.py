import re
from difflib import Differ
from io import StringIO
from itertools import chain
from time import sleep

# noinspection PyUnresolvedReferences
from gotoh import align_it_aa
import requests

from micall.core import aln2counts
from micall.g2p.fastq_g2p import GAP_OPEN_COST, GAP_EXTEND_COST, USE_TERMINAL_COST
from micall.core.project_config import ProjectConfig

# From pssmlib
PSSM_V3LOOP = "CTRPNXNNTXXRKSIRIXXXGPGQXXXAFYATXXXXGDIIGDIXXRQAHC"
PSSM_V3LOOP = PSSM_V3LOOP.replace('X', '')


def find_best_sequence(sequences):
    best_match = (0, None, None)  # (score, name, sequence)
    for name, nuc_seq in sorted(sequences.items()):
        nuc_seq = nuc_seq.replace('\n', '')
        nuc_seq = nuc_seq.replace('*', '-').replace('?', '-')
        for frame in range(3):
            offset_seq = '-' * frame + nuc_seq
            aa_seq = aln2counts.translate(offset_seq)
            aseq, aref, score = align_it_aa(aa_seq,
                                            PSSM_V3LOOP,
                                            GAP_OPEN_COST,
                                            GAP_EXTEND_COST,
                                            USE_TERMINAL_COST)
            best_match = max(best_match, (score, name, aseq, aref))

            if score >= 272:
                pairs = zip(aseq, aref)
                diffs = [' ' if a == b else '*' for a, b in pairs]
                print('score', score, name)
                print('result ', aseq)
                print('diffs  ', ''.join(diffs) if aseq != aref else 'no diffs')
                print('compare', aref)
    return best_match[:2]


def main():
    # find_best_match_for_pssm()
    sequences = fetch_alignment_sequences(2004,
                                          'CON',  # Consensus/Ancestral
                                          'ENV')
    consensus = sequences['CON_OF_CONS'].replace('-', '').upper()

    project_config = ProjectConfig.loadDefault()
    ref_names = set(project_config.getAllReferences().keys())
    new_sequences = fetch_alignment_sequences('2015', 'COM')
    consensus_accession = 'Consensus'
    assert consensus_accession not in new_sequences, sorted(new_sequences.keys())
    new_sequences[consensus_accession] = consensus

    for line in compare_config('HIV', project_config, new_sequences, ref_names):
        print(line, end='')

    print('Unchecked refs: ' + ', '.join(sorted(ref_names)))


def find_best_match_for_pssm():
    scenarios = (('ENV', '1997 1999 2000 2001 2002 2002JAN 2002AUG 2004'),
                 ('GENOME', '2000 2001 2002'))
    for region, years in scenarios:
        for year in years.split():
            sequences = fetch_alignment_sequences(year,
                                                  'CON',  # Consensus/Ancestral
                                                  region)
            score, name = find_best_sequence(sequences)
            print(score, name, year, region)
            sleep(1)


def fetch_alignment_sequences(year, alignment_type, region='GENOME'):
    response = requests.post(
        'https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi',
        data={'ORGANISM': 'HIV',
              'ALIGN_TYPE': alignment_type,
              'SUBORGANISM': 'HIV1',
              'PRE_USER': 'predefined',
              'REGION': region,
              'GENO_SUB': 'All',
              'BASETYPE': 'DNA',
              'YEAR': year,
              'FORMAT': 'fasta',
              'submit': 'Get Alignment'})
    match = re.search(r'<a href="([^"]*)">Download', response.text)
    assert match
    response = requests.get('https://www.hiv.lanl.gov' + match.group(1))
    sequences = parse_compendium(StringIO(response.text))
    return sequences


def fetch_by_accession(accession):
    response = requests.get(
        'https://www.hiv.lanl.gov/components/sequence/HIV/asearch/'
        'query_one.comp?se_id=' + accession)
    return extract_sequence(response.text)


def extract_sequence(source):
    """ Extract a reference's sequence from its search page.

    :param str source: HTML from the search result page
    """
    match = re.search(r'ORIGIN(.*?)//', source, re.DOTALL)
    assert match, source
    sections = []
    for line in match.group(1).splitlines():
        if not line:
            continue
        line_match = re.match(r'^\s*\d+\s*([acgtryswkmbdhvn\s]+)$', line)
        assert line_match, line
        section = line_match.group(1)
        sections.append(section.replace(' ', '').upper())
    return ''.join(sections)


def parse_compendium(fasta):
    key = None
    sequences = {}
    section = None
    for line in chain(fasta, [None]):
        is_header = line is None or line.startswith('>')
        if not is_header:
            section.append(line.strip())
        else:
            if key is not None:
                sequences[key] = ''.join(section)
            section = []
            if line is not None:
                key = line[1:].strip().split('.')[-1]

    return sequences


def compare_config(project, project_config, sequences, ref_names=None):
    differ = Differ()
    diff_count = 0
    for ref_name in sorted(project_config.getProjectSeeds(project)):
        yield ref_name + '\n'
        if ref_names is not None:
            ref_names.remove(ref_name)
        old_sequence = project_config.getReference(ref_name)
        name_parts = ref_name.split('-')
        accession = name_parts[3]
        new_sequence = sequences.get(accession)
        new_sequence = new_sequence and new_sequence.replace('-', '')
        if new_sequence is None:
            diff_count += 1
            yield 'missing\n'
        elif new_sequence == old_sequence:
            yield '==\n'
        else:
            diff_count += 1
            diff = differ.compare([old_sequence+'\n'], [new_sequence+'\n'])
            yield from diff
            
        yield '\n'
    yield '{} differences\n'.format(diff_count)


if __name__ == '__main__':
    main()
