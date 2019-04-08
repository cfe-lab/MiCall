""" Report where references come from, and validate that they match the source.

Running this script reports a source for each group of references, downloads
them again, and checks that they match the source.
Finally, it reports any references that aren't documented yet.
"""

import re
from difflib import Differ
from io import StringIO
from itertools import chain
from textwrap import fill
from time import sleep

# noinspection PyUnresolvedReferences
from gotoh import align_it_aa, align_it

from micall.utils.translation import translate

try:
    import requests
except ImportError:
    # Allow tests to run without requests module
    requests = None

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
    project_config = ProjectConfig.loadDefault()
    error_count = 0
    unchecked_ref_names = set(project_config.getAllReferences().keys())
    # find_best_match_for_pssm()
    error_count += check_hiv_seeds(project_config, unchecked_ref_names)
    error_count += check_hiv_coordinates(project_config, unchecked_ref_names)

    print(fill_report(f'Unchecked refs: {", ".join(sorted(unchecked_ref_names))}'))
    error_count += len(unchecked_ref_names)
    print(f'Total errors: {error_count}.')


def check_hiv_seeds(project_config, unchecked_ref_names: set):
    print("""\
HIV seed references come from the HIV Sequence Database Compendium. They are
downloaded from https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html
The HIV1-CON-XX-Consensus-seed is a special case: the consensus of consensuses,
downloaded from the same page with the Consensus/Ancestral alignment type.
""")
    sequences = fetch_alignment_sequences(2004,
                                          'CON',  # Consensus/Ancestral
                                          'ENV')
    consensus = sequences['CON_OF_CONS'].replace('-', '').upper()

    source_sequences = fetch_alignment_sequences('2015', 'COM')
    consensus_accession = 'Consensus'
    assert consensus_accession not in source_sequences, sorted(source_sequences.keys())
    source_sequences[consensus_accession] = consensus
    ref_names = project_config.getProjectSeeds('HIV')
    unchecked_ref_names.difference_update(ref_names)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences,
                                         name_part=3)
    print(report)
    return error_count


def check_hiv_coordinates(project_config, unchecked_ref_names: set):
    print("""\
HIV coordinate references come from HXB2 (K03455). That sequence is trimmed
to each of the gene regions, then translated to amino acids. The gene positions
are documented at https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html
We replace HXB2's premature stop codon in Nef with the codon from Consensus B,
and drop a nucleotide of Vpr to avoid the frame shift.
V3LOOP is a subsection of Env from Consensus B with one nucleotide replaced by
the nucleotide from the overall consensus.
""")
    hxb2 = fetch_by_accession('K03455')

    consensus_sequences = fetch_alignment_sequences(2004,
                                                    'CON',  # Consensus/Ancestral
                                                    'NEF')
    consensus_b_nef = consensus_sequences['CONSENSUS_B'].replace('-', '').upper()

    consensus_sequences = fetch_alignment_sequences('2004',
                                                    'CON',  # Consensus/Ancestral
                                                    'ENV')
    consensus_b_env = consensus_sequences['CONSENSUS_B'].replace('-', '').upper()
    overall_consensus_env = \
        consensus_sequences['CON_OF_CONS'].replace('-', '').upper()

    region_boundaries = {'HIV1B-gag': (789, 2289),
                         'PR': (2252, 2549),
                         'RT': (2549, 3869),  # p51 only
                         'INT': (4229, 5093),
                         'HIV1B-vif': (5040, 5616),
                         'HIV1B-vpr': (5558, 5847),
                         'HIV1B-vpu': (6061, 6307),
                         'GP41': (7757, 8792),
                         'HIV1B-nef': (8796, 9414)}
    source_sequences = {}
    for region, (start, stop) in region_boundaries.items():
        source_nuc_sequence = hxb2[start:stop]
        if region == 'HIV1B-nef':
            # Replace premature stop codon with Consensus B.
            source_nuc_sequence = (source_nuc_sequence[:369] +
                                   consensus_b_nef[369:372] +
                                   source_nuc_sequence[372:])
        elif region == 'HIV1B-vpr':
            # Drop nucleotide to avoid frame shift.
            source_nuc_sequence = (source_nuc_sequence[:213] +
                                   source_nuc_sequence[214:])
        source_sequences[region] = translate(source_nuc_sequence)

    consensus_b_v3_nuc_seq = consensus_b_env[876:981]
    overall_v3_nuc_seq = overall_consensus_env[855:960]
    consensus_b_v3_amino_seq = translate(consensus_b_v3_nuc_seq)
    overall_v3_amino_seq = translate(overall_v3_nuc_seq)
    assert 'CTRPNNNTRKSIHIGPGRAFYTTGEIIGDIRQAHC' == consensus_b_v3_amino_seq, consensus_b_v3_amino_seq
    assert 'CTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHC' == overall_v3_amino_seq, overall_v3_amino_seq
    # Difference is here:        ^

    combined_v3_nuc_seq = (consensus_b_v3_nuc_seq[:63] +
                           overall_v3_nuc_seq[63] +
                           consensus_b_v3_nuc_seq[64:])
    source_sequences['V3LOOP'] = translate(combined_v3_nuc_seq)

    hiv_project = project_config.config['projects']['HIV']
    ref_names = {project_region['coordinate_region']
                 for project_region in hiv_project['regions']}
    unchecked_ref_names.difference_update(ref_names)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


def find_v3_match(env_seq, v3_loop_seq):
    for frame in range(3):
        env_aminos = translate('-' * frame + env_seq)
        aligned_env, aligned_v3, score = align_it_aa(env_aminos,
                                                     v3_loop_seq,
                                                     GAP_OPEN_COST,
                                                     GAP_EXTEND_COST,
                                                     USE_TERMINAL_COST)
        if score < 200:
            continue
        offset = (len(aligned_v3) - len(aligned_v3.lstrip('-'))) * 3 + frame
        end = offset + len(v3_loop_seq)*3
        print(translate(env_seq[offset:end]))
        env_section = env_seq[offset + 60:offset + 69]
        if env_section == 'TATGCAACA':
            print('**G**')
        elif env_section == 'TATACAACA':
            print('**A**')
        else:
            print(f'??{env_section}??')
        print(env_section)
        print(translate(env_section))
        print(score)


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


def fill_report(report):
    return fill(report, break_long_words=False, break_on_hyphens=False) + '\n'


def compare_config(ref_names, project_config, source_sequences, name_part=None):
    differ = Differ()
    error_count = 0
    matching_references = []
    missing_references = []
    diff_report = StringIO()
    for ref_name in sorted(ref_names):
        project_sequence = project_config.getReference(ref_name)
        if name_part is None:
            source_name = ref_name
        else:
            name_parts = ref_name.split('-')
            source_name = name_parts[name_part]
        source_sequence = source_sequences.get(source_name)
        source_sequence = source_sequence and source_sequence.replace('-', '')
        if source_sequence is None:
            error_count += 1
            missing_references.append(ref_name)
        elif source_sequence == project_sequence:
            matching_references.append(ref_name)
        else:
            error_count += 1
            diff_report.write(ref_name + '\n')
            for line in differ.compare([source_sequence+'\n'],
                                       [project_sequence+'\n']):
                diff_report.write(line)
            diff_report.write('\n')

    report = StringIO()
    if matching_references:
        matching_references.sort()
        match_report = f'Matching references: {", ".join(matching_references)}'
        match_report = fill_report(match_report)
        report.write(match_report)
    if missing_references:
        missing_references.sort()
        missing_report = (f'ERROR: references missing from source: '
                          f'{", ".join(missing_references)}')
        missing_report = fill_report(missing_report)
        report.write(missing_report)
    if diff_report.tell():
        report.write('ERROR: changed references:\n')
        report.write(diff_report.getvalue())
    return report.getvalue(), error_count


if __name__ == '__main__':
    main()
