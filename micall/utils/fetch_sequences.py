""" Report where references come from, and validate that they match the source.

Running this script reports a source for each group of references, downloads
them again, and checks that they match the source.
Finally, it reports any references that aren't documented yet.
"""

import re
import typing
from difflib import Differ, SequenceMatcher
from io import StringIO
from itertools import chain
from textwrap import fill, wrap

from yaml import safe_load

from micall.resistance.asi_algorithm import WILD_TYPES_PATH
from micall.utils.translation import translate

try:
    import requests
except ImportError:
    # Allow tests to run without requests module
    requests = None

from micall.core.project_config import ProjectConfig


def main():
    project_config = ProjectConfig.loadDefault()
    error_count = 0
    unchecked_ref_names = set(project_config.getAllReferences().keys())
    error_count += check_hcv_seeds(project_config, unchecked_ref_names)
    error_count += check_hcv_coordinates(project_config, unchecked_ref_names)
    error_count += check_hiv_seeds(project_config, unchecked_ref_names)
    error_count += check_hiv_coordinates(project_config, unchecked_ref_names)
    error_count += check_hiv_wild_types(project_config)
    error_count += check_hla_seeds(project_config, unchecked_ref_names)
    error_count += check_hla_coordinates(project_config, unchecked_ref_names)

    if not unchecked_ref_names:
        print('No unchecked refs.')
    else:
        print(fill_report(f'Unchecked refs: '
                          f'{", ".join(sorted(unchecked_ref_names))}'))
        error_count += len(unchecked_ref_names)
    print(f'Total errors: {error_count}.')


def check_hiv_seeds(project_config, unchecked_ref_names: set):
    print("""\
HIV seed references come from the HIV Sequence Database Compendium. They are
downloaded from https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html
Recombinant references are excluded (names start with a digit).
The HIV1-CON-XX-Consensus-seed is a special case: the consensus of consensuses,
downloaded from the same page with the Consensus/Ancestral alignment type.
See the project_seeds_from_compendium.py script for full details.
""")
    sequences = fetch_alignment_sequences(2004,
                                          'CON',  # Consensus/Ancestral
                                          'ENV')
    consensus = sequences['CON_OF_CONS'].replace('-', '').upper()

    compendium_sequences = fetch_alignment_sequences('2015', 'COM')
    source_sequences = {ref_name.split('.')[-1]: ref_sequence
                        for ref_name, ref_sequence in compendium_sequences.items()
                        if not ref_name[0].isdigit()}  # Exclude recombinants.
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


def check_hiv_wild_types(project_config):
    print("""\
HIV wild types for resistance reports are extracted from Consensus B.
""")
    sequences = fetch_alignment_sequences(2004,
                                          'CON',  # Consensus/Ancestral
                                          'POL')
    consensus_b = sequences['CONSENSUS_B'].upper()

    with open(WILD_TYPES_PATH) as wild_types_file:
        wild_types = safe_load(wild_types_file)
    boundaries = {'PR': (171, 468),
                  'RT': (468, 1788),
                  'INT': (2148, 3014)}
    ref_names = sorted(boundaries.keys())
    source_wild_types = {}
    for ref_name, (start, end) in boundaries.items():
        source_nuc_sequence = consensus_b[start:end]
        source_wild_types[ref_name] = translate(source_nuc_sequence)
    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_wild_types,
                                         reference_overrides=wild_types)
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
Our G2P reference is the same V3LOOP with a bunch of X's added to leave space
for insertions in the scoring matrix. See pssm_lib.py for more details.
""")
    hxb2 = fetch_hiv_by_accession('K03455')

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
                         'GP120': (6224, 7757),
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


def check_hcv_seeds(project_config, unchecked_ref_names: set):
    print("""\
Most HCV seed references come from the HCV Sequence Genotype Reference,
downloaded from https://hcv.lanl.gov/content/sequence/NEWALIGN/align.html
This script contains a complete list of the reference accession numbers.
""")
    source_ids = {'HCV-1a': 'Ref.1a.US.77.H77.NC_004102',
                  'HCV-1b': 'Ref.1b.BR.03.BR1427_P1_10-7-03.EF032892',
                  'HCV-1c': 'Ref.1c.ID.x.HC-G9.D14853',
                  'HCV-1e': 'KJ439769',
                  'HCV-1g': 'Ref.1g.ES.x.1804.AM910652',
                  'HCV-2a': 'Ref.2a.JP.x.AY746460.AY746460',
                  'HCV-2b': 'Ref.2b.JP.x.HC-J8.D10988',
                  'HCV-2c': 'Ref.2c.x.x.BEBE1.D50409',
                  'HCV-2i': 'Ref.2i.VN.x.D54.DQ155561',
                  'HCV-2j': 'Ref.2j.VE.05.C1292.HM777359',
                  'HCV-2k': 'Ref.2k.MD.x.VAT96.AB031663',
                  'HCV-2q': 'Ref.2q.ES.01.852.FN666429',
                  'HCV-3a': 'Ref.3a.DE.x.HCVCENS1.X76918',
                  'HCV-3b': 'Ref.3b.JP.x.HCV-Tr.D49374',
                  'HCV-3i': 'Ref.3i.IN.02.IND-HCV-3i.FJ407092',
                  'HCV-3k': 'Ref.3k.ID.x.JK049.D63821',
                  'HCV-4a': 'Ref.4a.EG.x.ED43.Y11604',
                  'HCV-4b': 'Ref.4b.CA.x.QC264.FJ462435',
                  'HCV-4c': 'Ref.4c.CA.x.QC381.FJ462436',
                  'HCV-4d': 'Ref.4d.CA.x.QC382.FJ462437',
                  'HCV-4f': 'Ref.4f.FR.x.IFBT84.EF589160',
                  'HCV-4g': 'Ref.4g.CA.x.QC193.FJ462432',
                  'HCV-4k': 'Ref.4k.CA.x.QC383.FJ462438',
                  'HCV-4l': 'Ref.4l.CA.x.QC274.FJ839870',
                  'HCV-4m': 'Ref.4m.CA.x.QC249.FJ462433',
                  'HCV-4n': 'Ref.4n.CA.x.QC97.FJ462441',
                  'HCV-4o': 'Ref.4o.CA.x.QC93.FJ462440',
                  'HCV-4p': 'Ref.4p.CA.x.QC139.FJ462431',
                  'HCV-4q': 'Ref.4q.CA.x.QC262.FJ462434',
                  'HCV-4r': 'Ref.4r.CA.x.QC384.FJ462439',
                  'HCV-4t': 'Ref.4t.CA.x.QC155.FJ839869',
                  'HCV-4v': 'Ref.4v.CY.05.CYHCV048.HQ537008',
                  'HCV-5a': 'Ref.5a.GB.x.EUH1480.NC_009826',
                  'HCV-6a': 'Ref.6a.HK.x.6a33.AY859526',
                  'HCV-6b': 'Ref.6b.x.x.Th580.NC_009827',
                  'HCV-6c': 'Ref.6c.TH.x.Th846.EF424629',
                  'HCV-6d': 'Ref.6d.VN.x.VN235.D84263',
                  'HCV-6e': 'Ref.6e.CN.x.GX004.DQ314805',
                  'HCV-6f': 'Ref.6f.TH.x.C-0044.DQ835760',
                  'HCV-6g': 'Ref.6g.HK.x.HK6554.DQ314806',
                  'HCV-6h': 'Ref.6h.VN.x.VN004.D84265',
                  'HCV-6i': 'Ref.6i.TH.x.C-0159.DQ835762',
                  'HCV-6j': 'Ref.6j.TH.x.C-0667.DQ835761',
                  'HCV-6k': 'Ref.6k.CN.x.KM41.DQ278893',
                  'HCV-6l': 'Ref.6l.US.x.537796.EF424628',
                  'HCV-6m': 'Ref.6m.TH.x.C-0185.DQ835765',
                  'HCV-6n': 'Ref.6n.CN.x.KM42.DQ278894',
                  'HCV-6o': 'Ref.6o.CA.x.QC227.EF424627',
                  'HCV-6p': 'Ref.6p.CA.x.QC216.EF424626',
                  'HCV-6q': 'Ref.6q.CA.x.QC99.EF424625',
                  'HCV-6r': 'Ref.6r.CA.x.QC245.EU408328',
                  'HCV-6s': 'Ref.6s.CA.x.QC66.EU408329',
                  'HCV-6t': 'Ref.6t.VN.x.D49.EU246939',
                  'HCV-6u': 'EU408330',
                  'HCV-6v': 'EU158186',
                  'HCV-6w': 'Ref.6w.TW.x.HCV-6-D140.EU643834',
                  'HCV-7a': 'Ref.7a.CA.x.QC69.EF108306'}
    source_sequences = {}
    for ref_name, source_id in source_ids.items():
        accession_number = source_id.split('.')[-1]
        source_sequences[ref_name] = fetch_hcv_by_accession(accession_number)

    ref_names = project_config.getProjectSeeds('HCV')
    unchecked_ref_names.difference_update(ref_names)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


def check_hcv_coordinates(project_config, unchecked_ref_names: set):
    print("""\
Most HCV coordinate references were listed in the FDA guidance:
https://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM340712.pdf
This script contains a complete list of the reference accession numbers.
""")
    accession_numbers = {'HCV1A': 'NC_004102',
                         'HCV1B': 'AJ238799',
                         'HCV2': 'AB047639',
                         'HCV3': 'GU814263',
                         'HCV4': 'GU814265',
                         'HCV5': 'AF064490',
                         'HCV6': 'Y12083',
                         'HCV7': 'EF108306'}
    source_nuc_sequences = {
        genotype: fetch_hcv_by_accession(accession_number)
        for genotype, accession_number in accession_numbers.items()}

    gene_names = [
        'Core', 'E1', 'E2', 'p7', 'NS2', 'NS3', 'NS4a', 'NS4b', 'NS5a', 'NS5b']
    # Boundary positions are from the European HCV database records.
    # https://euhcvdb.ibcp.fr/euHCVdb/do/displayHCVEntry?primaryAC=AF009606
    # That is the original H77 accession number for HCV1A. NC_004102 is the
    # curated and annotated version that was derived from the AF009606 entry.
    # All the other genotypes can be found by their regular accession numbers.

    genotype_boundaries = {
        #         Core E1   E2    p7    NS2   NS3   NS4a  NS4b  NS5a  NS5b
        'HCV1A': [342, 915, 1491, 2580, 2769, 3420, 5313, 5475, 6258, 7602, 9375],
        'HCV1B': [342, 915, 1491, 2580, 2769, 3420, 5313, 5475, 6258, 7599, 9372],
        'HCV2': [341, 914, 1490, 2591, 2780, 3431, 5324, 5486, 6269, 7667, 9440],
        'HCV3': [340, 913, 1489, 2596, 2785, 3436, 5329, 5491, 6274, 7630, 9403],
        'HCV4': [341, 914, 1490, 2579, 2768, 3419, 5312, 5474, 6257, 7592, 9365],
        'HCV5': [247, 820, 1396, 2488, 2677, 3328, 5221, 5383, 6166, 7516, 9289],
        'HCV6': [284, 857, 1433, 2534, 2723, 3374, 5267, 5429, 6212, 7565, 9338],
        'HCV7': [309, 882, 1458, 2559, 2748, 3399, 5292, 5454, 6237, 7575, 9348]}

    hcv_project = project_config.config['projects']['HCV']
    ref_names = {project_region['coordinate_region']
                 for project_region in hcv_project['regions']}
    unchecked_ref_names.difference_update(ref_names)

    source_sequences = {}
    for ref_name in sorted(ref_names):
        ref_parts = ref_name.split('-')
        genotype = ref_parts[0]
        gene_name = ref_parts[-1]
        gene_index = gene_names.index(gene_name)
        boundaries = genotype_boundaries[genotype]
        start, stop = boundaries[gene_index:gene_index+2]
        nuc_seq_ref_trimmed = source_nuc_sequences[genotype][start-1:stop-1]
        source_sequences[ref_name] = translate(nuc_seq_ref_trimmed)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


def check_hla_seeds(project_config, unchecked_ref_names: set):
    print("""\
HLA seed reference is downloaded from
https://www.ncbi.nlm.nih.gov/nuccore/AJ458991.3?report=fasta
""")

    response = requests.get('https://www.ncbi.nlm.nih.gov/nuccore/AJ458991.3?'
                            'report=fasta&format=text')
    response.raise_for_status()
    match = re.search(r'<div id="viewercontent1".*val="(\d+)"', response.text)
    response_id = match.group(1)

    response = requests.get(f'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?'
                            f'id={response_id}&db=nuccore&report=fasta')
    response.raise_for_status()
    hla_fasta = StringIO(response.text)
    source_sequence, = parse_fasta(hla_fasta).values()
    source_sequences = {'HLA-B-seed': source_sequence[301:1869]}

    ref_names = project_config.getProjectSeeds('HLA-B')
    unchecked_ref_names.difference_update(ref_names)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


def check_hla_coordinates(project_config, unchecked_ref_names: set):
    print("""\
HLA coordinate references are translated from the seed reference.
""")
    boundaries = {'HLA-B-exon2': (200, 470),
                  'HLA-B-exon3': (716, 992)}

    seed_sequence = project_config.getReference('HLA-B-seed')
    ref_names = sorted(boundaries.keys())
    unchecked_ref_names.difference_update(ref_names)

    source_sequences = {}
    for ref_name, (start, end) in boundaries.items():
        source_nuc_sequence = seed_sequence[start:end]
        source_sequences[ref_name] = translate(source_nuc_sequence)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


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
    sequences = parse_fasta(StringIO(response.text))
    return sequences


def fetch_hcv_sequences(year, alignment_type='REF', region='GENOME'):
    response = requests.post(
        'https://hcv.lanl.gov/cgi-bin/NEWALIGN/align.cgi',
        data={'ORGANISM': 'HCV',
              'ALIGN_TYPE': alignment_type,
              'SUBORGANISM': 'HCV',
              'PRE_USER': 'predefined',
              'REGION': region,
              'START': '',
              'END': '',
              'GENO_SUB': 'ALL',
              'BASETYPE': 'DNA',
              'YEAR': year,
              'FORMAT': 'fasta',
              'submit': 'Get Alignment'})
    match = re.search(r'<a href="([^"]*)">Download', response.text)
    assert match
    response = requests.get('https://www.hcv.lanl.gov' + match.group(1))
    sequences = parse_fasta(StringIO(response.text))
    return sequences


def fetch_hiv_by_accession(accession):
    response = requests.get(
        'https://www.hiv.lanl.gov/components/sequence/HIV/asearch/'
        'query_one.comp?se_id=' + accession)
    return extract_sequence(response.text)


def fetch_hcv_by_accession(accession):
    response = requests.get(
        'https://www.hcv.lanl.gov/components/sequence/HCV/asearch/'
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


def parse_fasta(fasta):
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
                key = line[1:].strip()

    return sequences


def fill_report(report):
    return fill(report, break_long_words=False, break_on_hyphens=False) + '\n'


def compare_config(ref_names,
                   project_config,
                   source_sequences,
                   name_part=None,
                   reference_overrides=None):
    source_sequences = dict(source_sequences)  # Make a copy.
    source_count = len(source_sequences)
    differ = Differ()
    error_count = 0
    matching_references = []
    missing_references = []
    diff_report = StringIO()
    for ref_name in sorted(ref_names):
        if reference_overrides:
            project_sequence = reference_overrides.get(ref_name)
        else:
            project_sequence = None
        if project_sequence is None:
            project_sequence = project_config.getReference(ref_name)
        if name_part is None:
            source_name = ref_name
        else:
            name_parts = ref_name.split('-')
            source_name = name_parts[name_part]
        try:
            source_sequence = source_sequences.pop(source_name)
        except KeyError:
            source_sequence = None
        source_sequence = source_sequence and source_sequence.replace('-', '')
        if source_sequence is None:
            error_count += 1
            missing_references.append(ref_name)
        elif source_sequence == project_sequence:
            matching_references.append(ref_name)
        else:
            error_count += 1
            diff_report.write(ref_name + '\n')
            source_sections = []
            project_sections = []
            matcher = SequenceMatcher(a=source_sequence,
                                      b=project_sequence,
                                      autojunk=False)
            for tag, i1, i2, j1, j2 in matcher.get_opcodes():
                if tag != 'equal' or i2 - i1 < 210:
                    source_sections.append(source_sequence[i1:i2])
                    project_sections.append(project_sequence[j1:j2])
                else:
                    section = source_sequence[i1:i2]
                    section_size = i2-i1
                    hidden_size = section_size - 140
                    label = f'[{hidden_size} matching characters]'
                    source_sections.append(section[:70])
                    source_sections.append(label)
                    source_sections.append(section[-70:])
                    project_sections.append(section[:70])
                    project_sections.append(label)
                    project_sections.append(section[-70:])
            source_lines = join_and_wrap(source_sections)
            project_lines = join_and_wrap(project_sections)
            for line in differ.compare(source_lines, project_lines):
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
    if source_sequences:
        unused_count = len(source_sequences)
        unused_report = f'ERROR: Left {unused_count} out of {source_count} sequences unused: '
        unused_report += ', '.join(sorted(source_sequences))
        unused_report = fill_report(unused_report)
        report.write(unused_report)
        error_count += len(source_sequences)

    if diff_report.tell():
        report.write('ERROR: changed references:\n')
        report.write(diff_report.getvalue())
    return report.getvalue(), error_count


def join_and_wrap(sections: typing.Sequence[str]) -> typing.List[str]:
    joined = ''.join(sections)
    wrapped = wrap(joined)
    return [line + '\n' for line in wrapped]


if __name__ == '__main__':
    main()
