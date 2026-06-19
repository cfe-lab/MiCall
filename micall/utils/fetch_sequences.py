""" Report where references come from, and validate that they match the source.

Running this script reports a source for each group of references, downloads
them again, and checks that they match the source.
Finally, it reports any references that aren't documented yet.
"""

import re
import typing
from collections import defaultdict
from difflib import Differ, SequenceMatcher
from io import StringIO
from itertools import chain, groupby, zip_longest
from textwrap import fill, wrap

from yaml import safe_load

from micall.data.landmark_reader import LandmarkReader
from micall.resistance.asi_algorithm import WILD_TYPES_PATH
from micall.utils.translation import translate

try:
    # noinspection PyPackageRequirements
    import requests
except ImportError:
    # Allow tests to run without requests module
    requests = None

from micall.core.project_config import ProjectConfig


def main():
    project_config = ProjectConfig.loadDefault()
    scoring_config = ProjectConfig.loadScoring()
    error_count = 0
    unchecked_ref_names = set(project_config.getAllReferences().keys())
    unchecked_ref_names.update(set(project_config.getAllGenotypeReferences().keys()))
    error_count += check_hcv_seeds(project_config, unchecked_ref_names)
    error_count += check_hcv_coordinates(project_config, unchecked_ref_names)
    error_count += check_hiv_seeds(project_config, unchecked_ref_names)
    error_count += check_hiv_coordinates(project_config, unchecked_ref_names)
    error_count += check_hiv_wild_types(project_config)
    error_count += check_hla_seeds(project_config, unchecked_ref_names)
    error_count += check_hla_coordinates(project_config, unchecked_ref_names)
    error_count += check_sars_seeds(project_config, unchecked_ref_names)
    error_count += check_sars_coordinates(project_config, unchecked_ref_names)
    error_count += check_scoring_config(scoring_config, project_config)

    if not unchecked_ref_names:
        print('No unchecked refs.')
    else:
        print(fill_report(f'Unchecked refs: '
                          f'{", ".join(sorted(unchecked_ref_names))}'))
        error_count += len(unchecked_ref_names)
    print(f'Total errors: {error_count}.')


def check_length_and_translate(ref_name, sequence):
    length_error = 0
    if len(sequence) % 3 != 0:
        print(f"Warning: region {ref_name} has length {len(sequence)}, "
              f"not an even multiple of 3.")
        length_error = 1
    translated_seq = translate(sequence)
    return translated_seq, length_error


def check_hiv_seeds(project_config, unchecked_ref_names: set):
    print("""\
HIV seed references come from the HIV Sequence Database Compendium. They are
downloaded from https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html
Recombinant references are excluded (names start with a digit).
The HIV1-CON-XX-Consensus-seed is a special case: the consensus of consensuses,
downloaded from the same page with the Consensus/Ancestral alignment type.
See the project_seeds_from_compendium.py script for full details.
More special cases were added for the HIVGHA project to include circulating
recombinant forms common in Ghana. The CRF references are excluded in the
pipeline for all but the HIVGHA project code.
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

    ghana_source_names = ["HIV1-CRF02_AG-GH-AB286855-seed",
                          "HIV1-CRF06_CPX-GH-AB286851-seed"]
    for source_name in ghana_source_names:
        parts = source_name.split('-')
        accession_number = parts[3]
        assert source_name not in source_sequences, sorted(source_sequences.keys())
        source_sequences[accession_number] = fetch_by_accession(accession_number)

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

    # Fetch GAG sequences for CA region
    gag_sequences = fetch_alignment_sequences(2004,
                                              'CON',  # Consensus/Ancestral
                                              'GAG')
    consensus_b_gag = gag_sequences['CONSENSUS_B'].upper()

    with open(WILD_TYPES_PATH) as wild_types_file:
        wild_types = safe_load(wild_types_file)
    boundaries = {'PR': (171, 468),
                  'RT': (468, 2148),
                  'INT': (2148, 3017)}
    gag_boundaries = {'CA': (342, 1089)}  # CA/p24 region in gag
    ref_names = sorted(boundaries.keys()) + sorted(gag_boundaries.keys())
    source_wild_types = {}
    length_errors = 0
    for ref_name, (start, end) in boundaries.items():
        source_nuc_sequence = consensus_b[start:end]
        source_wild_type, length_error = check_length_and_translate(ref_name, source_nuc_sequence)
        length_errors += length_error
        source_wild_types[ref_name] = source_wild_type
        mapping_ref = project_config.getReference(ref_name)
        wild_type_length = len(source_wild_type)
        mapping_length = len(mapping_ref)
        if wild_type_length != mapping_length:
            print(f'{ref_name} expected length {mapping_length}, '
                  f'was {wild_type_length}')
            length_errors += 1
    for ref_name, (start, end) in gag_boundaries.items():
        source_nuc_sequence = consensus_b_gag[start:end]
        source_wild_type, length_error = check_length_and_translate(ref_name, source_nuc_sequence)
        length_errors += length_error
        source_wild_types[ref_name] = source_wild_type
        mapping_ref = project_config.getReference(ref_name)
        wild_type_length = len(source_wild_type)
        mapping_length = len(mapping_ref)
        if wild_type_length != mapping_length:
            print(f'{ref_name} expected length {mapping_length}, '
                  f'was {wild_type_length}')
            length_errors += 1
    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_wild_types,
                                         reference_overrides=wild_types)
    print(report)
    return error_count + length_errors


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

    region_names = ['HIV1B-gag',
                    'CA',
                    'PR',
                    'RT',  # p51 only
                    'INT',
                    'HIV1B-tat',
                    'HIV1B-vif',
                    'HIV1B-vpr',
                    'HIV1B-vpu',
                    'GP120',
                    'GP41',
                    'HIV1B-nef']
    nucleotide_region_names = ["HIV1B-5' LTR",
                               'HIV1B-635-789',
                               'HIV1B-sl1',
                               'HIV1B-sl2',
                               'HIV1B-sl3',
                               'HIV1B-sl4',
                               'HIV1B-D4',
                               'HIV1B-6046-6061',
                               'HIV1B-8796',
                               "HIV1B-3' LTR"]
    source_sequences = {}
    ref_positions = set()
    landmark_reader = LandmarkReader.load()
    length_errors = 0
    for region_name in region_names + nucleotide_region_names:
        source_nuc_sequence = extract_region(landmark_reader,
                                             'HIV1-B-FR-K03455-seed',
                                             hxb2,
                                             region_name,
                                             ref_positions)
        if region_name == 'HIV1B-nef':
            # Replace premature stop codon with Consensus B.
            source_nuc_sequence = (source_nuc_sequence[:369] +
                                   consensus_b_nef[369:372] +
                                   source_nuc_sequence[372:])
        if region_name in nucleotide_region_names:
            source_sequences[region_name] = source_nuc_sequence
        else:
            source_sequences[region_name], length_error = check_length_and_translate(
                region_name,
                source_nuc_sequence)
            length_errors += length_error

    consensus_b_v3_nuc_seq = consensus_b_env[876:981]
    overall_v3_nuc_seq = overall_consensus_env[855:960]
    consensus_b_v3_amino_seq, length_error = check_length_and_translate(
        'consensus_b_v3_nuc_seq',
        consensus_b_v3_nuc_seq)
    length_errors += length_error
    overall_v3_amino_seq, length_error = check_length_and_translate(
        'overall_v3_nuc_seq',
        overall_v3_nuc_seq)
    length_errors += length_error
    assert 'CTRPNNNTRKSIHIGPGRAFYTTGEIIGDIRQAHC' == consensus_b_v3_amino_seq, consensus_b_v3_amino_seq
    assert 'CTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHC' == overall_v3_amino_seq, overall_v3_amino_seq
    # Difference is here:        ^

    combined_v3_nuc_seq = (consensus_b_v3_nuc_seq[:63] +
                           overall_v3_nuc_seq[63] +
                           consensus_b_v3_nuc_seq[64:])
    source_sequences['V3LOOP'], length_error = check_length_and_translate('V3LOOP', combined_v3_nuc_seq)
    length_errors += length_error

    hiv_project = project_config.config['projects']['HIV']
    ref_names = {project_region['coordinate_region']
                 for project_region in hiv_project['regions']}
    unchecked_ref_names.difference_update(ref_names)

    report_missing_positions(ref_positions, hxb2, 'HIV')
    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count + length_errors


def extract_region(landmark_reader: LandmarkReader,
                   coordinate_name: str,
                   coordinate_seq: str,
                   region_name: str,
                   ref_positions: set = None) -> str:
    """ Extract a section of reference sequence, based on a gene landmark.

    :param landmark_reader: source of landmark definitions
    :param coordinate_name: name of coordinate reference for the landmarks
    :param coordinate_seq: coordinate sequence to extract from
    :param region_name: gene region to extract
    :param ref_positions: track which positions have been extracted, so we can
        check for gaps.
    :return: extracted sequence for the gene region
    """
    region = landmark_reader.get_gene(coordinate_name,
                                      region_name)
    full_region = landmark_reader.get_gene(coordinate_name,
                                           region_name,
                                           drop_stop_codon=False)
    check_stop_codons = True
    if check_stop_codons:
        region = full_region
    start = region['start']
    end = region['end']
    source_nuc_sequence = coordinate_seq[start - 1:end]
    duplicated_base = full_region.get('duplicated_pos')
    skipped_base = full_region.get('skipped_pos')
    if ref_positions is not None:
        ref_positions.update(range(full_region['start'], full_region['end'] + 1))
    if duplicated_base is not None:
        source_nuc_sequence = (
                source_nuc_sequence[:duplicated_base - start + 1] +
                source_nuc_sequence[duplicated_base - start:])
    if skipped_base is not None:
        source_nuc_sequence = (source_nuc_sequence[:skipped_base - start] +
                               source_nuc_sequence[skipped_base - start + 1:])
    return source_nuc_sequence


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
                  # EF108306.2 is available, but only extends 5' and 3'.
                  'HCV-7a': 'Ref.7a.CA.x.QC69.EF108306.1'}
    source_sequences = {}
    for ref_name, source_id in source_ids.items():
        match = re.search(r'[A-Z_]+[0-9.]+$', source_id)
        assert match is not None, source_id
        accession_number = match.group(0)
        source_sequences[ref_name] = fetch_by_accession(accession_number)

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
                         # EF108306.2 is available, but only extends 5' and 3'.
                         'HCV7': 'EF108306.1'}
    source_nuc_sequences = {
        genotype: fetch_by_accession(accession_number)
        for genotype, accession_number in accession_numbers.items()}

    # Boundary positions in landmarks are from the European HCV database records.
    # https://euhcvdb.ibcp.fr/euHCVdb/do/displayHCVEntry?primaryAC=AF009606
    # That is the original H77 accession number for HCV1A. NC_004102 is the
    # curated and annotated version that was derived from the AF009606 entry.
    # All the other genotypes can be found by their regular accession numbers.

    report, errors_genotype = compare_config(accession_numbers.keys(),
                                             project_config,
                                             source_nuc_sequences)
    unchecked_ref_names.difference_update(accession_numbers.keys())
    print(report)

    hcv_project = project_config.config['projects']['HCV']
    ref_names = {project_region['coordinate_region']
                 for project_region in hcv_project['regions']}
    unchecked_ref_names.difference_update(ref_names)
    landmark_reader = LandmarkReader.load()

    source_sequences = {}
    genotype_ref_positions = defaultdict(set)
    length_errors = 0
    for ref_name in sorted(ref_names):
        ref_parts = ref_name.split('-')
        genotype = ref_parts[0]
        seed_name = f'HCV-{genotype[3:].lower()}'
        if len(seed_name) == 5:
            seed_name += 'a'

        coordinate_seq = source_nuc_sequences[genotype]
        ref_positions = genotype_ref_positions[genotype]
        nuc_seq_ref_trimmed = extract_region(landmark_reader,
                                             genotype,
                                             coordinate_seq,
                                             ref_name,
                                             ref_positions)
        source_sequences[ref_name], length_error = check_length_and_translate(
            ref_name,
            nuc_seq_ref_trimmed)
        length_errors += length_error

    # for genotype, ref_positions in genotype_ref_positions.items():
    #     seed_sequence = source_nuc_sequences[genotype]
    #     report_missing_positions(ref_positions, seed_sequence, genotype)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count + length_errors + errors_genotype


def check_hla_seeds(project_config, unchecked_ref_names: set):
    print("""\
HLA seed reference is downloaded from accession AJ458991.3
""")

    source_sequence = fetch_by_accession('AJ458991.3')
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
HLA coordinate references are translated from the seed reference. Coordinates
come from the features described at https://www.ncbi.nlm.nih.gov/nuccore/AJ458991
""")
    ref_names = ('HLA-B-exon2', 'HLA-B-exon3')

    seed_sequence = project_config.getReference('HLA-B-seed')
    unchecked_ref_names.difference_update(ref_names)
    landmark_reader = LandmarkReader.load()

    source_sequences = {}
    length_errors = 0
    for ref_name in ref_names:
        source_nuc_sequence = extract_region(landmark_reader,
                                             'HLA-B-seed',
                                             seed_sequence,
                                             ref_name)
        source_sequences[ref_name], length_error = check_length_and_translate(ref_name, source_nuc_sequence)
        length_errors += length_error

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count + length_errors


def check_sars_seeds(project_config, unchecked_ref_names: set):
    print("""\
SARS-CoV-2 seed reference is downloaded from accession MN908947
""")

    source_sequence = fetch_by_accession('MN908947')
    source_sequences = {'SARS-CoV-2-seed': source_sequence}

    ref_names = project_config.getProjectSeeds('SARSCOV2')
    unchecked_ref_names.difference_update(ref_names)

    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count


def check_sars_coordinates(project_config, unchecked_ref_names: set):
    print("""\
SARS-CoV-2 coordinate references are translated from the seed reference.
Some boundaries come from https://www.ncbi.nlm.nih.gov/nuccore/MN908947
Others come from supplemental data #2 in:
https://www.nature.com/articles/s41467-021-22905-7#Sec26
TRS-B entries, or Transcription Regulatory Sequences (Body) were created to
cover any sections that weren't included elsewhere.
""")
    ref_names = ["SARS-CoV-2-5'UTR",
                 'SARS-CoV-2-ORF1ab',
                 'SARS-CoV-2-S',
                 'SARS-CoV-2-ORF3a',
                 'SARS-CoV-2-E',
                 'SARS-CoV-2-M',
                 'SARS-CoV-2-ORF6',
                 'SARS-CoV-2-ORF7a',
                 'SARS-CoV-2-ORF7b',
                 'SARS-CoV-2-ORF8',
                 'SARS-CoV-2-N',
                 'SARS-CoV-2-ORF10',
                 "SARS-CoV-2-3'UTR",
                 'SARS-CoV-2-nsp1',
                 'SARS-CoV-2-nsp2',
                 'SARS-CoV-2-nsp3',
                 'SARS-CoV-2-nsp4',
                 'SARS-CoV-2-nsp5',
                 'SARS-CoV-2-nsp6',
                 'SARS-CoV-2-nsp7',
                 'SARS-CoV-2-nsp8',
                 'SARS-CoV-2-nsp9',
                 'SARS-CoV-2-nsp10',
                 'SARS-CoV-2-nsp11',
                 'SARS-CoV-2-nsp12',
                 'SARS-CoV-2-nsp13',
                 'SARS-CoV-2-nsp14',
                 'SARS-CoV-2-nsp15',
                 'SARS-CoV-2-nsp16',
                 'SARS-CoV-2-TRS-B-1',
                 'SARS-CoV-2-TRS-B-2',
                 'SARS-CoV-2-TRS-B-3',
                 'SARS-CoV-2-TRS-B-4',
                 'SARS-CoV-2-TRS-B-5',
                 'SARS-CoV-2-TRS-B-6',
                 'SARS-CoV-2-TRS-B-7',
                 'SARS-CoV-2-TRS-B-8',
                 'SARS-CoV-2-TRS-B-9']

    seed_name = 'SARS-CoV-2-seed'
    seed_sequence = project_config.getReference(seed_name)
    coordinate_references = project_config.getCoordinateReferences(seed_name)
    unconfigured_references = set(ref_names).difference(coordinate_references)
    if unconfigured_references:
        print('Unconfigured references:')
        print('  ' + ', '.join(sorted(unconfigured_references)))
    unchecked_ref_names.difference_update(ref_names)
    landmark_reader = LandmarkReader.load()

    source_sequences = {}
    ref_positions = set()
    length_errors = 0
    for ref_name in ref_names:
        source_nuc_sequence = extract_region(landmark_reader,
                                             seed_name,
                                             seed_sequence,
                                             ref_name,
                                             ref_positions,)
        if 'TRS-B' in ref_name or ref_name.endswith('UTR'):
            source_sequences[ref_name] = source_nuc_sequence
        else:
            source_sequences[ref_name], length_error = check_length_and_translate(ref_name, source_nuc_sequence)
            length_errors += length_error
        print(ref_name, len(source_sequences[ref_name]))

    report_missing_positions(ref_positions, seed_sequence, 'SARS-COV-2')
    report, error_count = compare_config(ref_names,
                                         project_config,
                                         source_sequences)
    print(report)
    return error_count + length_errors


def check_scoring_config(scoring_config: ProjectConfig,
                         project_config: ProjectConfig) -> int:
    error_count = 0
    protein_lengths = {}
    for ref_name in project_config.getAllReferences():
        if not project_config.isAmino(ref_name):
            continue
        ref = project_config.getReference(ref_name)
        protein_lengths[ref_name] = len(ref)
    unscored_protein_names = set(protein_lengths)
    for project_name, project in scoring_config.config['projects'].items():
        for region in project['regions']:
            coordinate_name = region['coordinate_region']
            coordinate_length = region['coordinate_region_length']
            try:
                expected_length = protein_lengths[coordinate_name]
            except KeyError:
                error_count += 1
                print(f'{coordinate_name} not a protein in {project_name}?')
                continue
            if coordinate_length != expected_length:
                error_count += 1
                print(f'{coordinate_name} expected length {expected_length} '
                      f'in {project_name} but found {coordinate_length}.')
            unscored_protein_names.discard(coordinate_name)

    if unscored_protein_names:
        print('Unscored protein regions:')
        print(', '.join(sorted(unscored_protein_names)))

    return error_count


def report_missing_positions(ref_positions: typing.Set[int],
                             seed_sequence: str,
                             genome_name: str):
    """ Print positions that were not covered by reported regions.

    :param ref_positions: all 1-based positions that were covered
    :param seed_sequence: reference sequence to show how long it is.
    :param genome_name: label for displaying errors
    """
    all_positions = set(range(1, len(seed_sequence) + 1))
    print_positions(f'Missing positions from {genome_name} genome:',
                    sorted(all_positions - ref_positions))

    print_positions(f'Extra positions from {genome_name} genome:',
                    sorted(ref_positions - all_positions))


def print_positions(label, positions):
    if not positions:
        return
    print(label)
    for offset, items in groupby(enumerate(positions),
                                 lambda item: item[1] - item[0]):
        positions = [x for i, x in items]
        if len(positions) == 1:
            print(*positions)
        else:
            print(f'{positions[0]}-{positions[-1]}')


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


def fetch_by_accession(accession):
    response = requests.get('https://www.ncbi.nlm.nih.gov/nuccore/' +
                            accession + '?report=fasta&format=text')
    response.raise_for_status()
    match = re.search(r'<div id="viewercontent1".*val="(\d+)"', response.text)
    if not match:
        raise ValueError(f'Accession not found: {accession!r}.')
    response_id = match.group(1)

    response = requests.get(f'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?'
                            f'id={response_id}&db=nuccore&report=fasta')
    response.raise_for_status()
    fasta = StringIO(response.text)
    source_sequence, = parse_fasta(fasta).values()
    return source_sequence


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
            try:
                project_sequence = project_config.getReference(ref_name)
            except KeyError:
                try:
                    project_sequence = project_config.getGenotypeReference(ref_name)
                except KeyError:
                    project_sequence = ''
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
        is_compressed = False
        for source_name, source_sequence in sorted(source_sequences.items()):
            report.write(source_name + '\n')
            if is_compressed:
                if 75 < len(source_sequence):
                    source_sequence = f'{source_sequence[:35]}[...]' \
                                      f'{source_sequence[-35:]}'
                report.write('  ' + source_sequence + '\n')
            else:
                iters = [iter(source_sequence)] * 65
                lines = []
                for row in zip_longest(*iters, fillvalue=''):
                    line = ''.join(row)
                    lines.append(f'"{line}"')
                report.write(',\n'.join(lines) + '\n')
        error_count += len(source_sequences)

    if diff_report.tell():
        report.write('ERROR: changed references:\n')
        report.write(diff_report.getvalue())
    return report.getvalue(), error_count


def join_and_wrap(sections: typing.Sequence[str]) -> typing.List[str]:
    joined = ''.join(sections)
    wrapped = wrap(joined, width=65)
    return [line + '\n' for line in wrapped]


if __name__ == '__main__':
    main()
