from collections import namedtuple
from gotoh import align_it_aa
from random import randrange

from micall.core.project_config import ProjectConfig
from micall.utils.translation import translate, reverse_and_complement

FastqSection = namedtuple('FastqSection', 'coord_name start_pos end_pos count mutations')
FastqSection.__new__.__defaults__ = (None, None, None, 1, ())
CodonMutation = namedtuple('CodonMutation', 'pos codon')
FastqFile = namedtuple('FastqFile', 'name extract_num is_reversed sections')


# noinspection PyArgumentList
def main():
    projects = ProjectConfig.loadDefault()
    sections_2100hcv_1, sections_2100hcv_2 = make_random_sections(
        'HCV1A-H77-NS5a',
        1,
        300,
        projects,
        400)
    sections_2100v3_1, sections_2100v3_2 = (
        [FastqSection('HIV1-B-FR-K03455-seed', 7056, 7312, 50),
         FastqSection('HIV1-B-FR-K03455-seed', 7062, 7312, 50)],
        [FastqSection('HIV1-B-FR-K03455-seed', 7123, 7373, 50),
         FastqSection('HIV1-B-FR-K03455-seed', 7123, 7376, 50)])
    sections_2100hiv_1, sections_2100hiv_2 = make_random_sections(
        'RT',
        1,
        300,
        projects,
        400)
    sections_2160_1, sections_2160_2 = make_random_sections(
        'HCV2-JFH-1-NS5b',
        1,
        230,
        projects,
        mutations=(CodonMutation(159, 'GTC'),))
    sections_2160midi_1, sections_2160midi_2 = make_random_sections(
        'HCV2-JFH-1-NS5b',
        231,
        561,
        projects,
        mutations=(CodonMutation(316, 'AGC'),))
    sections_2170_1a_1, sections_2170_1a_2 = make_random_sections('HCV-1a',
                                                                  6258,
                                                                  9375)
    sections_2170_2_1, sections_2170_2_2 = make_random_sections('HCV-2a',
                                                                6269,
                                                                9440)
    sections_2180_1, sections_2180_2 = make_random_sections(
        'HIV1-B-FR-K03455-seed',
        6225,
        7757)
    hxb2_ref = projects.getReference('HIV1-B-FR-K03455-seed')

    sections_2220_mix1a_1, sections_2220_mix1a_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                    6225,
                                                                    7757,
                                                                    projects,
                                                                    9000
                                                                    )
    sections_2220_mix1b_1, sections_2220_mix1b_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                   6225,
                                                                   7757,
                                                                   projects,
                                                                   1000,
                                                                   mutations = (CodonMutation(7000, 'AAA'),)
                                                                   )
    sections_2220_mix2a_1, sections_2220_mix2a_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                    6225,
                                                                    7757,
                                                                    projects,
                                                                    8000
                                                                    )
    sections_2220_mix2b_1, sections_2220_mix2b_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                   6225,
                                                                   7757,
                                                                   projects,
                                                                   2000,
                                                                   mutations = (CodonMutation(7000, 'AAA'),)
                                                                   )
    sections_2220_mix3a_1, sections_2220_mix3a_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                    6225,
                                                                    7757,
                                                                    projects,
                                                                    7000
                                                                    )
    sections_2220_mix3b_1, sections_2220_mix3b_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                   6225,
                                                                   7757,
                                                                   projects,
                                                                   3000,
                                                                   mutations = (CodonMutation(7000, 'AAA'),)
                                                                   )
    sections_2220_mix4a_1, sections_2220_mix4a_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                    6225,
                                                                    7757,
                                                                    projects,
                                                                    6000
                                                                    )
    sections_2220_mix4b_1, sections_2220_mix4b_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                   6225,
                                                                   7757,
                                                                   projects,
                                                                   4000,
                                                                   mutations = (CodonMutation(7000, 'AAA'),)
                                                                   )
    sections_2220_mix5a_1, sections_2220_mix5a_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                    6225,
                                                                    7757,
                                                                    projects,
                                                                    5000
                                                                    )
    sections_2220_mix5b_1, sections_2220_mix5b_2 = make_random_sections('HIV1-B-FR-K03455-seed',
                                                                   6225,
                                                                   7757,
                                                                   projects,
                                                                   5000,
                                                                   mutations = (CodonMutation(7000, 'AAA'),)
                                                                   )

    projects.config['regions']['HXB2-with-deletion'] = dict(
        reference=hxb2_ref[617:928] + hxb2_ref[9358:9652],
        is_nucleotide=True,
        seed_group=None)
    sections_2210_1, sections_2210_2 = make_random_sections(
        'HXB2-with-deletion',
        projects=projects)
    fastq_files = [FastqFile('2010A-V3LOOP_S3_L001_R1_001.fastq',
                             '2010',
                             False,
                             (FastqSection('HIV1-CON-XX-Consensus-seed', 855, 906, 10),
                              FastqSection('HIV1-CON-XX-Consensus-seed', 912, 960, 10))),
                   FastqFile('2010A-V3LOOP_S3_L001_R2_001.fastq',
                             '2010',
                             True,
                             (FastqSection('HIV1-CON-XX-Consensus-seed', 855, 906, 10),
                              FastqSection('HIV1-CON-XX-Consensus-seed', 912, 960, 10))),
                   FastqFile('2020A-GP41_S4_L001_R1_001.fastq',
                             '2020',
                             False,
                             (FastqSection('HIV1-B-FR-KF716496-seed',
                                           6957,
                                           7065,
                                           10,
                                           (CodonMutation(6981, 'GGGATA'),)),)),
                   FastqFile('2020A-GP41_S4_L001_R2_001.fastq',
                             '2020',
                             True,
                             (FastqSection('HIV1-B-FR-KF716496-seed',
                                           6957,
                                           7065,
                                           10,
                                           (CodonMutation(6981, 'GGGATA'),)),)),
                   FastqFile('2040A-HLA-B_S6_L001_R1_001.fastq',
                             '2040',
                             False,
                             (FastqSection('HLA-B-seed', 201, 315, 80),
                              FastqSection('HLA-B-seed',
                                           201,
                                           315,
                                           20,
                                           (CodonMutation(207, 'TCT'),)))),
                   FastqFile('2040A-HLA-B_S6_L001_R2_001.fastq',
                             '2040',
                             True,
                             (FastqSection('HLA-B-seed', 201, 315, 80),
                              FastqSection('HLA-B-seed',
                                           201,
                                           315,
                                           20,
                                           (CodonMutation(207, 'TCT'),)))),
                   FastqFile('2070A-PR_S9_L001_R1_001.fastq',
                             '2070',
                             False,
                             (FastqSection('PR',
                                           40,
                                           80,
                                           12,
                                           (CodonMutation(45, ''),)),
                              FastqSection('PR',
                                           40,
                                           80,
                                           3,
                                           (CodonMutation(45, ''),
                                            CodonMutation(64, ''))))),
                   FastqFile('2070A-PR_S9_L001_R2_001.fastq',
                             '2070',
                             True,
                             (FastqSection('PR',
                                           40,
                                           80,
                                           12,
                                           (CodonMutation(45, ''),)),
                              FastqSection('PR',
                                           40,
                                           80,
                                           3,
                                           (CodonMutation(45, ''),
                                            CodonMutation(64, ''))))),
                   FastqFile('2100A-HCV-1337B-V3LOOP-PWND-HIV_S12_L001_R1_001.fastq',
                             '2100',
                             False,
                             sections_2100hcv_1 +
                             sections_2100v3_1 +
                             sections_2100hiv_1),
                   FastqFile('2100A-HCV-1337B-V3LOOP-PWND-HIV_S12_L001_R2_001.fastq',
                             '2100',
                             True,
                             sections_2100hcv_2 +
                             sections_2100v3_2 +
                             sections_2100hiv_2),
                   FastqFile('2130A-HCV_S15_L001_R1_001.fastq',
                             '2130',
                             False,
                             (FastqSection('HCV2-JFH-1-NS5b', 1, 66, 100),
                              FastqSection('HCV2-JFH-1-NS5b',
                                           115,
                                           181,
                                           100,
                                           (CodonMutation(159, 'GTC'),)))),
                   FastqFile('2130A-HCV_S15_L001_R2_001.fastq',
                             '2130',
                             True,
                             (FastqSection('HCV2-JFH-1-NS5b', 51, 114, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 165, 230, 100))),
                   FastqFile('2130AMIDI-MidHCV_S16_L001_R1_001.fastq',
                             '2130',
                             False,
                             (FastqSection('HCV2-JFH-1-NS5b', 231, 315, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 398, 485, 100))),
                   FastqFile('2130AMIDI-MidHCV_S16_L001_R2_001.fastq',
                             '2130',
                             True,
                             (FastqSection('HCV2-JFH-1-NS5b',
                                           305,
                                           397,
                                           100,
                                           (CodonMutation(316, 'AGC'),)),
                              FastqSection('HCV2-JFH-1-NS5b', 470, 561, 100))),
                   FastqFile('2140A-HIV_S17_L001_R1_001.fastq',
                             '2140',
                             False,
                             (FastqSection('PR',
                                           1,
                                           80,
                                           100,
                                           (CodonMutation(24, 'ATA'),)),)),
                   FastqFile('2140A-HIV_S17_L001_R2_001.fastq',
                             '2140',
                             True,
                             (FastqSection('PR',
                                           20,
                                           99,
                                           100,
                                           (CodonMutation(24, 'ATA'),)),)),
                   # Simplify with one_contig.
                   FastqFile('2160A-HCV_S19_L001_R1_001.fastq',
                             '2160',
                             False,
                             sections_2160_1),
                   FastqFile('2160A-HCV_S19_L001_R2_001.fastq',
                             '2160',
                             True,
                             sections_2160_2),
                   # Simplify with one_contig.
                   FastqFile('2160AMIDI-MidHCV_S20_L001_R1_001.fastq',
                             '2160',
                             False,
                             sections_2160midi_1),
                   FastqFile('2160AMIDI-MidHCV_S20_L001_R2_001.fastq',
                             '2160',
                             True,
                             sections_2160midi_2),
                   # Simplify with two_long_contigs.
                   FastqFile('2170A-HCV_S21_L001_R1_001.fastq',
                             '2170',
                             False,
                             sections_2170_1a_1 + sections_2170_2_1),
                   FastqFile('2170A-HCV_S21_L001_R2_001.fastq',
                             '2170',
                             True,
                             sections_2170_1a_2 + sections_2170_2_2),
                   FastqFile('2180A-HIV_S22_L001_R1_001.fastq',
                             '2180',
                             False,
                             sections_2180_1),
                   FastqFile('2180A-HIV_S22_L001_R2_001.fastq',
                             '2180',
                             True,
                             sections_2180_2),
                   FastqFile('2190A-SARSCOV2_S23_L001_R1_001.fastq',
                             '2190',
                             False,
                             (FastqSection('SARS-CoV-2-ORF1ab',
                                           4393,
                                           4429,
                                           50,
                                           (CodonMutation(4400, 'TCA'),)),
                              FastqSection('SARS-CoV-2-ORF1ab',
                                           4393,
                                           4430,
                                           50,
                                           (CodonMutation(4400, 'TCA'),)))),
                   FastqFile('2190A-SARSCOV2_S23_L001_R2_001.fastq',
                             '2190',
                             True,
                             (FastqSection('SARS-CoV-2-ORF1ab',
                                           4393,
                                           4429,
                                           50,
                                           (CodonMutation(4400, 'TCA'),)),
                              FastqSection('SARS-CoV-2-ORF1ab',
                                           4393,
                                           4430,
                                           50,
                                           (CodonMutation(4400, 'TCA'),)))),
                   FastqFile('2200A-SARSCOV2_S24_L001_R1_001.fastq',
                             '2200',
                             False,
                             (FastqSection('SARS-CoV-2-nsp1', 20, 66, 100),)),
                   FastqFile('2200A-SARSCOV2_S24_L001_R2_001.fastq',
                             '2200',
                             True,
                             (FastqSection('SARS-CoV-2-nsp1', 56, 102, 100),)),
                   FastqFile('2210A-NFLHIVDNA_S25_L001_R1_001.fastq',
                             '2210',
                             False,
                             sections_2210_1),
                   FastqFile('2210A-NFLHIVDNA_S25_L001_R2_001.fastq',
                             '2210',
                             True,
                             sections_2210_2),
                   FastqFile('2220-HIV-mixture10_S26_L001_R1_001.fastq',
                             '2220',
                             False,
                             sections_2220_mix1a_1 + sections_2220_mix1b_1),
                   FastqFile('2220-HIV-mixture10_S26_L001_R2_001.fastq',
                             '2220',
                             True,
                             sections_2220_mix1a_2 + sections_2220_mix1b_2),
                   FastqFile('2220-HIV-mixture20_S27_L001_R1_001.fastq',
                             '2220',
                             False,
                             sections_2220_mix2a_1 + sections_2220_mix2b_1),
                   FastqFile('2220-HIV-mixture20_S27_L001_R2_001.fastq',
                             '2220',
                             True,
                             sections_2220_mix2a_2 + sections_2220_mix2b_2),
                   FastqFile('2220-HIV-mixture30_S28_L001_R1_001.fastq',
                             '2220',
                             False,
                             sections_2220_mix3a_1 + sections_2220_mix3b_1),
                   FastqFile('2220-HIV-mixture30_S28_L001_R2_001.fastq',
                             '2220',
                             True,
                             sections_2220_mix3a_2 + sections_2220_mix3b_2),
                   FastqFile('2220-HIV-mixture40_S29_L001_R1_001.fastq',
                             '2220',
                             False,
                             sections_2220_mix4a_1 + sections_2220_mix4b_1),
                   FastqFile('2220-HIV-mixture40_S29_L001_R2_001.fastq',
                             '2220',
                             True,
                             sections_2220_mix4a_2 + sections_2220_mix4b_2),
                   FastqFile('2220-HIV-mixture50_S30_L001_R1_001.fastq',
                             '2220',
                             False,
                             sections_2220_mix5a_1 + sections_2220_mix5b_1),
                   FastqFile('2220-HIV-mixture50_S30_L001_R2_001.fastq',
                             '2220',
                             True,
                             sections_2220_mix5a_2 + sections_2220_mix5b_2)
                   ]
    for fastq_file in fastq_files:
        if not fastq_file.name.startswith('2220-HIV-mixture30'):
            continue
        with open(fastq_file.name, 'w') as f:
            next_cluster = 1
            for section in fastq_file.sections:
                ref_name, ref_start, ref_end = find_coord_pos(projects,
                                                              section.coord_name,
                                                              section.start_pos,
                                                              section.end_pos)

                ref_nuc_seq = projects.getReference(ref_name)
                ref_nuc_section = list(ref_nuc_seq[ref_start:ref_end])
                is_nucleotide = ((ref_start, ref_end) ==
                                 (section.start_pos, section.end_pos))
                for mutation in section.mutations:
                    if section.start_pos <= mutation.pos <= section.end_pos:
                        section_pos = mutation.pos - section.start_pos
                        if not is_nucleotide:
                            section_pos *= 3
                        ref_nuc_section[section_pos:section_pos+3] = list(mutation.codon)
                ref_nuc_section = ''.join(ref_nuc_section)
                if fastq_file.is_reversed:
                    ref_nuc_section = reverse_and_complement(ref_nuc_section)
                phred_scores = 'A' * len(ref_nuc_section)
                file_num = '2' if fastq_file.is_reversed else '1'
                # noinspection PyTypeChecker
                for cluster in range(section.count):
                    f.write('@M01234:01:000000000-AAAAA:1:1101:{}:{:04} {}:N:0:1\n'.format(
                        fastq_file.extract_num,
                        cluster + next_cluster,
                        file_num))
                    f.write(ref_nuc_section+'\n')
                    f.write('+\n')
                    f.write(phred_scores+'\n')
                next_cluster += section.count


def make_random_sections(coord_name,
                         min_start: int = None,
                         max_end: int = None,
                         projects=None,
                         read_count=10000,
                         mutations=()):
    if projects is None:
        ref_name = coord_name
        ref_start = min_start
        ref_end = max_end
    else:
        ref_name, ref_start, ref_end = find_coord_pos(projects,
                                                      coord_name,
                                                      min_start,
                                                      max_end)
    sections_1 = []
    sections_2 = []
    min_read_length = 20
    max_read_length = 250
    max_pair_length = 600
    for _ in range(read_count):
        start = randrange(ref_start, ref_end - min_read_length)
        end_stop = min(ref_end+1, start+max_pair_length)
        end = randrange(start + min_read_length, end_stop)
        end1 = min(end, start + max_read_length)
        start2 = max(start, end - max_read_length)
        sections_1.append(FastqSection(ref_name, start, end1, 1, mutations))
        sections_2.append(FastqSection(ref_name, start2, end, 1, mutations))
    return sections_1, sections_2


def find_coord_pos(projects: ProjectConfig,
                   coord_name: str,
                   start_pos: int = None,
                   end_pos: int = None):
    coord_seq = projects.getReference(coord_name)
    if start_pos is None:
        start_pos = 1
    if end_pos is None:
        end_pos = len(coord_seq) + 1
    if projects.config['regions'][coord_name]['is_nucleotide']:
        # Already have a nucleotide sequence, nothing to do.
        return coord_name, start_pos, end_pos
    gap_open = 40
    gap_extend = 10
    use_terminal_gap_penalty = 1
    highest_score = 0
    best_match = None
    ref_names = set()
    for project in projects.config['projects'].values():
        for region in project['regions']:
            if coord_name == region['coordinate_region']:
                ref_names.update(region['seed_region_names'])

    for ref_name in sorted(ref_names):
        ref_nuc_seq = projects.getReference(ref_name)
        for nuc_offset in range(3):
            ref_amino_seq = translate(ref_nuc_seq, nuc_offset)
            aligned_coord, aligned_ref, score = align_it_aa(
                coord_seq,
                ref_amino_seq,
                gap_open,
                gap_extend,
                use_terminal_gap_penalty)
            if score > highest_score:
                highest_score = score
                best_match = (ref_name, nuc_offset, aligned_coord, aligned_ref)
    ref_name, nuc_offset, aligned_coord, aligned_ref = best_match
    coord_pos = ref_pos = 0
    ref_start = ref_end = None
    for coord_amino, ref_amino in zip(aligned_coord, aligned_ref):
        if ref_amino != '-':
            ref_pos += 1
        if coord_amino != '-':
            coord_pos += 1
            if start_pos == coord_pos:
                ref_start = ref_pos * 3 - nuc_offset - 3
            if coord_pos == end_pos:
                ref_end = ref_pos * 3 - nuc_offset
    assert ref_start is not None
    assert ref_end is not None
    return ref_name, ref_start, ref_end


main()
