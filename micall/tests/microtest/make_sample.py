from collections import namedtuple
from gotoh import align_it_aa
from random import randrange

from micall.core.project_config import ProjectConfig
from micall.utils.translation import translate, reverse_and_complement

FastqSection = namedtuple('FastqSection', 'coord_name start_pos end_pos count')
CodonMutation = namedtuple('CodonMutation', 'pos codon')
FastqFile = namedtuple('FastqFile', 'name extract_num is_reversed sections mutations')


def main():
    projects = ProjectConfig.loadDefault()
    sections_2150_1, sections_2150_2 = make_random_sections('HCV-1a',
                                                            8000,
                                                            9000)
    sections_2160_1, sections_2160_2 = make_random_sections('HCV2-JFH-1-NS5b',
                                                            1,
                                                            230,
                                                            projects)
    sections_2160midi_1, sections_2160midi_2 = make_random_sections(
        'HCV2-JFH-1-NS5b',
        231,
        561,
        projects)
    sections_2170_1a_1, sections_2170_1a_2 = make_random_sections('HCV-1a',
                                                                  6258,
                                                                  9375)
    sections_2170_2_1, sections_2170_2_2 = make_random_sections('HCV-2a',
                                                                6269,
                                                                9440)
    fastq_files = [FastqFile('2130A-HCV_S15_L001_R1_001.fastq',
                             '2130',
                             False,
                             (FastqSection('HCV2-JFH-1-NS5b', 1, 60, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 117, 176, 100)),
                             (CodonMutation(159, 'GTC'),)),
                   FastqFile('2130A-HCV_S15_L001_R2_001.fastq',
                             '2130',
                             True,
                             (FastqSection('HCV2-JFH-1-NS5b', 57, 116, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 171, 230, 100)),
                             (CodonMutation(159, 'GTC'),)),
                   FastqFile('2130AMIDI-MidHCV_S16_L001_R1_001.fastq',
                             '2130',
                             False,
                             (FastqSection('HCV2-JFH-1-NS5b', 231, 313, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 396, 478, 100)),
                             (CodonMutation(316, 'AGC'),)),
                   FastqFile('2130AMIDI-MidHCV_S16_L001_R2_001.fastq',
                             '2130',
                             True,
                             (FastqSection('HCV2-JFH-1-NS5b', 313, 395, 100),
                              FastqSection('HCV2-JFH-1-NS5b', 479, 561, 100)),
                             (CodonMutation(316, 'AGC'),)),
                   FastqFile('2140A-HIV_S17_L001_R1_001.fastq',
                             '2140',
                             False,
                             (FastqSection('PR', 1, 80, 100),),
                             (CodonMutation(24, 'ATA'),)),
                   FastqFile('2140A-HIV_S17_L001_R2_001.fastq',
                             '2140',
                             True,
                             (FastqSection('PR', 20, 99, 100),),
                             (CodonMutation(24, 'ATA'),)),
                   FastqFile('2150A-HCV_S18_L001_R1_001.fastq',
                             '2150',
                             False,
                             sections_2150_1,
                             tuple()),
                   FastqFile('2150A-HCV_S18_L001_R2_001.fastq',
                             '2150',
                             True,
                             sections_2150_2,
                             tuple()),
                   FastqFile('2160A-HCV_S19_L001_R1_001.fastq',
                             '2160',
                             False,
                             sections_2160_1,
                             (CodonMutation(159, 'GTC'),)),
                   FastqFile('2160A-HCV_S19_L001_R2_001.fastq',
                             '2160',
                             True,
                             sections_2160_2,
                             (CodonMutation(159, 'GTC'),)),
                   FastqFile('2160AMIDI-MidHCV_S20_L001_R1_001.fastq',
                             '2160',
                             False,
                             sections_2160midi_1,
                             (CodonMutation(316, 'AGC'),)),
                   FastqFile('2160AMIDI-MidHCV_S20_L001_R2_001.fastq',
                             '2160',
                             True,
                             sections_2160midi_2,
                             (CodonMutation(316, 'AGC'),)),
                   FastqFile('2170A-HCV_S21_L001_R1_001.fastq',
                             '2170',
                             False,
                             sections_2170_1a_1 + sections_2170_2_1,
                             ()),
                   FastqFile('2170A-HCV_S21_L001_R2_001.fastq',
                             '2170',
                             True,
                             sections_2170_1a_2 + sections_2170_2_2,
                             ())]
    for fastq_file in fastq_files:
        with open(fastq_file.name, 'w') as f:
            next_cluster = 1
            for section in fastq_file.sections:
                ref_name, ref_start, ref_end = find_coord_pos(projects,
                                                              section.coord_name,
                                                              section.start_pos,
                                                              section.end_pos)

                ref_nuc_seq = projects.getReference(ref_name)
                ref_nuc_section = list(ref_nuc_seq[ref_start:ref_end])
                for mutation in fastq_file.mutations:
                    if section.start_pos <= mutation.pos <= section.end_pos:
                        section_pos = (mutation.pos - section.start_pos) * 3
                        ref_nuc_section[section_pos:section_pos+3] = list(mutation.codon)
                ref_nuc_section = ''.join(ref_nuc_section)
                if fastq_file.is_reversed:
                    ref_nuc_section = reverse_and_complement(ref_nuc_section)
                phred_scores = 'A' * len(ref_nuc_section)
                file_num = '2' if fastq_file.is_reversed else '1'
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
                         min_start,
                         max_end,
                         projects=None,
                         read_count=10000):
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
        end = randrange(start + min_read_length, start + max_pair_length)
        end1 = min(end, start + max_read_length)
        start2 = max(start, end - max_read_length)
        sections_1.append(FastqSection(ref_name, start, end1, 1))
        sections_2.append(FastqSection(ref_name, start2, end, 1))
    return sections_1, sections_2


def find_coord_pos(projects, coord_name, start_pos, end_pos):
    if projects.config['regions'][coord_name]['is_nucleotide']:
        # Already have a nucleotide sequence, nothing to do.
        return coord_name, start_pos, end_pos
    coord_seq = projects.getReference(coord_name)
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
