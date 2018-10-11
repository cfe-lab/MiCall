from collections import namedtuple
from gotoh import align_it_aa
from random import randrange

from micall.core.project_config import ProjectConfig
from micall.utils.translation import translate, reverse_and_complement

FastqSection = namedtuple('FastqSection', 'coord_name start_pos end_pos count')
CodonMutation = namedtuple('CodonMutation', 'pos codon')
FastqFile = namedtuple('FastqFile', 'name extract_num is_reversed sections mutations')


def main():
    sections_2140_1, sections_2140_2 = make_sections_2140()
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
                   FastqFile('2140A-HCV_S17_L001_R1_001.fastq',
                             '2140',
                             False,
                             sections_2140_1,
                             tuple()),
                   FastqFile('2140A-HCV_S17_L001_R2_001.fastq',
                             '2140',
                             True,
                             sections_2140_2,
                             tuple())]
    projects = ProjectConfig.loadDefault()
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
                phred_scores = 'A' * (ref_end-ref_start)
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


def make_sections_2140():
    sections_2140_1 = []
    sections_2140_2 = []
    ref_name = 'HCV-1a'
    read_count = 10000
    min_read_length = 20
    max_read_length = 250
    max_pair_length = 600
    min_start = 8000
    max_end = 9000
    for _ in range(read_count):
        start = randrange(min_start, max_end - min_read_length)
        end = randrange(start + min_read_length, start + max_pair_length)
        end1 = min(end, start + max_read_length)
        start2 = max(start, end - max_read_length)
        sections_2140_1.append(FastqSection(ref_name, start, end1, 1))
        sections_2140_2.append(FastqSection(ref_name, start2, end, 1))
    return sections_2140_1, sections_2140_2


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
    for ref_name in sorted(projects.getProjectSeeds('HCV')):
        if not ref_name.startswith('HCV-2'):
            continue
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
        if coord_amino != '-':
            coord_pos += 1
        if ref_amino != '-':
            ref_pos += 1
        if start_pos == coord_pos:
            ref_start = ref_pos * 3 - nuc_offset - 3
        if coord_pos == end_pos:
            ref_end = ref_pos * 3 - nuc_offset
    return ref_name, ref_start, ref_end


main()
