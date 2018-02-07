from gotoh import align_it_aa

from micall.core.project_config import ProjectConfig
from micall.utils.translation import translate, reverse_and_complement


def main():
    coord_name = 'HCV2-JFH-1-NS5b'
    extract_num = '2130'
    start_pos, end_pos = 396, 561
    read_count = 100
    ref_name = None  # 'HCV-2a'
    ref_start, ref_end = 0, 0
    is_reversed = True
    projects = ProjectConfig.loadDefault()
    if ref_name is None:
        ref_name, ref_start, ref_end = find_coord_pos(projects,
                                                      coord_name,
                                                      start_pos,
                                                      end_pos)
    ref_nuc_seq = projects.getReference(ref_name)
    ref_nuc_section = list(ref_nuc_seq[ref_start:ref_end])
    # ref_nuc_section[(316-231)*3] = 'A'
    ref_nuc_section = ''.join(ref_nuc_section)
    # print(ref_nuc_section[(316-231)*3:(317-231)*3])
    if is_reversed:
        ref_nuc_section = reverse_and_complement(ref_nuc_section)
    phred_scores = 'A' * (ref_end-ref_start)
    file_num = '2' if is_reversed else '1'
    for cluster in range(read_count):
        print('@M01234:01:000000000-AAAAA:1:1101:{}:{:04} {}:N:0:1'.format(
            extract_num,
            cluster+1,
            file_num))
        print(ref_nuc_section)
        print('+')
        print(phred_scores)


def find_coord_pos(projects, coord_name, start_pos, end_pos):
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
    print(highest_score, ref_name, nuc_offset)
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
    print(ref_start, ref_end)
    return ref_name, ref_start, ref_end


main()
