from mappy import Aligner
from enum import IntEnum
from pathlib import Path
from pprint import pprint
from collections import Counter
import yaml
import argparse
from micall.core.project_config import ProjectConfig
from micall.data.landmark_reader import LandmarkReader

CigarActions = IntEnum(
    'CigarActions',
    'MATCH INSERT DELETE SKIPPED SOFT_CLIPPED HARD_CLIPPED',
    start=0)
CigarActionsToLetters = {0: 'M', 1: 'I', 2: 'D'}


class AlignmentEvaluator:
    def __init__(self,
                 ref_name='HIV1-B-FR-K03455-seed',
                 project='HIV',
                 aligner_preset='asm10',
                 verbose=0,
                 seed_group=None):
        self.max_deletion = 21
        self.max_insertion = 21
        self.min_coverage = 0.8
        self.ref_name = ref_name
        self.project_code = project
        self.projects = ProjectConfig.loadDefault()
        landmarks_path = (Path(__file__).parent.parent / 'data' / 'landmark_references.yaml')
        landmarks_yaml = landmarks_path.read_text()
        landmarks = yaml.safe_load(landmarks_yaml)
        self.landmark_reader = LandmarkReader(landmarks)
        self.main_reference = self.projects.getReference(ref_name)
        self.seeds = sorted(self.projects.getProjectSeeds(project))
        self.seed_group = seed_group
        if seed_group is not None:
            group_seeds = []
            for seed in self.seeds:
                split_name = seed.split('-')
                current_seed_group = '-'.join(split_name[0:2])
                if current_seed_group == seed_group:
                    group_seeds.append(seed)
            self.seeds = group_seeds
        self.aligner_preset = aligner_preset
        self.aligner = Aligner(seq=self.main_reference, preset=aligner_preset)
        self.num_warning_seeds = 0
        self.verbose = verbose
        self.warning_seeds = {'insertions': [], 'deletions': [], 'frameshift': [], 'coverage': [], 'alignment': []}

    def evaluate_alignment(self):
        for seed in self.seeds:
            seed_alignment = SeedAlignment(seed,
                                           self.ref_name,
                                           self.aligner_preset,
                                           self.project_code,
                                           self.verbose,
                                           self.seed_group)
            warned = seed_alignment.evaluate_seed_alignment()
            self.update_summary(seed_alignment)
            self.num_warning_seeds += warned
        if self.verbose == 1:
            self.print_summary()
        elif self.verbose == 0:
            self.seed_group_summary()
        print(f"There were warnings for {self.num_warning_seeds} of {len(self.seeds)} seeds")

    def update_summary(self, seed_alignment):
        for (warning, list) in seed_alignment.warning_dict.items():
            if len(list) > 0:
                self.warning_seeds[warning].append(seed_alignment.seed_name)

    def print_summary(self):
        warning_start = 'Total number of references with '
        print(warning_start + f' warnings: {self.num_warning_seeds} of {len(self.seeds)}')
        for (warning, seeds) in self.warning_seeds.items():
            if warning == 'insertions':
                print(warning_start + f'insertions larger than {self.max_insertion}: {len(seeds)}')
            elif warning == 'deletions':
                print(warning_start + f'deletions larger than {self.max_deletion}: {len(seeds)}')
            elif warning == 'coverage':
                print(warning_start + f'coverage smaller than {int(100*self.min_coverage)}%: {len(seeds)}')
            elif warning == 'alignment':
                print(warning_start + f'more than one alignment: {len(seeds)}')
            elif warning == 'frameshift':
                print(warning_start + f'frameshifts in translated regions: {len(seeds)}')
            if len(seeds) > 0:
                print(f'    Those were: ')
                pprint(seeds)

    def seed_group_summary(self):
        seed_groups_counts = Counter()
        seed_group_warnings = {key: Counter() for key in self.warning_seeds.keys()}
        for warning in self.warning_seeds.keys():
            for seed in self.warning_seeds[warning]:
                split_name = seed.split('-')
                seed_group = '-'.join(split_name[0:2])
                seed_group_warnings[warning].update([seed_group])
        for seed in self.seeds:
            split_name = seed.split('-')
            seed_group = '-'.join(split_name[0:2])
            seed_groups_counts.update([seed_group])
        for seed_group in seed_groups_counts.keys():
            print(f'Seed group {seed_group}')
            for warning in seed_group_warnings.keys():
                if warning == 'insertions':
                    message = f'had a large insertion (over {self.max_insertion})'
                elif warning == 'deletions':
                    message = f'had a large deletion (over {self.max_deletion})'
                elif warning == 'alignment':
                    message = 'had more than one alignment'
                elif warning == 'coverage':
                    message = f'had insufficient coverage (less than {int(100*self.min_coverage)}%)'
                else:  # warning = 'frameshift':
                    message = 'had a frame shift in a translated region'
                num_warning = seed_group_warnings[warning].get(seed_group)
                num_seeds = seed_groups_counts[seed_group]
                if num_warning is not None:
                    print(f'{num_warning} of {num_seeds} ' + message)


class SeedAlignment(AlignmentEvaluator):
    def __init__(self,
                 seed_name,
                 ref_name='HIV1-B-FR-K03455-seed',
                 project='HIV',
                 aligner_preset='asm10',
                 verbose=0,
                 seed_group = None):
        AlignmentEvaluator.__init__(self, ref_name, aligner_preset, project, verbose, seed_group)
        self.seed_name = seed_name
        self.seed_ref = self.projects.getReference(seed_name)
        self.alignments = list(self.aligner.map(self.seed_ref))
        self.seed_coverage = 0
        self.warning_dict = {key: [] for key in self.warning_seeds.keys()}
        self.warnings = ''

    def evaluate_seed_alignment(self):
        if len(self.alignments) > 1:
            self.warning_dict['alignment'].append(self.seed_name)
            self.warnings += 'Warning: more than one alignment for this seed\n'
            self.alignments.sort(key=lambda item: item.r_st)
        for alignment in self.alignments:
            self.seed_coverage += self.parse_alignment(alignment, 0, len(self.main_reference), 'overall', False)
            for entry in self.projects.config['projects'][self.project_code]['regions']:
                region_name = entry['coordinate_region']
                region_info = self.landmark_reader.get_gene(
                    self.ref_name,
                    region_name,
                    drop_stop_codon=False)
                region_start = region_info['start']
                region_end = region_info['end']
                is_amino = self.projects.isAmino(region_name)
                self.parse_alignment(alignment, region_start, region_end, region_name, is_amino)
        self.seed_coverage = self.seed_coverage / len(self.main_reference)
        if self.seed_coverage < self.min_coverage:
            self.warning_dict['coverage'].append(self.seed_name)
            self.warnings += f'Insufficient genome coverage of aligned seed: {int(self.seed_coverage * 100)}%\n'
        if self.verbose >= 2:
            self.print_seed_details()
        if self.warnings != '':
            return 1
        else:
            return 0

    def print_seed_details(self):
        if self.warnings != '':
            print(f'Warnings for seed {self.seed_name}:')
            print(f'Length of seed: {len(self.seed_ref)}')
            for i, alignment in enumerate(self.alignments):
                print(f'Alignment {i + 1} from ref {alignment.r_st} to {alignment.r_en}, '
                      f'query {alignment.q_st} to {alignment.q_en}')
                print(alignment.cigar_str)
            print(self.warnings)
        else:
            print(f'No warnings for seed {self.seed_name}.')
            print(f'Length of seed: {len(self.seed_ref)}\n')

    def parse_alignment(self, alignment, region_start, region_end, region_name, warn_indels):
        ref_index = alignment.r_st
        removed_skipped_pos = False
        region_cigar = ''
        indels = 0
        coverage = 0
        for cigar_index, (size, action) in enumerate(alignment.cigar):
            cigar_start = ref_index
            cigar_end = ref_index + size if (action == CigarActions.DELETE or action == CigarActions.MATCH)\
                else ref_index
            if cigar_start < region_start:
                # not started in the region yet, or very first cigar (can be ignored). Only advance ref_index
                ref_index = cigar_end
                if cigar_end > region_start:
                    # first cigar, write partial to region cigar
                    if action != CigarActions.INSERT:
                        region_cigar += str(cigar_end-region_start) + CigarActionsToLetters.get(action) + '.'
                    if action == CigarActions.MATCH:
                        coverage += cigar_end-region_start
                continue
            if region_start <= cigar_start:
                if cigar_end >= region_end:
                    # last cigar, write partial to region cigar
                    if action != CigarActions.INSERT:
                        region_cigar += str(region_end-cigar_start) + CigarActionsToLetters.get(action) + '.'
                    if action == CigarActions.MATCH:
                        coverage += region_end-cigar_start
                    break
            region_cigar += str(size) + CigarActionsToLetters.get(action) + '.'
            if action == CigarActions.DELETE:
                ref_index += size
                if 5770 <= ref_index <= 5780 and removed_skipped_pos is False and self.project_code == 'HIV':
                    # don't count deleted skipped position
                    removed_skipped_pos = True
                    size -= 1
                if size > self.max_deletion and region_name == 'overall':
                    regions_del_start = self.landmark_reader.get_region(self.ref_name, ref_index)
                    regions_del_end = self.landmark_reader.get_region(self.ref_name, ref_index-size)
                    regions = list(set(regions_del_start) | set(regions_del_end))
                    self.warning_dict['deletions'] += regions
                    self.warnings += f'Warning: large deletion of size {size} at ref position {ref_index},' \
                                     f'in the following regions: {regions}\n'
                indels -= size
            elif action == CigarActions.INSERT:
                if size > self.max_insertion and region_name == 'overall':
                    regions_ins = self.landmark_reader.get_region(self.ref_name, ref_index)
                    self.warning_dict['insertions'] += regions_ins
                    self.warnings += f'Warning: large insertion of size {size} at ref position {ref_index},' \
                                     f'in the following regions: {regions_ins}\n'
                indels += size
            elif action == CigarActions.MATCH:
                ref_index += size
                coverage += size
            else:
                self.warnings += 'Weird alignment stuff?\n'
        if indels % 3 != 0 and warn_indels:
            self.warning_dict['frameshift'].append(region_name)
            self.warnings += f'Warning: deletions and insertions not a multiple of 3 in region {region_name}\n'
            self.warnings += f'The sum of all insertions and deletions is {indels}.\n'
            self.warnings += f'Cigar string of {region_name}: {region_cigar}\n'
        if region_name in self.warning_dict['insertions']:
            self.warnings += f'Region {region_name} has a large insertion (see above).\n'
            self.warnings += f'Cigar string of {region_name}: {region_cigar}\n'
        if region_name in self.warning_dict['deletions']:
            self.warnings += f'Region {region_name} has a large deletion (see above).\n'
            self.warnings += f'Cigar string of {region_name}: {region_cigar}\n'
        return coverage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aligner_preset',
                        type=str,
                        default='asm10',
                        help="Minimap aligner preset option")
    parser.add_argument('--verbose', '-v',
                        action='count',
                        default=0,
                        help="Verbosity: -vv gives individual seed details,"
                             " -v gives individual seed summary, "
                             "and nothing will result in seed group summary.")
    parser.add_argument('--seed_group',
                        type=str,
                        default=None,
                        help="Use this option to investigate the seeds of a specific seed group.")
    parser.add_argument('--ref_name',
                        type=str,
                        default='HIV1-B-FR-K03455-seed',
                        help="Reference to compare seeds to")
    parser.add_argument('--project_code',
                        type=str,
                        default='HIV',
                        help='Project code for which to investigate seeds')
    args = parser.parse_args()

    alignment_evaluator = AlignmentEvaluator(aligner_preset=args.aligner_preset,
                                             verbose=args.verbose,
                                             seed_group=args.seed_group,
                                             ref_name=args.ref_name,
                                             project=args.project_code)
    # 3 verbosity levels: 0 - seed group summary, 1 - individual seed summary, 2 - individual seed details
    alignment_evaluator.evaluate_alignment()


if __name__ == '__main__':
    main()
