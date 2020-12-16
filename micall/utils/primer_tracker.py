import typing
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from micall.core.project_config import ProjectConfig
from micall.data.landmark_reader import LandmarkReader
from micall.utils.alignment_wrapper import align_nucs


class PrimerTracker:
    @staticmethod
    def load_locations() -> typing.Dict[str, typing.Dict[str, typing.Set[int]]]:
        """ Load ignored primer locations for each HCV gene.

        :return: {seed_name: {gene_name: {ignored_pos}}}
        """
        locations = defaultdict(lambda: defaultdict(set))

        projects = ProjectConfig.loadDefault()
        all_refs = projects.getAllReferences()
        hcv_refs = {name: ref
                    for name, ref in all_refs.items()
                    if name.startswith('HCV-')}
        primers = {}
        core_path = Path(__file__).parent / 'micall' / 'core'
        left_primers_path = core_path / 'primers_sarscov2_left.fasta'
        right_primers_path = core_path / 'primers_sarscov2_right_end.fasta'
        for primers_path in (left_primers_path, right_primers_path):
            with primers_path.open() as f:
                for primer in SeqIO.parse(f, 'fasta'):
                    primer_name = primer.name
                    if not primer_name.startswith('HCV'):
                        continue
                    if 'dA20' in primer_name or 'TIM' in primer_name:
                        continue
                    primer_seq = str(primer.seq).replace('X', '')
                    primers[primer_name] = primer_seq
        landmark_reader = LandmarkReader.load()
        for seed_name, seed_seq in hcv_refs.items():
            region_names = list(projects.getCoordinateReferences(seed_name))
            coordinate_name = landmark_reader.get_coordinates(seed_name)
            if coordinate_name != seed_name:
                locations[seed_name] = locations[coordinate_name]
                continue
            for primer_name, primer_seq in primers.items():
                aseed, aprimer, score = align_nucs(seed_seq, primer_seq)
                primer_end = seed_pos = 0
                primer_start = len(aseed)
                for seed_nuc, primer_nuc in zip(aseed, aprimer):
                    if seed_nuc != '-':
                        seed_pos += 1
                    if primer_nuc != '-':
                        primer_start = min(primer_start, seed_pos)
                        primer_end = max(primer_end, seed_pos)
                for region in region_names:
                    region_details = landmark_reader.get_gene(seed_name, region)
                    region_start = region_details['start']
                    region_end = region_details['end']
                    if region_start <= primer_end and primer_start <= region_end:
                        start = max(1, primer_start - region_start + 1)
                        end = min(region_end-region_start+1,
                                  primer_end-region_start+1)
                        # print(f'Adding {seed_name}, {region}: {start}->{end} '
                        #       f'({region_start=}, {region_end=},'
                        #       f'{primer_start=}, {primer_end=})')
                        locations[seed_name][region] |= set(range(start, end+1))
        return locations

    locations = None

    def __init__(self):
        pass

    def is_ignored(self, seed: str, region: str, pos: int):
        if PrimerTracker.locations is None:
            PrimerTracker.locations = self.load_locations()
        ignored_positions = self.locations[seed][region]
        return pos in ignored_positions
