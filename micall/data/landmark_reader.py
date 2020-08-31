import re
import typing
from operator import itemgetter
from pathlib import Path

from yaml import safe_load

DEFAULT_PATH = Path(__file__).parent / 'landmark_references.yaml'


class LandmarkReader:
    @classmethod
    def load(cls, f: typing.TextIO = None):
        """ Load an instance of this class from an open JSON file.

        :param f: The file to load from, or None to load from the default.
        """
        if f is None:
            with DEFAULT_PATH.open() as f:
                return LandmarkReader.load(f)
        return LandmarkReader(safe_load(f))

    def __init__(self, landmarks: list):
        self.landmarks = landmarks

    def get_gene(self, coordinates: str, gene_name: str):
        """ Get the details for a gene in the given set of coordinates.

        :param coordinates: the group of seeds that are displayed together
        :param gene_name: the gene within the group
        """
        genotype_landmarks, = (entry
                               for entry in self.landmarks
                               if entry['coordinates'] == coordinates)
        regions = sorted(genotype_landmarks['landmarks'], key=itemgetter('start'))
        prefix = genotype_landmarks.get('prefix', '')
        if not gene_name.startswith(prefix):
            raise ValueError(f'Gene name {gene_name!r} does not start with '
                             f'prefix {prefix!r}.')
        gene_name = gene_name[len(prefix):]
        for i, region in enumerate(regions):
            full_name = region.get('full_name')
            if full_name is None:
                full_name = region.get('name')
            if full_name != gene_name:
                continue
            if 'end' not in region:
                region['end'] = regions[i+1]['start']
            break
        else:
            raise ValueError(f'Landmarks not found for gene {gene_name!r} in '
                             f'{coordinates}.')
        return region

    def get_coordinates(self, seed_name: str) -> str:
        for genotype_landmarks in self.landmarks:
            seed_pattern = genotype_landmarks['seed_pattern']
            if re.fullmatch(seed_pattern, seed_name):
                return genotype_landmarks['coordinates']
        raise ValueError(f'No landmarks match {seed_name!r}.')
