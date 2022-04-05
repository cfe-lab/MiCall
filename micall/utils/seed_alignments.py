from micall.utils.consensus_aligner import determine_region_positions
from micall.core.project_config import ProjectConfig
from micall.data.landmark_reader import LandmarkReader
from pathlib import Path
import yaml

landmarks_path = (Path(__file__).parent.parent / 'data' / 'landmark_references.yaml')
landmarks_yaml = landmarks_path.read_text()
landmark_reader = LandmarkReader(yaml.safe_load(landmarks_yaml))
projects = ProjectConfig.loadDefault()
project_codes = ['HCV', 'HIVB', 'SARSCOV2']
for code in project_codes:
    seeds = projects.getProjectSeeds(code)
    for seed in seeds:
        print(f"--------- Seed: {seed} -------------")
        seed_ref = projects.getReference(seed)
        for region in projects.getCoordinateReferences(seed):
            if not projects.isAmino(region):
                continue
            region_ref = projects.getReference(region)
            coordinate_name = landmark_reader.get_coordinates(seed)
            gene_info = landmark_reader.get_gene(coordinate_name, region)
            start_pos, end_pos = determine_region_positions(seed_ref, region_ref)
            ref_start = gene_info['start']
            ref_end = gene_info['end']
            ref_len = ref_end - ref_start
            len = end_pos - start_pos
            if abs(start_pos - ref_start) > 1000 or abs(end_pos - ref_end) > 1000 or start_pos == -1 or end_pos == -1 \
                    or not (1.2*ref_len > len > 0.8*ref_len):
                print(f"region: {region}, coordinates {coordinate_name}")
                print(f"start pos: {start_pos}, ref start {ref_start}")
                print(f"end pos: {end_pos}, ref end {ref_end}")