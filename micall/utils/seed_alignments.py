from micall.utils.consensus_aligner import determine_region_positions
from micall.core.project_config import ProjectConfig
from micall.data.landmark_reader import LandmarkReader
from pathlib import Path
import yaml

landmarks_path = (Path(__file__).parent.parent / 'data' / 'landmark_references.yaml')
landmarks_yaml = landmarks_path.read_text()
landmark_reader = LandmarkReader(yaml.safe_load(landmarks_yaml))
regions_dict = {}
projects = ProjectConfig.loadDefault()
project_codes = ['HCV', 'HIVB', 'SARSCOV2']
for code in project_codes:
    regions_dict[code] = {}
    seeds = projects.getProjectSeeds(code)
    print(f"----------- Project code: {code} ---------------")
    num_problems = 0
    num_seeds = 0
    for seed in seeds:
        regions_dict[code][seed] = {}
        print(f"--------- Seed: {seed} -------------")
        seed_ref = projects.getReference(seed)
        num_seeds += 1
        # only do V3LOOP for consensus seed:
        if seed == 'HIV1-CON-XX-Consensus-seed':
            regions = ["V3LOOP"]
        else:
            regions = projects.getCoordinateReferences(seed)
            try:
                regions.pop("V3LOOP")
            except KeyError:
                pass
        for region in regions:
            if not projects.isAmino(region):
                continue
            regions_dict[code][seed][region] = {}
            region_ref = projects.getReference(region)
            coordinate_name = landmark_reader.get_coordinates(seed)
            gene_info = landmark_reader.get_gene(coordinate_name, region)
            start_pos, end_pos = determine_region_positions(seed_ref, region_ref)
            ref_start = gene_info['start']
            ref_end = gene_info['end']
            ref_length = ref_end - ref_start
            length = end_pos - start_pos
            seed_length = len(seed_ref)
            if abs(start_pos - ref_start) > 1000 or abs(end_pos - ref_end) > 1000 or start_pos == -1 or end_pos == -1 \
                    or not (1.2*ref_length > length > 0.8*ref_length):
                if region == "HIV1B-nef" and seed_length-3 <= end_pos <= len(seed_ref):
                    # nef is cut off at the end - that's fine
                    continue
                if region == "V3LOOP" and (1.2*ref_length > length > 0.8*ref_length):
                    # consensus seed has different coordinates - as long as the region is long enough, it's fine
                    continue
                print(f"region: {region}, coordinates {coordinate_name}")
                print(f"start pos: {start_pos}, ref start {ref_start}")
                print(f"end pos: {end_pos}, ref end {ref_end}")
                num_problems += 1
            else:
                regions_dict[code][seed][region]['start'] = start_pos
                regions_dict[code][seed][region]['end'] = end_pos
    print(f"{num_seeds} seeds, {num_problems} seeds with issues")

output_path = (Path(__file__).parent.parent / 'data' / 'seed_coordinates.yaml')
with open(output_path, 'w') as file:
    yaml.dump(regions_dict, file)