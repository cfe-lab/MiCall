from collections import defaultdict
import json

import matplotlib.pyplot as plt
import Levenshtein
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('reference_distances')


def plot_distances(projects_filename):
    with open(projects_filename, 'rU') as f:
        config = json.load(f)
    groups = defaultdict(list)
    for name, region in config['regions'].iteritems():
        seed_group = region['seed_group']
        if seed_group and seed_group.startswith('HCV-'):
            groups[seed_group].append(name)
    del groups['HCV-seeds']
    group_names = groups.keys()
    group_names.sort()
    source_seed_names = []
    all_seeds = {}  # {name: (group_index, reference)}
    for group_index, group_name in enumerate(group_names):
        seed_names = groups[group_name]
        seed_names.sort()
        source_seed_names.append(seed_names[0])
        for seed_name in seed_names:
            reference = ''.join(config['regions'][seed_name]['reference'])
            all_seeds[seed_name] = (group_index, reference)
    config = None

    intragroup_source_groups = []
    intragroup_distances = []
    intergroup_source_groups = []
    intergroup_distances = []

    for source_seed_name in source_seed_names:
        logger.info('Processing %s.', source_seed_name)
        source_index, source_reference = all_seeds[source_seed_name]
        for dest_seed_name, (dest_index, dest_reference) in all_seeds.iteritems():
            if dest_seed_name != source_seed_name:
                distance = Levenshtein.distance(source_reference, dest_reference)
                if source_index == dest_index:
                    intragroup_source_groups.append(source_index)
                    intragroup_distances.append(distance)
                else:
                    intergroup_source_groups.append(source_index)
                    intergroup_distances.append(distance)

    plt.plot(intragroup_source_groups, intragroup_distances, 'go')
    plt.plot(intergroup_source_groups, intergroup_distances, 'ro')
    plt.xticks(range(len(source_seed_names)), source_seed_names, rotation='vertical')
    plt.margins(0.1)
    plt.subplots_adjust(bottom=0.15)
    plt.show()

plot_distances('../projects.json')
