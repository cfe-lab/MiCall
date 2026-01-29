from collections import defaultdict
import json
import logging
import re

from gotoh import align_it  # @UnresolvedImport
import Levenshtein
import matplotlib.pyplot as plt

from micall.utils.externals import ProjectsFile

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('reference_distances')
key_references = []
key_sections = {}  # {reference: [key]}


def populate_key_references(all_regions):
    DEFINITIONS = [('HCV1A-H77-NS3-seed', 701, 950),
                   ('HCV1A-H77-NS5a-seed', 1, 250),
                   ('HCV1A-H77-NS5b-seed', 101, 350)]
    for name, start, end in DEFINITIONS:
        reference = ''.join(all_regions[name]['reference'])
        key_references.append(reference[start-1:end])


def calculate_keys(reference):
    keys = key_sections.get(reference, None)
    if keys is not None:
        return keys

    GAP_INIT_PENALTY = 10
    GAP_EXTEND_PENALTY = 10
    USE_TERMINAL_GAP_PENALTY = False
    keys = []
    for key in key_references:
        # s1 is large sequence, s2 is key region
        aligned_source, aligned_key, _score = align_it(reference,
                                                       key,
                                                       GAP_INIT_PENALTY,
                                                       GAP_EXTEND_PENALTY,
                                                       USE_TERMINAL_GAP_PENALTY)
        match = re.match('^-*(.*?)-*$', aligned_key)
        excerpt = aligned_source[match.start(1):match.end(1)].replace('-', '')
        keys.append(excerpt)
    key_sections[reference] = keys
    return keys


def calculate_distance(source, destination):
    source_keys = calculate_keys(source)
    dest_keys = calculate_keys(destination)
    distance = 0
    for source_key, dest_key in zip(source_keys, dest_keys):
        distance += Levenshtein.distance(source_key, dest_key)
    return distance


def plot_distances(projects_filename):
    with open(projects_filename, 'r') as f:
        config = json.load(f)
    populate_key_references(config['regions'])
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
    median_references = []
    group_labels = []
    for group_index, group_name in enumerate(group_names):
        logger.info('Grouping %s.', group_name)
        seed_names = groups[group_name]
        seed_names.sort()
        source_seed_names.append(seed_names[0])
        references = []
        for seed_name in seed_names:
            reference = ''.join(config['regions'][seed_name]['reference'])
            all_seeds[seed_name] = (group_index, reference)
            references.append(reference)
        median_references.append(Levenshtein.median(references))
        group_labels.append(group_name[4:-6])  # trim HCV- and -seeds
    config = None

    intragroup_source_groups = []
    intragroup_distances = []
    intergroup_source_groups = []
    intergroup_distances = []

    for source_index, source_group_name in enumerate(group_names):
        logger.info('Processing %s.', source_group_name)
        source_reference = median_references[source_index]
        for dest_index, dest_reference in all_seeds.itervalues():
            distance = calculate_distance(source_reference, dest_reference)
            if source_index == dest_index:
                intragroup_source_groups.append(source_index)
                intragroup_distances.append(distance)
            else:
                intergroup_source_groups.append(source_index)
                intergroup_distances.append(distance)

    fig = plt.figure()
    ax = fig.add_subplot(111,
                         title='Distance From Genotype Median Reference in Key Regions',
                         xlabel='genotype',
                         ylabel='Levenshtein distance',
                         xticks=range(len(group_labels)),
                         xticklabels=group_labels)
    ax.plot(intragroup_source_groups, intragroup_distances, 'go', alpha=0.4)
    ax.plot(intergroup_source_groups, intergroup_distances, 'ro', alpha=0.4)
    ax.margins(0.1)
    plt.show()

with ProjectsFile().path() as projects_file_path:
    plot_distances(projects_file_path)
