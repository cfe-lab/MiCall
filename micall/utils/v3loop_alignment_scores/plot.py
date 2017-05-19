from collections import Counter
from csv import DictReader, DictWriter
from operator import itemgetter
import os

from matplotlib import pyplot as plt

from ...core.project_config import ProjectConfig
from ...g2p.fastq_g2p import HIV_SEED_NAME, COORDINATE_REF_NAME, extract_target, \
    FastqReader, merge_reads, trim_reads

CURRENT_MINIMUM_SCORE = 52
NEW_MINIMUM_SCORE = 200


def align_reads(fastq1, fastq2):
    project_config = ProjectConfig.loadDefault()
    hiv_seed = project_config.getReference(HIV_SEED_NAME)
    coordinate_ref = project_config.getReference(COORDINATE_REF_NAME)
    v3loop_ref = extract_target(hiv_seed, coordinate_ref)
    reader = FastqReader(fastq1, fastq2)
    merged_reads = merge_reads(reader)
    score_counts = Counter()
    for _ in trim_reads(merged_reads, v3loop_ref, score_counts):
        pass
    return score_counts


def plot_file(filename1):
    base_name = os.path.basename(filename1)
    name_parts = base_name.split('_')
    work_path = os.path.dirname(__file__)
    scores_filename = os.path.join(
        work_path,
        '_'.join(name_parts[:2] + ['v3loop_scores.csv']))
    if os.path.exists(scores_filename):
        with open(scores_filename) as f:
            reader = DictReader(f)
            score_rows = [list(map(int, row))
                          for row in map(itemgetter('score', 'count'), reader)]
    else:
        source1 = os.path.join('micall/tests/working/v3loop_alignment_scores/',
                               filename1)
        source2 = source1.replace('_R1_', '_R2_')
        with open(source1) as fastq1, open(source2) as fastq2:
            score_counts = align_reads(fastq1, fastq2)
        score_rows = sorted(score_counts.items())
        with open(scores_filename, 'w') as scores_csv:
            writer = DictWriter(scores_csv,
                                ('score', 'count'),
                                lineterminator=os.linesep)
            writer.writeheader()
            for score, count in score_rows:
                writer.writerow(dict(score=score, count=count))
    scores = [row[0] for row in score_rows]
    counts = [row[1] for row in score_rows]
    total_count = sum(counts)
    fractions = [count/total_count for count in counts]
    plt.plot(scores, fractions, label=base_name.split('_')[0], alpha=0.7)

    
def main():
    # plot_file('1234A-V3LOOP_S1_L001_R1_001.fastq')
    plot_file('63899A-V3-2-V3LOOP_S14_L001_R1_001.fastq')
    plot_file('63901A-V3-1-V3LOOP_S25_L001_R1_001.fastq')
    plot_file('73053BMIDI-MidHCV_S7_L001_R1_001.fastq')
    plot_file('73092A-HCV_S69_L001_R1_001.fastq')
    plt.plot((CURRENT_MINIMUM_SCORE, CURRENT_MINIMUM_SCORE), (0, 1), label='current minimum')
    plt.plot((NEW_MINIMUM_SCORE, NEW_MINIMUM_SCORE), (0, 1), label='new minimum')
    plt.legend()
    plt.title('Comparing V3LOOP alignment scores')
    plt.xlabel('Gotoh alignment score')
    plt.ylabel('Fraction of reads')
    plt.savefig('docs/images/v3loop_alignment_scores.png')

main()
