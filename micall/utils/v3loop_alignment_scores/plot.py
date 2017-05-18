from csv import DictReader
from operator import itemgetter

from matplotlib import pyplot as plt

CURRENT_MINIMUM_SCORE = 52
NEW_MINIMUM_SCORE = 200

def plot_file(filename):
    with open(filename) as f:
        reader = DictReader(f)
        rows = [list(map(int, row))
                for row in map(itemgetter('score', 'count'), reader)]
        scores = [row[0] for row in rows]
        counts = [row[1] for row in rows]
        total_count = sum(counts)
        fractions = [count/total_count for count in counts]
        plt.plot(scores, fractions, label=filename.split('_')[0], alpha=0.7)

plot_file('63899A-V3-2-V3LOOP_S14_v3loop_scores.csv')
plot_file('63901A-V3-1-V3LOOP_S25_v3loop_scores.csv')
plot_file('73053BMIDI-MidHCV_S7_v3loop_scores.csv')
plot_file('73092A-HCV_S69_v3loop_scores.csv')
plt.plot((CURRENT_MINIMUM_SCORE, CURRENT_MINIMUM_SCORE), (0, 1), label='current minimum')
plt.plot((NEW_MINIMUM_SCORE, NEW_MINIMUM_SCORE), (0, 1), label='new minimum')
plt.legend()
plt.title('Comparing V3LOOP alignment scores')
plt.xlabel('Gotoh alignment score')
plt.ylabel('Fraction of reads')
plt.savefig('../../../docs/images/v3loop_alignment_scores.png')

