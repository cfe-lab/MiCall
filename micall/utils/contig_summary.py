from argparse import ArgumentParser
from csv import DictReader
from io import StringIO
from pathlib import Path

from micall.utils.fasta_to_csv import default_database
from micall.utils.externals import Blastn

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa

DEFAULT_SCRATCH_PATH = (Path(__file__).parent.parent / "tests" / "working" /
                        "basespace_linked_data" / "scratch")


def parse_args():
    parser = ArgumentParser(
        description='Scan contig files, and write a summary.')
    parser.add_argument('scratch',
                        nargs='?',
                        type=Path,
                        default=DEFAULT_SCRATCH_PATH)
    return parser.parse_args()


def main():
    args = parse_args()
    sample_dirs = [d for d in args.scratch.iterdir() if d.is_dir()]
    contig_plots_path = sample_dirs[0].parent / 'contig_plots'
    contig_plots_path.mkdir(exist_ok=True)
    sample_dirs.sort()
    incomplete_count = empty_count = 0
    for sample_dir in sample_dirs:
        if sample_dir == contig_plots_path:
            continue
        contigs_path = sample_dir / 'contigs.csv'
        contigs_size = contigs_path.stat().st_size
        if contigs_size == 0:
            incomplete_count += 1
        elif contigs_size == 22:
            empty_count += 1
        else:
            contigs_fasta_paths = list(
                sample_dir.glob('assembly_*/iva/contigs.fasta'))
            contigs_fasta_paths += list(
                sample_dir.glob('assembly_*/contigs_stage_c.fasta'))
            if len(contigs_fasta_paths) != 1:
                print(sample_dir, contigs_fasta_paths)
                continue
            contigs_fasta_path, = contigs_fasta_paths
            with default_database() as DEFAULT_DATABASE:
                stdout = Blastn().genotype(
                    contigs_fasta=contigs_fasta_path,
                    database=DEFAULT_DATABASE,
                )
            plot_contigs(sample_dir, stdout)
            plot_path = contig_plots_path / (sample_dir.name + '.png')
            plt.savefig(str(plot_path))
            plt.close()

    print(empty_count, 'empty')
    print(incomplete_count, 'incomplete')
    print('=' * 20)
    print(len(sample_dirs), 'Total')


def plot_contigs(sample_dir, contigs_csv):
    reader = DictReader(StringIO(contigs_csv), Blastn.BLAST_COLUMNS)
    rows = list(reader)
    rows.sort(reverse=True, key=lambda row: int(row['score']))
    contig_names = sorted({row['qaccver'] for row in rows})
    if not contig_names:
        plt.title('No contigs')
        plt.yticks([])
        plt.xticks([])
        return
    y_contig = 0
    # noinspection PyTypeChecker
    fig, axes_list = plt.subplots(nrows=len(contig_names), sharex=True, squeeze=False)
    for i, contig_name in enumerate(contig_names):
        contig_rows = [row for row in rows if row['qaccver'] == contig_name]
        ax = axes_list[i][0]
        if i == 0:
            ax.set_title(sample_dir.name)
        ref_name = contig_rows[0]['saccver']
        gt_rows = [row for row in contig_rows if row['saccver'] == ref_name]
        ax.set_yticks([y_contig, len(gt_rows)])
        ax.set_yticklabels([contig_name, ref_name])
        ax.plot([1, int(contig_rows[0]['qlen'])], [y_contig, y_contig], 'k')
        for match_num, row in enumerate(gt_rows):
            match_y = len(gt_rows) - match_num
            contig_start = int(row['qstart'])
            ref_start = int(row['sstart'])
            contig_end = int(row['qend'])
            ref_end = int(row['send'])
            start_format, end_format = 'gr'
            if ref_end < ref_start:
                start_format += ':'
                end_format += ':'
            ax.plot([contig_start, ref_start], [y_contig, match_y], start_format)
            ax.plot([contig_end, ref_end], [y_contig, match_y], end_format)
            ax.plot([ref_start, ref_end], [match_y, match_y], 'k')
    plt.tight_layout()


def live_plot():
    # qaccver,saccver,pident,score,qcovhsp,qstart,qend,sstart,send
    contigs_csv = """\
contig.00001,HCV-3a,9405,91.717,4747,75,2291,9389,1523,8621
contig.00001,HCV-3a,9405,93.964,930,13,990,2215,222,1447
contig.00001,HCV-3a,9405,100.000,58,1,172,229,223,280
contig.00001,HCV-3a,9405,100.000,58,1,228,285,223,280
contig.00001,HCV-3a,9405,100.000,55,1,284,338,223,277
contig.00001,HCV-3a,9405,100.000,38,0,136,173,243,280
contig.00001,HCV-3a,9405,100.000,32,0,960,991,221,252
contig.00001,HCV-2q,9405,100.000,31,0,40,70,258,288
contig.00001,HCV-2q,9405,100.000,31,0,72,102,258,288
contig.00001,HCV-2q,9405,100.000,31,0,104,134,258,288
contig.00001,HCV-2q,9405,100.000,28,0,9362,9389,8645,8672
contig.00001,HCV-2c,9405,84.679,235,6,990,1596,264,870
contig.00001,HCV-2c,9405,98.276,54,1,172,229,265,322
contig.00001,HCV-2c,9405,98.276,54,1,228,285,265,322
contig.00002,HCV-6k,812,100.000,31,4,1,31,292,262
contig.00002,HCV-6k,812,100.000,29,4,68,96,262,290
contig.00002,HCV-6k,812,100.000,29,4,100,128,262,290
contig.00002,HCV-6k,812,100.000,29,4,197,225,290,262
contig.00002,HCV-6e,812,100.000,31,4,1,31,292,262
contig.00002,HCV-6e,812,100.000,29,4,68,96,262,290
contig.00002,HCV-6e,812,100.000,29,4,100,128,262,290
contig.00003,HCV-2i,658,100.000,32,5,3,34,8649,8680
contig.00003,HCV-2i,658,100.000,32,5,62,93,8649,8680
contig.00003,HCV-2i,658,100.000,32,5,505,536,8649,8680
contig.00003,HCV-2i,658,100.000,32,5,564,595,8649,8680
contig.00003,HCV-2i,658,100.000,32,5,623,654,8649,8680
contig.00003,HCV-6k,658,100.000,29,4,447,475,290,262
contig.00003,HCV-6k,658,100.000,29,4,415,443,290,262
contig.00003,HCV-6k,658,100.000,29,4,383,411,290,262
contig.00003,HCV-6k,658,100.000,29,4,351,379,290,262
"""
    plot_contigs(Path('/tmp/1234A-HCV'), contigs_csv)
    plt.show()


if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    live_plot()
