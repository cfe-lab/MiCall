from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter
from csv import DictReader
from io import StringIO
from itertools import groupby
from operator import itemgetter
from pathlib import Path

import yaml
from genetracks import Figure, Track, Multitrack, Label, Coverage


def hcv_region_tracks(f):
    # H77 coords
    first = [
            (342, 915, "C", "salmon"),
            (915, 1491, "E1", "red"),
            (1491, 2580, "E2", "plum"),
            (2580, 2769, "P7", "yellowgreen"),
            (2769, 3420, "NS2", "green"),
            (3420, 5313, "NS3", "turquoise"),
            (5313, 5475, "NS4A", "yellowgreen"),
            (5475, 6258, "NS4B", "darkgrey"),
            (6258, 7602, "NS5A", "lightblue"),
            (7602, 9378, "NS5B", "steelblue"),
            ]
    for reading_frame in [first]:
        f.add(Multitrack([Track(l,
                                r,
                                label=Label(0, text, offset=1),
                                color=color)
                          for l, r, text, color in reading_frame]),
              gap=0)


def hiv_region_tracks(f):
    third = [
        (2085, 5096, "pol", "orange"),
        (5559, 5850, "vpr", "turquoise"),
        (5970, 6045, "rev", 'yellowgreen'),
        (6225, 8795, "env", 'salmon'),
    ]
    second = [
        (5831, 6045, "tat", "plum"),
        (6062, 6310, "vpu", "red"),
        (8379, 8653, "rev", 'yellowgreen'),
        (9086, 9719, "3' LTR", 'darkgrey'),
    ]

    first = [
        (0, 634, "5' LTR", "darkgrey"),
        (790, 2292, "gag", "lightblue"),
        (5041, 5619, "vif", 'steelblue'),
        (8379, 8469, "tat", 'plum'),
        (8797, 9417, "nef", 'mediumaquamarine'),
    ]

    for reading_frame in [first, second, third]:
        f.add(Multitrack([Track(l, r, label=Label(0, text, offset=1), color=color)
                          for l, r, text, color in reading_frame]),
              gap=0)


def plot_contigs(blast_csv, blast_svg_path):
    blast_rows = DictReader(blast_csv)
    f = build_contigs_figure(blast_rows, centre_x=4000)
    f.show(w=970).saveSvg(blast_svg_path)


def build_blast_key(org_scores, ref_scores, row):
    ref_name = row['ref_name']
    ref_score = ref_scores[ref_name]
    org_name = row['organism']
    org_score = org_scores[org_name]
    return -org_score, org_name, -ref_score, ref_name, row['ref_start']


def build_contigs_figure(blast_rows, limit=10, centre_x=0):
    all_rows = list(blast_rows)
    org_scores = Counter()
    ref_scores = Counter()
    for row in all_rows:
        row['match'] = float(row['match'])
        row['organism'] = row['ref_name'][:3]
        for field in ('score', 'start', 'end', 'ref_start', 'ref_end'):
            field_value = row[field]
            row[field] = int(field_value) if field_value else None
        org_scores[row['organism']] += row['score']
        ref_scores[row['ref_name']] += row['score']
    header_writers = {'HIV': hiv_region_tracks,
                      'HCV': hcv_region_tracks}

    all_rows.sort(key=itemgetter('score'), reverse=True)
    all_rows = all_rows[:limit]
    all_rows.sort(key=lambda r: build_blast_key(org_scores, ref_scores, r))
    f = Figure()
    for organism, org_rows in groupby(all_rows, itemgetter('organism')):
        header_writer = header_writers.get(organism)
        if not header_writer:
            continue
        header_writer(f)
        for ref_name, ref_rows in groupby(org_rows, itemgetter('ref_name')):
            f.add(Track(centre_x, centre_x, label=ref_name + ':'))
            for row in ref_rows:
                f.add(Track(row['ref_start'],
                            row['ref_end'],
                            label=row['contig_name']))
    if not f.elements:
        f.add(Track(centre_x, centre_x, label='No contigs assembled.'))
    return f


def plot_contig_coverage(contig_coverage_csv, contig_coverage_svg_path):
    f = build_coverage_figure(contig_coverage_csv)
    f.show(w=970).saveSvg(contig_coverage_svg_path)


def build_coverage_figure(contig_coverage_csv):
    min_position, max_position = 1, 500
    contig_headers = set()
    reader = DictReader(contig_coverage_csv)
    for row in reader:
        query_nuc_pos = int(row['query_nuc_pos'])
        if row['refseq_nuc_pos']:
            refseq_nuc_pos = int(row['refseq_nuc_pos'])
        else:
            refseq_nuc_pos = min_position
        min_position = min(min_position, refseq_nuc_pos, query_nuc_pos)
        max_position = max(max_position, refseq_nuc_pos, query_nuc_pos)
        contig_headers.add((row['coordinates'] == '', row['coordinates']))
    position_offset = -min_position + 1

    landmarks_path = (Path(__file__).parent.parent / "data" /
                      "landmark_references.yaml")
    landmark_groups = yaml.safe_load(landmarks_path.read_text())
    f = Figure()
    for _, coordinates_name in sorted(contig_headers):
        for reference_set in landmark_groups:
            if coordinates_name != reference_set['coordinates']:
                continue
            prev_landmark = None
            for i, landmark in enumerate(sorted(reference_set['landmarks'],
                                                key=itemgetter('start'))):
                landmark.setdefault('frame', 0)
                if prev_landmark and 'end' not in prev_landmark:
                    prev_landmark['end'] = landmark['start'] - 1
                prev_landmark = landmark
            for frame, frame_landmarks in groupby(reference_set['landmarks'],
                                                  itemgetter('frame')):
                subtracks = []
                for landmark in frame_landmarks:
                    subtracks.append(Track(landmark['start']+position_offset,
                                           landmark['end']+position_offset,
                                           label=landmark['name'],
                                           color=landmark['colour']))
                    max_position = max(max_position,
                                       landmark['end'] + position_offset)
                f.add(Multitrack(subtracks))
            break
        else:
            f.add(Track(min_position,
                        max_position,
                        label='Partial Blast Results',
                        color='none'))
        contig_coverage_csv.seek(0)
        reader = DictReader(contig_coverage_csv)
        for contig_name, contig_rows in groupby(reader, itemgetter('contig')):
            contig_rows = list(contig_rows)
            if contig_rows[0]['coordinates'] != coordinates_name:
                continue
            if coordinates_name:
                pos_field = 'refseq_nuc_pos'
            else:
                pos_field = 'query_nuc_pos'
            for contig_row in contig_rows:
                for field_name in (pos_field, 'coverage'):
                    contig_row[field_name] = int(contig_row[field_name])
            start = contig_rows[0][pos_field]
            end = contig_rows[-1][pos_field]
            coverage = [0] * (end-start+1)
            for contig_row in contig_rows:
                coverage[contig_row[pos_field] - start] = contig_row['coverage']
            if max(coverage) > 0:
                f.add(Coverage(start+position_offset,
                               end+position_offset,
                               coverage),
                      gap=-4)
            track_label = f"{contig_name} - depth {max(coverage)}"
            f.add(Multitrack([Track(start+position_offset, end+position_offset),
                              Track(1,
                                    max_position,
                                    label=track_label,
                                    color='none')]))

    if not f.elements:
        f.add(Track(1, max_position, label='No contigs found.', color='none'))
    return f


def summarize_figure(figure: Figure):
    """ Summarize the contents of a figure to text.

    Useful for testing.
    """
    figure.show()  # Test that all the display math works.

    summary = StringIO()
    for padding, track in figure.elements:
        spans = getattr(track, 'tracks', [track])
        for i, span in enumerate(spans):
            if i:
                summary.write(', ')
            ys = getattr(span, 'ys', None)
            if ys is not None:
                summary.write('Coverage ')
                summary.write(', '.join(map(str, ys)))
                continue
            span_text = getattr(span.label, 'text', span.label) or ''
            summary.write(span_text)
            color = getattr(span, 'color')
            if span.a or span.b:
                if color != 'none':
                    summary.write(f'[{span.a}-{span.b}]')
                else:
                    summary.write(f'({span.a}-{span.b})')
        summary.write('\n')
    return summary.getvalue()


def main():
    parser = ArgumentParser(
        description='Plot assembled contigs against a reference.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('contig_coverage_csv',
                        help='CSV file with coverage counts for each contig',
                        type=FileType())
    parser.add_argument('contig_coverage_svg',
                        help='SVG file to plot coverage counts for each contig')
    args = parser.parse_args()

    plot_contig_coverage(args.contig_coverage_csv, args.contig_coverage_svg)
    print('Wrote', args.contig_coverage_svg)


if __name__ == '__main__':
    main()
