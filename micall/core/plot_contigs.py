from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictReader
from io import StringIO
from itertools import groupby
from math import log10
from operator import itemgetter
from pathlib import Path

import yaml
from genetracks import Figure, Track, Multitrack, Label, Coverage
# noinspection PyPep8Naming
import drawSvg as draw
from matplotlib import cm, colors
from matplotlib.colors import Normalize


class SmoothCoverage(Coverage):
    def __init__(self, a, b, ys, height=10, color='blue', opacity='1.0'):
        groups = []
        group_y = None
        group_size = 0
        for y in ys + [0]:
            if group_size != 0 and not group_y*0.9 <= y <= group_y*1.1:
                groups.append((group_y, group_size))
                group_size = 0
            if group_size == 0:
                group_y = y
            group_size += 1
        self.coverage_groups = groups
        ys = None
        super(SmoothCoverage, self).__init__(a, b, ys, height, color, opacity)

    def draw(self, x=0, y=0, xscale=1.0):
        a = self.a * xscale
        x = x * xscale
        d = draw.Group(transform="translate({} {})".format(x, y))
        yscale = self.h / max(y for y, count in self.coverage_groups)
        pos = 0
        for y, count in self.coverage_groups:
            if y != 0:
                d.append(draw.Rectangle(a+(pos*xscale),
                                        0,
                                        count*xscale,
                                        y*yscale,
                                        fill=self.get_color(y),
                                        fill_opacity=self.opacity,
                                        shape_rendering='crispEdges'))
            pos += count
        return d

    def get_color(self, coverage):
        return self.color


class ShadedCoverage(SmoothCoverage):
    def __init__(self, a, b, ys, height=10, color='blue', opacity='1.0'):
        super().__init__(a, b, ys, height, color, opacity)
        # noinspection PyUnresolvedReferences
        self.cm = cm.viridis_r
        self.normalize = Normalize(0, 6)

    def get_color(self, coverage):
        log_coverage = log10(coverage)
        rgba = self.cm(self.normalize(log_coverage))

        return colors.to_hex(rgba)


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
    coordinate_depths = Counter()
    contig_depths = Counter()
    contig_groups = defaultdict(set)
    reader = DictReader(contig_coverage_csv)
    for row in reader:
        query_nuc_pos = int(row['query_nuc_pos'])
        if row['refseq_nuc_pos']:
            refseq_nuc_pos = int(row['refseq_nuc_pos'])
        else:
            refseq_nuc_pos = min_position
        min_position = min(min_position, refseq_nuc_pos, query_nuc_pos)
        max_position = max(max_position, refseq_nuc_pos, query_nuc_pos)
        coordinates_name = row['coordinates']
        row_coverage = int(row['coverage'])
        coordinate_depths[coordinates_name] = max(
            coordinate_depths[coordinates_name],
            row_coverage)
        contig_name = row['contig']
        contig_depths[contig_name] = max(contig_depths[contig_name],
                                         row_coverage)
        contig_groups[coordinates_name].add(contig_name)
    if '' in coordinate_depths:
        # Force partial contigs to come last.
        coordinate_depths[''] = -1
    position_offset = -min_position + 1
    max_position += position_offset

    landmarks_path = (Path(__file__).parent.parent / "data" /
                      "landmark_references.yaml")
    landmark_groups = yaml.safe_load(landmarks_path.read_text())
    f = Figure()
    for _, coordinates_name in sorted((-depth, name)
                                      for name, depth in coordinate_depths.items()):
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
            add_partial_banner(f, position_offset, max_position)
        for _, contig_name in sorted((-contig_depths[name], name)
                                     for name in contig_groups[coordinates_name]):
            contig_coverage_csv.seek(0)
            reader = DictReader(contig_coverage_csv)
            build_contig(reader, f, contig_name, max_position, position_offset)

    if not f.elements:
        f.add(Track(1, max_position, label='No contigs found.', color='none'))
    return f


def build_contig(reader, f, contig_name, max_position, position_offset):
    insertion_size = 0
    insertion_ranges = []  # [(start, end)]
    for contig_name2, contig_rows in groupby(reader, itemgetter('contig')):
        if contig_name2 != contig_name:
            continue
        contig_rows = list(contig_rows)
        coordinates_name = contig_rows[0]['coordinates']
        if coordinates_name:
            pos_field = 'refseq_nuc_pos'
        else:
            pos_field = 'query_nuc_pos'
        for contig_row in contig_rows:
            for field_name in (pos_field, 'coverage'):
                field_text = contig_row[field_name]
                field_value = None if field_text == '' else int(field_text)
                contig_row[field_name] = field_value
        start = contig_rows[0][pos_field]
        end = contig_rows[-1][pos_field]
        coverage = [0] * (end - start + 1)
        for contig_row in contig_rows:
            pos = contig_row[pos_field]
            if pos is None:
                insertion_size += 1
            else:
                if insertion_size:
                    insertion_ranges.append((pos, pos+insertion_size-1))
                    insertion_size = 0
                coverage[pos - start] = contig_row['coverage']
        subtracks = []
        for has_coverage, group_positions in groupby(
                enumerate(coverage, 0),
                lambda item: item[1] != 0):
            if has_coverage:
                group_positions = list(group_positions)
                group_start, _ = group_positions[0]
                group_end, _ = group_positions[-1]
                subtracks.append(Track(start + group_start + position_offset,
                                       start + group_end + position_offset))
        if max(coverage) > 0:
            f.add(ShadedCoverage(start + position_offset,
                                 end + position_offset,
                                 coverage),
                  gap=-4)
        track_label = f"{contig_name} - depth {max(coverage)}"
        subtracks.append(Track(1,
                               max_position,
                               label=track_label,
                               color='none',
                               regions=[(a+position_offset,
                                         b+position_offset,
                                         'lightgreen')
                                        for a, b in insertion_ranges]))
        f.add(Multitrack(subtracks))
        break


def add_partial_banner(f, position_offset, max_position):
    """ Build a dashed line with dashes 500 wide. """
    dash_width = 500
    banner_width = max_position - position_offset
    subtracks = [Track(i*dash_width + position_offset + 1,
                       min((i+1)*dash_width + position_offset, max_position))
                 for i in range((banner_width + dash_width) // dash_width)
                 if not i % 2]
    subtracks.append(Track(position_offset+1,
                           max_position,
                           label='Partial Blast Results',
                           color='none'))
    f.add(Multitrack(subtracks))


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
            coverage_groups = getattr(span, 'coverage_groups', None)
            if coverage_groups is not None:
                summary.write('Coverage ')
                for j, (y, count) in enumerate(coverage_groups):
                    if j:
                        summary.write(', ')
                    summary.write(str(y))
                    if count > 1:
                        summary.write(f'x{count}')
                continue
            span_text = getattr(span.label, 'text', span.label) or ''
            summary.write(span_text)
            color = getattr(span, 'color')
            if span.a or span.b:
                if color != 'none':
                    summary.write(f'[{span.a}-{span.b}]')
                else:
                    summary.write(f'({span.a}-{span.b})')
            regions = getattr(span, 'regions', [])
            for start, end, colour in regions:
                summary.write(f', {colour}{{{start}-{end}}}')
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
