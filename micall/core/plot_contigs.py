import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictReader
from io import StringIO
from itertools import groupby
from math import log10, copysign
from operator import itemgetter, attrgetter
from pathlib import Path

import yaml
from genetracks import Figure, Track, Multitrack, Coverage
# noinspection PyPep8Naming
import drawSvg as draw
from genetracks.elements import Element
from matplotlib import cm, colors
from matplotlib.colors import Normalize

from micall.core.project_config import ProjectConfig
from micall.utils.alignment_wrapper import align_nucs


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
                                        self.h//2-1,
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


class Arrow(Element):
    def __init__(self, start, end, h=20, elevation=0, label=None):
        x = start
        w = end-start
        self.direction = copysign(1, w)
        if w < 0:
            x = end
            w = -w
        super().__init__(x=x, y=0, w=w, h=h)
        self.elevation = elevation
        self.label = label

    def __repr__(self):
        if self.direction >= 0:
            start = self.x
            end = self.x + self.w
        else:
            end = self.x
            start = self.x + self.w
        return f'Arrow({start}, {end}, label={self.label!r})'

    def draw(self, x=0, y=0, xscale=1.0):
        h = self.h
        a = self.x * xscale
        b = (self.x + self.w) * xscale
        x = x * xscale
        r = h/2
        font_size = h * 0.55
        arrow_size = 7
        if self.direction >= 0:
            line_start = a
            arrow_end = b
            arrow_start = max(arrow_end-arrow_size, line_start)
        else:
            line_start = b
            arrow_end = a
            arrow_start = min(arrow_end+arrow_size, line_start)
        centre = (a + b)/2
        arrow_y = h/2 + self.elevation*r
        group = draw.Group(transform="translate({} {})".format(x, y))
        group.append(draw.Line(line_start, arrow_y,
                               arrow_start, arrow_y,
                               stroke='black'))
        group.append(draw.Circle(centre, h/2, r, fill='ivory', stroke='black'))
        group.append(draw.Lines(arrow_end, arrow_y,
                                arrow_start, arrow_y + arrow_size/2,
                                arrow_start, arrow_y - arrow_size/2,
                                arrow_end, arrow_y,
                                fill='black'))
        group.append(draw.Text(self.label,
                               font_size,
                               centre, h / 2,
                               text_anchor='middle',
                               dy="0.35em"))
        return group


class ArrowGroup(Element):
    def __init__(self, arrows: typing.Sequence[Arrow], gap=3):
        arrows = sorted(arrows, key=attrgetter('x', 'w'))
        self.arrows = []
        self.y_coordinates = []
        x_coordinates = []  # [(start, end, index)]
        max_x = None
        for i, arrow in enumerate(arrows):
            self.arrows.append(arrow)
            self.y_coordinates.append(0)
            x2 = arrow.x + arrow.w
            x_coordinates.append((arrow.x, x2, i))
            if arrow.w < 100:
                # Extra padding for label.
                x2 += 100
            if max_x is None:
                max_x = x2
            else:
                max_x = max(max_x, x2)
        h = 0
        while x_coordinates:
            if h > 0:
                h += gap
            row_end = x_coordinates[0][0]-1
            row_height = 0
            for key in x_coordinates[:]:
                x1, x2, i = key
                if x1 < row_end:
                    continue
                x_coordinates.remove(key)
                arrow = self.arrows[i]
                arrow_h = arrow.h
                self.y_coordinates[i] = h + arrow_h
                row_height = max(row_height, arrow_h)
                row_end = x2
            h += row_height
        self.y_coordinates = [y-h for y in self.y_coordinates]
        super().__init__(0, 0, w=max_x, h=h)

    def draw(self, x=0, y=0, xscale=1.0):
        group = draw.Group(transform="translate({} {})".format(x, y))
        for i, (child_y, arrow) in enumerate(zip(self.y_coordinates,
                                                 self.arrows)):
            group.append(arrow.draw(y=child_y, xscale=xscale))
        return group


class ContigMatcher:
    def __init__(self, contig_name):
        self.name = contig_name
        self.num, self.ref = contig_name.split('-', 1)
        if self.num == 'contig':
            self.num, self.ref = self.ref.split('-', 1)
        try:
            int(self.num)
        except ValueError:
            self.num = None
            self.ref = self.name

    def is_match(self, row):
        if self.num is None:
            return False
        row_contig = row.get('contig')
        if row_contig is not None:
            if row_contig != self.name:
                return False
        else:
            if row['contig_num'] != self.num:
                return False
            if row['ref_name'] != self.ref:
                return False
        return True


def plot_genome_coverage(genome_coverage_csv,
                         blast_csv,
                         genome_coverage_svg_path):
    f = build_coverage_figure(genome_coverage_csv, blast_csv)
    f.show(w=970).saveSvg(genome_coverage_svg_path)


def build_coverage_figure(genome_coverage_csv, blast_csv=None):
    min_position, max_position = 1, 500
    coordinate_depths = Counter()
    contig_depths = Counter()
    contig_groups = defaultdict(set)  # {coordinates_name: {contig_name}}
    reader = DictReader(genome_coverage_csv)
    for row in reader:
        query_nuc_pos = int(row['query_nuc_pos'])
        if row['refseq_nuc_pos']:
            refseq_nuc_pos = int(row['refseq_nuc_pos'])
        else:
            refseq_nuc_pos = min_position
        min_position = min(min_position, refseq_nuc_pos, query_nuc_pos)
        max_position = max(max_position, refseq_nuc_pos, query_nuc_pos)
        coordinates_name = row['coordinates']
        contig_name = row['contig']
        if row['coverage'] != '':
            row_coverage = int(row['coverage']) - int(row['dels'])
            coordinate_depths[coordinates_name] = max(
                coordinate_depths[coordinates_name],
                row_coverage)
            contig_depths[contig_name] = max(contig_depths[contig_name],
                                             row_coverage)
        contig_groups[coordinates_name].add(contig_name)
    if '' in coordinate_depths:
        # Force partial contigs to come last.
        coordinate_depths[''] = -1
    position_offset = -min_position + 1
    max_position += position_offset

    blast_rows = []
    if blast_csv is not None:
        for blast_row in DictReader(blast_csv):
            for field_name in ('start', 'end', 'ref_start', 'ref_end'):
                # noinspection PyTypeChecker
                blast_row[field_name] = int(blast_row[field_name])
            blast_rows.append(blast_row)
    blast_rows.sort(key=itemgetter('start', 'ref_start'))

    landmarks_path = (Path(__file__).parent.parent / "data" /
                      "landmark_references.yaml")
    landmark_groups = yaml.safe_load(landmarks_path.read_text())
    projects = ProjectConfig.loadDefault()
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
                    landmark_colour = landmark.get('colour')
                    if landmark_colour is None:
                        continue
                    subtracks.append(Track(landmark['start'] + position_offset,
                                           landmark['end'] + position_offset,
                                           label=landmark['name'],
                                           color=landmark_colour))
                    max_position = max(max_position,
                                       landmark['end'] + position_offset)
                f.add(Multitrack(subtracks))
            break
        else:
            add_partial_banner(f, position_offset, max_position)
        contig_names = contig_groups[coordinates_name]
        sorted_contig_names = sort_contig_names(contig_names, contig_depths)
        ref_arrows = []
        for contig_name in sorted_contig_names:
            if contig_name.startswith('contig-'):
                # No arrows on original contig tracks.
                continue
            contig_matcher = ContigMatcher(contig_name)
            ref_positions = None
            arrow_count = 0
            for blast_row in blast_rows:
                if not contig_matcher.is_match(blast_row):
                    continue
                if (ref_positions is None and
                        coordinates_name != '' and
                        blast_row['ref_name'] != coordinates_name):
                    ref_positions = map_references(blast_row['ref_name'],
                                                   coordinates_name,
                                                   projects)
                arrow_count += 1
                ref_start = int(blast_row['ref_start'])
                ref_end = int(blast_row['ref_end'])
                if ref_positions is None:
                    coordinate_start = ref_start
                    coordinate_end = ref_end
                else:
                    coordinate_start = ref_positions[ref_start]
                    coordinate_end = ref_positions[ref_end]
                ref_arrows.append(
                    Arrow(coordinate_start+position_offset,
                          coordinate_end+position_offset,
                          elevation=1,
                          label=f'{contig_matcher.num}.{arrow_count}'))
        if ref_arrows:
            f.add(ArrowGroup(ref_arrows))
        for contig_name in sorted_contig_names:
            genome_coverage_csv.seek(0)
            reader = DictReader(genome_coverage_csv)
            build_contig(reader,
                         f,
                         contig_name,
                         max_position,
                         position_offset,
                         blast_rows)

    if not f.elements:
        f.add(Track(1, max_position, label='No contigs found.', color='none'))
    return f


def map_references(contig_ref_name: str,
                   coordinates_name: str,
                   projects: ProjectConfig) -> typing.Mapping[int, int]:
    ref_seq = projects.getReference(contig_ref_name)
    coordinates_seq = projects.getReference(coordinates_name)
    aligned_coordinates, aligned_ref, _ = align_nucs(coordinates_seq, ref_seq)
    mapped_positions = {}
    coordinate_pos = ref_pos = 0
    for coordinate_nuc, ref_nuc in zip(aligned_coordinates, aligned_ref):
        if coordinate_nuc != '-':
            coordinate_pos += 1
        if ref_nuc != '-':
            ref_pos += 1
            mapped_positions[ref_pos] = coordinate_pos
    return mapped_positions


def sort_contig_names(contig_names, contig_depths):
    unused_names = set(contig_names)
    sorted_contig_names = [
        contig_name
        for _, contig_name in sorted(
            (-contig_depths[name], name)
            for name in contig_names)]
    final_contig_names = []
    for main_name in sorted_contig_names:
        if main_name not in unused_names:
            continue
        final_contig_names.append(main_name)
        contig_nums, ref_name = main_name.split('-', 1)
        for contig_num in contig_nums.split('_'):
            contig_name = f'contig-{contig_num}-{ref_name}'
            try:
                unused_names.remove(contig_name)
                final_contig_names.append(contig_name)
            except KeyError:
                pass
    return final_contig_names


def build_contig(reader,
                 f,
                 contig_name,
                 max_position,
                 position_offset,
                 blast_rows):
    contig_matcher = ContigMatcher(contig_name)
    blast_ranges = []  # [[start, end, blast_num]]
    blast_starts = defaultdict(set)  # {start: {blast_num}}
    blast_ends = defaultdict(set)  # {end: {blast_num}}
    if not contig_name.startswith('contig-'):
        for blast_row in blast_rows:
            if not contig_matcher.is_match(blast_row):
                continue
            blast_num = len(blast_ranges) + 1
            blast_ranges.append([None, None, blast_num])
            blast_starts[blast_row['start']].add(blast_num)
            blast_ends[blast_row['end']].add(blast_num)
    event_positions = set(blast_starts)
    event_positions.update(blast_ends)
    event_positions = sorted(event_positions, reverse=True)

    insertion_size = 0
    insertion_ranges = []  # [(start, end)]
    unmatched_ranges = []  # [[start, end]]
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
            for field_name in (pos_field, 'coverage', 'dels'):
                field_text = contig_row[field_name]
                field_value = None if field_text == '' else int(field_text)
                contig_row[field_name] = field_value
        start = contig_rows[0][pos_field]
        end = contig_rows[-1][pos_field]
        coverage = [0] * (end - start + 1)
        pos = 0
        for contig_row in contig_rows:
            pos = contig_row[pos_field]
            if pos is None:
                insertion_size += 1
            else:
                if insertion_size:
                    insertion_ranges.append((pos, pos+insertion_size-1))
                    insertion_size = 0
                if contig_row['coverage'] is not None:
                    coverage[pos - start] = (contig_row['coverage'] -
                                             contig_row['dels'])
                contig_pos = int(contig_row['query_nuc_pos'])
                while event_positions and event_positions[-1] <= contig_pos:
                    event_pos = event_positions.pop()
                    for blast_num in blast_starts[event_pos]:
                        blast_ranges[blast_num-1][0] = pos
                    for blast_num in blast_ends[event_pos]:
                        blast_ranges[blast_num-1][1] = pos
            link = contig_row.get('link')
            if link == 'U':
                # Position is unmatched, add to list.
                if not unmatched_ranges or unmatched_ranges[-1][-1] != pos-1:
                    unmatched_ranges.append([pos, pos])
                else:
                    unmatched_ranges[-1][-1] = pos
        while event_positions:
            # Use up any events that went past the end of the contig.
            event_pos = event_positions.pop()
            for blast_num in blast_starts[event_pos]:
                blast_ranges[blast_num - 1][0] = pos
            for blast_num in blast_ends[event_pos]:
                blast_ranges[blast_num - 1][1] = pos

        arrows = []
        for arrow_start, arrow_end, blast_num in blast_ranges:
            arrows.append(Arrow(arrow_start+position_offset,
                                arrow_end+position_offset,
                                elevation=-1,
                                label=f'{contig_matcher.num}.{blast_num}'))
        if arrows:
            f.add(ArrowGroup(arrows))
        subtracks = []
        for has_coverage, group_positions in groupby(
                enumerate(coverage),
                lambda item: item[1] != 0):
            if has_coverage:
                group_positions = list(group_positions)
                group_start, _ = group_positions[0]
                group_end, _ = group_positions[-1]
                subtracks.append(Track(start + group_start + position_offset,
                                       start + group_end + position_offset))
        if not subtracks:
            group_start = prev_pos = None
            included_positions = [row[pos_field] for row in contig_rows]
            included_positions.append(None)  # Trigger final section.
            for pos in included_positions:
                if group_start is None:
                    group_start = pos
                else:
                    if pos != prev_pos + 1:
                        subtracks.append(Track(group_start + position_offset,
                                               prev_pos + position_offset))
                        group_start = pos
                prev_pos = pos
        if max(coverage) <= 0:
            track_label = contig_name
        else:
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
                                        for a, b in insertion_ranges] +
                                       [(a+position_offset,
                                         b+position_offset,
                                         'yellow')
                                        for a, b in unmatched_ranges]))
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
        spans = getattr(track, 'arrows', None)
        if spans is None:
            spans = getattr(track, 'tracks', [track])
        else:
            spans.sort(key=attrgetter('x', 'w', 'label'))
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
            direction = getattr(span, 'direction', None)
            if direction is not None and direction != '':
                if direction >= 0:
                    summary.write(f'{span.x}--{span.label}->{span.x+span.w}')
                else:
                    summary.write(f'{span.x}<-{span.label}--{span.x+span.w}')
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
    parser.add_argument('--blast',
                        help='CSV file with BLAST results for each contig',
                        type=FileType())
    parser.add_argument('genome_coverage_csv',
                        help='CSV file with coverage counts for each contig',
                        type=FileType())
    parser.add_argument('genome_coverage_svg',
                        help='SVG file to plot coverage counts for each contig')
    args = parser.parse_args()

    plot_genome_coverage(args.genome_coverage_csv,
                         args.blast,
                         args.genome_coverage_svg)
    print('Wrote', args.genome_coverage_svg)


if __name__ == '__main__':
    main()
