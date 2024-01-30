import typing
from typing import Dict, Tuple, List, Set, Iterable, NoReturn, Literal
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictReader
from io import StringIO
from itertools import groupby
from math import log10, copysign
from operator import itemgetter, attrgetter
from pathlib import Path
import dataclasses
import logging

import yaml
from genetracks import Figure, Track, Multitrack, Coverage
# noinspection PyPep8Naming
import drawsvg as draw
from genetracks.elements import Element, Label
from matplotlib import cm, colors
from matplotlib.colors import Normalize

from micall.core.project_config import ProjectConfig
from micall.utils.alignment_wrapper import align_nucs
from micall.core.contig_stitcher import Contig, GenotypedContig, AlignedContig, sliding_window
from micall.utils.cigar_tools import CigarHit
import micall.utils.contig_stitcher_events as events


logger = logging.getLogger(__name__)


class LeftLabel(Label):
    """Like Label, but anchored to the left, instead of the middle.
    """
    def draw(self, *args, **kwargs):
        d = super().draw(*args, **kwargs)
        assert len(d.children) == 1
        text = d.children[0]
        text.args['text-anchor'] = 'left'
        # text.args['fill'] = 'red' # works by the way
        return d


class SmoothCoverage(Coverage):
    def __init__(self, a, b, ys, height=10, color='blue', opacity='1.0'):
        groups = []
        group_y = None
        group_size = 0
        for y in ys + [-1]:
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


class ConcordanceLine(SmoothCoverage):
    def __init__(self, a, b, ys):
        self.a = a
        self.b = b
        self.ys = ys
        height = 10
        color = 'red'
        opacity = '1.0'
        super(ConcordanceLine, self).__init__(a, b, ys, height, color, opacity)

    def draw(self, x=0, y=0, xscale=1.0):
        a = self.a * xscale
        x = x * xscale
        d = draw.Group(transform="translate({} {})".format(x, y))
        p = draw.Path(stroke=self.color, stroke_width=1, fill='none')
        yscale = self.h / 100
        pos = 0
        for y, count in self.coverage_groups:
            if pos == 0:
                p.M(a, self.h//2-1 + y*yscale)
                pos += count
            else:
                # place points at the beginning and end of the group
                p.L(a + (pos * xscale), self.h // 2 - 1 + y * yscale)
                pos += count
                p.L(a + (pos * xscale), self.h // 2 - 1 + y * yscale)
        d.append(p)
        r = draw.Rectangle(self.a*xscale,
                           self.h//2-1,
                           (self.b-self.a)*xscale,
                           self.h,
                           fill=self.color,
                           fill_opacity='0.3',
                           shape_rendering='crispEdges')
        d.append(r)
        return d


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
        if self.label is not None:
            group.append(draw.Circle(centre, h/2, r, fill='ivory', stroke='black'))
        group.append(draw.Lines(arrow_end, arrow_y,
                                arrow_start, arrow_y + arrow_size/2,
                                arrow_start, arrow_y - arrow_size/2,
                                arrow_end, arrow_y,
                                fill='black'))
        if self.label is not None:
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
                         genome_coverage_svg_path,
                         use_concordance=False):
    f = build_coverage_figure(genome_coverage_csv, blast_csv, use_concordance)
    f.show(w=970).save_svg(genome_coverage_svg_path, context=draw.Context(invert_y=True))


def build_coverage_figure(genome_coverage_csv, blast_csv=None, use_concordance=False):
    min_position, max_position = 1, 500
    coordinate_depths = Counter()
    contig_depths = Counter()
    contig_groups = defaultdict(set)  # {coordinates_name: {contig_name}}
    reader = DictReader(genome_coverage_csv)
    for row in reader:
        if row['query_nuc_pos']:
            query_nuc_pos = int(row['query_nuc_pos'])
        else:
            query_nuc_pos = min_position
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
                         blast_rows,
                         use_concordance=use_concordance)

    if not f.elements:
        f.add(Track(1, max_position, label='No contigs found.', color='none'))
    return f


def plot_stitcher_coverage(logs: Iterable[events.EventType], genome_coverage_svg_path: str):
    f = build_stitcher_figure(logs)
    f.show(w=970).save_svg(genome_coverage_svg_path, context=draw.Context(invert_y=True))
    return f


def build_stitcher_figure(logs: Iterable[events.EventType]) -> Figure:
    complete_contig_map: Dict[str, GenotypedContig] = {}
    name_map: Dict[str, str] = {}
    complete_parent_graph: Dict[str, List[str]] = {}
    alive_set: Set[str] = set()
    morphism_graph: Dict[str, List[str]] = {}
    reduced_parent_graph: Dict[str, List[str]] = {}
    transitive_parent_graph: Dict[str, List[str]] = {}
    discarded: List[str] = []
    unknown: List[str] = []
    anomaly: List[str] = []
    unaligned: List[str] = []
    overlaps_list: List[str] = []
    overlap_leftparent_map: Dict[str, str] = {}
    overlap_rightparent_map: Dict[str, str] = {}
    overlap_lefttake_map: Dict[str, str] = {}
    overlap_righttake_map: Dict[str, str] = {}
    overlap_sibling_map: Dict[str, str] = {}
    combine_left_edge: Dict[str, str] = {}
    combine_right_edge: Dict[str, str] = {}
    children_join_points: List[str] = []
    query_position_map: Dict[str, Tuple[int, int]] = {}
    lstrip_map: Dict[str, str] = {}
    rstrip_map: Dict[str, str] = {}
    strip_set: Set[Tuple[str, int, int]] = set()

    def remove_intermediate_edges(graph):
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                if all(other not in graph.get(child, []) for other in children):
                    lst.append(child)
            ret[parent] = lst
        return ret

    def remove_transitive_edges(graph):
        tr_cl = transitive_closure(graph)
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                is_transitive = any(child in tr_cl.get(other_node, []) for other_node in children if other_node != child)
                if not is_transitive:
                    lst.append(child)
            ret[parent] = lst
        return ret

    def get_all_ancestors(recur, lst, graph, ancestor_name):
        if ancestor_name not in recur:
            recur = recur.copy()
            recur.add(ancestor_name)

            if ancestor_name not in lst:
                lst.append(ancestor_name)

            existing_ancestors = graph.get(ancestor_name, [])
            for existing in existing_ancestors:
                get_all_ancestors(recur, lst, graph, existing)

    def transitive_closure(graph):
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                get_all_ancestors(set(), lst, graph, child)
            ret[parent] = lst
        return ret

    def copy_graph(graph):
        ret = {}
        for parent, children in graph.items():
            ret[parent] = children[:]
        return ret

    def reflexive_closure(graph):
        ret = copy_graph(graph)
        for parent, children in ret.items():
            if parent not in children:
                children.append(parent)
            for child in children[:]:
                if child not in ret:
                    ret[child] = []
                lst = ret[child]
                if child not in lst:
                    ret[child].append(child)
        return ret

    def inverse_graph(graph):
        ret = {}
        for parent, children in graph.items():
            for child in children:
                if child not in ret:
                    ret[child] = []
                ret[child].append(parent)
        return ret

    def graph_sum(graph_a, graph_b):
        ret = copy_graph(graph_a)
        for key, values in graph_b.items():
            if key not in ret:
                ret[key] = []
            for value in values:
                lst = ret[key]
                if value not in lst:
                    lst.append(value)
        return ret

    def symmetric_closure(graph):
        return graph_sum(graph, inverse_graph(graph))

    def record_unaligned_parts(original: AlignedContig, q_st: int, r_st: int, length: int):
        key = (original.seq, q_st, q_st + length)
        if length > 0 and key not in strip_set:
            strip_set.add(key)
            alignment = CigarHit.from_default_alignment(q_st=q_st, q_ei=q_st + length - 1, r_st=r_st, r_ei=r_st-1)
            seq = 'A' * alignment.query_length
            query = dataclasses.replace(original, name=f"u{len(complete_contig_map)}", seq=seq)
            fake_aligned = AlignedContig.make(query, alignment, strand=original.strand)
            record_contig(fake_aligned, [original])
            record_bad_contig(fake_aligned, unaligned)
            record_alive(fake_aligned)
            return fake_aligned
        return None

    def record_regular_strip(result: AlignedContig, original: AlignedContig):
        length = abs(result.alignment.query_length - original.alignment.query_length)
        q_st = original.alignment.q_st
        r_st = original.alignment.r_st
        return record_unaligned_parts(original, q_st=q_st, r_st=r_st, length=length)

    def record_initial_strip(original: AlignedContig, q_st: int, q_ei: int):
        length = q_ei - q_st + 1
        contig = record_unaligned_parts(original, q_st, original.alignment.r_st, length)
        if contig:
            query_position_map[contig.name] = (q_st, q_ei)

    def record_contig(contig: GenotypedContig, parents: List[GenotypedContig]):
        complete_contig_map[contig.name] = contig
        if [contig.name] != [parent.name for parent in parents]:
            for parent in parents:
                complete_contig_map[parent.name] = parent
                if contig.name not in complete_parent_graph:
                    complete_parent_graph[contig.name] = []

                complete_parent_graph[contig.name].append(parent.name)

    def record_alive(contig: Contig):
        alive_set.add(contig.name)

    def record_bad_contig(contig: GenotypedContig, lst: List[str]):
        complete_contig_map[contig.name] = contig
        lst.append(contig.name)

    def record_lstrip(result: AlignedContig, original: AlignedContig):
        lstrip_map[result.name] = original.name
        unaligned = record_regular_strip(result, original)
        if unaligned:
            lstrip_map[unaligned.name] = result.name

    def record_rstrip(result: AlignedContig, original: AlignedContig):
        rstrip_map[result.name] = original.name
        unaligned = record_regular_strip(result, original)
        if unaligned:
            rstrip_map[unaligned.name] = result.name

    for event in logs:
        if isinstance(event, events.FinalCombine):
            record_contig(event.result, event.contigs)
            record_alive(event.result)
        elif isinstance(event, events.SplitGap):
            record_contig(event.left, [event.contig])
            record_contig(event.right, [event.contig])
            record_alive(event.left)
            record_alive(event.right)
        elif isinstance(event, events.Intro):
            record_contig(event.contig, [])
            record_alive(event.contig)
        elif isinstance(event, events.Hit):
            record_contig(event.part, [event.contig])
            record_alive(event.part)
        elif isinstance(event, events.NoRef):
            record_bad_contig(event.contig, unknown)
            record_alive(event.contig)
        elif isinstance(event, events.ZeroHits):
            record_bad_contig(event.contig, anomaly)
            record_alive(event.contig)
        elif isinstance(event, events.StrandConflict):
            record_bad_contig(event.contig, anomaly)
            record_alive(event.contig)
        elif isinstance(event, events.ReverseComplement):
            record_contig(event.result, [event.contig])
            record_alive(event.result)
        elif isinstance(event, events.HitNumber):
            record_alive(event.contig)
        elif isinstance(event, events.Munge):
            record_contig(event.result, [event.left, event.right])
        elif isinstance(event, events.LStrip):
            record_contig(event.result, [event.original])
            record_lstrip(event.result, event.original)
        elif isinstance(event, events.RStrip):
            record_contig(event.result, [event.original])
            record_rstrip(event.result, event.original)
        elif isinstance(event, events.InitialStrip):
            record_initial_strip(event.contig, event.q_st, event.q_ei)
        elif isinstance(event, events.Overlap):
            overlaps_list.append(event.left_overlap.name)
            overlaps_list.append(event.right_overlap.name)
            overlap_leftparent_map[event.left_remainder.name] = event.left.name
            overlap_rightparent_map[event.right_remainder.name] = event.right.name
            overlap_lefttake_map[event.left_remainder.name] = event.left_take.name
            overlap_righttake_map[event.right_remainder.name] = event.right_take.name
            overlap_sibling_map[event.left_remainder.name] = event.right_remainder.name
            overlap_sibling_map[event.right_remainder.name] = event.left_remainder.name
        elif isinstance(event, events.Drop):
            record_bad_contig(event.contig, discarded)
            record_alive(event.contig)
        elif isinstance(event, events.StitchCut):
            record_contig(event.left_overlap, [event.left])
            record_contig(event.left_remainder, [event.left])
            record_contig(event.right_overlap, [event.right])
            record_contig(event.right_remainder, [event.right])
        elif isinstance(event, events.Stitch):
            record_contig(event.result, [event.left, event.right])
            record_alive(event.result)
        elif isinstance(event, events.Cut):
            record_contig(event.left, [event.original])
            record_contig(event.right, [event.original])
        elif isinstance(event, events.Combine):
            record_alive(event.result)
            record_contig(event.result, event.contigs)
            if event.contigs:
                combine_left_edge[event.result.name] = event.contigs[0].name
                combine_right_edge[event.result.name] = event.contigs[-1].name
        elif isinstance(event, (events.IgnoreGap, events.NoOverlap)):
            pass
        else:
            x: NoReturn = event
            raise RuntimeError(f"Unrecognized action or event: {event}")

    nodup_parent_graph = remove_transitive_edges(complete_parent_graph)

    # Close alive set by parents
    def extend_alive(contig_name):
        alive_set.add(contig_name)
        for parent_name in nodup_parent_graph.get(contig_name, []):
            extend_alive(parent_name)

    for contig_name in alive_set.copy():
        extend_alive(contig_name)

    parent_graph: Dict[str, List[str]] = {}
    for contig_name in nodup_parent_graph:
        if contig_name in alive_set:
            parent_graph[contig_name] = nodup_parent_graph[contig_name]

    contig_map: Dict[str, GenotypedContig] = {k: v for k, v in complete_contig_map.items() if k in alive_set}
    bad_contigs = anomaly + discarded + unknown + unaligned
    group_refs = {contig.group_ref: len(contig.ref_seq) for contig in contig_map.values() if contig.ref_seq}
    children_graph = inverse_graph(parent_graph)
    transitive_parent_graph = transitive_closure(parent_graph)
    transitive_children_graph = transitive_closure(children_graph)
    reduced_parent_graph = remove_intermediate_edges(transitive_parent_graph)
    eqv_parent_graph = reflexive_closure(symmetric_closure(transitive_parent_graph))
    sorted_roots = list(sorted(parent_name for
                               parent_name in contig_map
                               if parent_name not in parent_graph))
    sorted_sinks = list(sorted(child_name for
                               child_name in contig_map
                               if child_name not in children_graph))

    for contig_name, parents in parent_graph.items():
        if len(parents) == 1:
            morphism_graph[parents[0]] = [contig_name]

    transitive_morphism_graph = transitive_closure(morphism_graph)
    reduced_morphism_graph = remove_intermediate_edges(transitive_morphism_graph)
    eqv_morphism_graph = reflexive_closure(symmetric_closure(transitive_morphism_graph))

    for contig_name, parents in parent_graph.items():
        if len(parents) > 1:
            children_join_points.append(contig_name)

    def set_query_position(contig_name: str) -> None:
        contig = contig_map[contig_name]
        children_names = children_graph.get(contig.name, [])

        def copy_from_parent(contig: AlignedContig, parent_name: str) -> None:
            if parent_name in query_position_map:
                (original_q_st, original_q_ei) = query_position_map[parent_name]
                (current_q_st, current_q_ei) = (contig.alignment.q_st, contig.alignment.q_ei)
                current_query_len = abs(current_q_st - current_q_ei)

                if contig_name in lstrip_map:
                    if contig_name in unaligned:
                        query_position_map[contig.name] = (original_q_st - current_query_len - 1, original_q_st - 1)
                    else:
                        query_position_map[contig.name] = (original_q_ei - current_query_len, original_q_ei)
                elif contig_name in rstrip_map:
                    if contig_name in unaligned:
                        query_position_map[contig.name] = (original_q_ei + 1, original_q_ei + 1 + current_query_len)
                    else:
                        query_position_map[contig.name] = (original_q_st, original_q_st + current_query_len)
                else:
                    query_position_map[contig_name] = query_position_map[parent_name]

        if contig_name not in query_position_map:
            if isinstance(contig, AlignedContig):
                regular_parents_names = parent_graph.get(contig_name, [])
                regular_parents_names = [name for name in regular_parents_names if name in query_position_map]
                strip_parents_names = lstrip_map.get(contig_name, None) or rstrip_map.get(contig_name, None)
                parents_names = (strip_parents_names and [strip_parents_names]) or regular_parents_names
                if parents_names:
                    for parent_name in parents_names:
                        copy_from_parent(contig, parent_name)
                else:
                    query_position_map[contig_name] = (contig.alignment.q_st, contig.alignment.q_ei)

        for child_name in children_names:
            set_query_position(child_name)

    for contig_name in sorted_roots:
        set_query_position(contig_name)

    def copy_takes_one_side(edge_table, overlap_xtake_map, overlap_xparent_map):
        for parent in edge_table:
            child_remainder = edge_table[parent]
            for child_remainder_morph in eqv_morphism_graph.get(child_remainder, [child_remainder]):
                if child_remainder_morph in overlap_xtake_map:
                    continue

                for parent_morph in eqv_morphism_graph.get(parent, [parent]):
                    for parent_remainder in overlap_xparent_map:
                        if overlap_xparent_map[parent_remainder] == parent_morph:
                            overlap_xtake_map[child_remainder_morph] = overlap_xtake_map[parent_remainder]
                            yield True

    # Closing `takes` by parents
    while list(copy_takes_one_side(combine_right_edge, overlap_lefttake_map, overlap_leftparent_map)): pass
    while list(copy_takes_one_side(combine_left_edge, overlap_righttake_map, overlap_rightparent_map)): pass

    final_nodes: List[str] = []
    final_parts: Dict[str, bool] = {}
    final_children_mapping: Dict[str, List[str]] = {}

    def add_join_parents(join_name):
        if join_name in children_join_points:
            for contig_name in parent_graph.get(join_name, [join_name]):
                add_join_parents(contig_name)
        else:
            final_nodes.append(join_name)

    for join_name in children_join_points + sorted_sinks:
        add_join_parents(join_name)

    def is_ancestor(contig_name, other_names):
        for other in other_names:
            if other == contig_name:
                continue

            if contig_name in transitive_children_graph.get(other, []):
                return True
        return False

    for contig_name in final_nodes[:]:
        if is_ancestor(contig_name, final_nodes):
            final_nodes.remove(contig_name)

    for contig_name in final_nodes:
        if any(contig_name in transitive_parent_graph.get(bad, []) for bad in bad_contigs):
            continue

        if any(contig_name in eqv_morphism_graph.get(temp_name, [temp_name]) for temp_name in overlaps_list):
            continue

        final_parts[contig_name] = True

    for contig_name in bad_contigs:
        final_parts[contig_name] = True

    for parent_name in sorted_roots:
        children = []
        for final_contig in final_parts:
            if final_contig == parent_name or \
               parent_name in reduced_parent_graph.get(final_contig, []):
                children.append(final_contig)

        final_children_mapping[parent_name] = children

    min_position, max_position = 1, 1
    position_offset = 100
    for _, contig in contig_map.items():
        if isinstance(contig, GenotypedContig) and contig.ref_seq is not None:
            max_position = max(max_position, len(contig.ref_seq) + 3 * position_offset)
        else:
            max_position = max(max_position, len(contig.seq) + 3 * position_offset)

    def overlaps(self, other) -> bool:
        def intervals_overlap(x, y):
            return x[0] <= y[1] and x[1] >= y[0]

        return intervals_overlap((self.alignment.q_st, self.alignment.q_ei),
                                 (other.alignment.q_st, other.alignment.q_ei))

    name_map = {}
    for i, (parent, children) in enumerate(sorted(final_children_mapping.items(), key=lambda p: p[0])):
        name_map[parent] = f"{i + 1}"

        unaligned_names = [name for name in children if name in unaligned]
        aligned_names = [name for name in children if name not in unaligned]

        todo_names = aligned_names
        for contig_name in unaligned_names:
            todo_names.append(contig_name)
            if contig_name not in discarded:
                discarded.append(contig_name)

        todo_names = list(sorted(todo_names, key=lambda name: query_position_map.get(name, (-1, -1))))
        for k, child_name in enumerate(todo_names):
            if len(todo_names) > 1:
                name_map[child_name] = f"{i + 1}.{k + 1}"
            else:
                name_map[child_name] = f"{i + 1}"

        for bad_name in bad_contigs:
            if bad_name not in children:
                if bad_name in transitive_parent_graph \
                   and parent in transitive_parent_graph[bad_name]:
                    k += 1
                    name_map[bad_name] = f"{i + 1}.{k + 1}"

    for contig_name, name in name_map.items():
        logger.debug(f"Contig name {contig_name!r} is displayed as {name!r}.")

    def get_neighbours(part, lookup):
        for clone in eqv_morphism_graph.get(part.name, [part.name]):
            maybe_name = lookup.get(clone, None)
            if maybe_name is not None:
                yield contig_map[maybe_name]

    def get_final_version(contig):
        name = reduced_morphism_graph.get(contig.name, [contig.name])[0] # FIXME: why 0???
        return contig_map[name]

    def get_neighbour(part, lookup):
        if not part: return None
        lst = list(get_neighbours(part, lookup))
        ret = max(map(get_final_version, lst), key=lambda contig: contig.alignment.ref_length, default=None)
        return ret

    aligned_size_map: Dict[str, Tuple[int, int]] = {}
    full_size_map: Dict[str, Tuple[int, int]] = {}

    for parent_name in sorted_roots:
        parts_names = final_children_mapping[parent_name]
        parts = [contig_map[part] for part in parts_names]

        for part in parts:
            if not isinstance(part, AlignedContig):
                continue

            prev_part = get_neighbour(part, overlap_righttake_map)
            next_part = get_neighbour(part, overlap_lefttake_map)

            if prev_part is not None:
                r_st = prev_part.alignment.r_st
            else:
                if part.name in bad_contigs:
                    start_delta = 0
                else:
                    start_delta = -1 * part.alignment.q_st
                r_st = part.alignment.r_st + start_delta

            if next_part is not None:
                r_ei = next_part.alignment.r_ei
            else:
                if part.name in bad_contigs:
                    end_delta = 0
                else:
                    end_delta = len(part.seq) - 1 - part.alignment.q_ei
                r_ei = part.alignment.r_ei + end_delta

            aligned_size_map[part.name] = (r_st, r_ei)

            sibling_name = ([overlap_sibling_map[name] for name in eqv_morphism_graph.get(part.name, [part.name]) if name in overlap_sibling_map] or [""])[0]
            sibling = sibling_name and contig_map[sibling_name]
            prev_part = get_neighbour(sibling, overlap_lefttake_map)
            next_part = get_neighbour(sibling, overlap_righttake_map)

            if prev_part is not None and prev_part.alignment.r_ei < part.alignment.r_st and prev_part:
                r_st = prev_part.alignment.r_st
            else:
                r_st = part.alignment.r_st

            if next_part is not None and next_part.alignment.r_st > part.alignment.r_ei and next_part:
                r_ei = next_part.alignment.r_ei
            else:
                r_ei = part.alignment.r_ei

            full_size_map[part.name] = (r_st, r_ei)

    def get_contig_coordinates(contig: GenotypedContig) -> Tuple[int, int, int, int]:
        if isinstance(contig, AlignedContig) and contig.alignment.ref_length > 0:
            r_st = contig.alignment.r_st
            r_ei = contig.alignment.r_ei
            if contig.name in aligned_size_map:
                a_r_st, a_r_ei = aligned_size_map[contig.name]
            else:
                a_r_st = r_st
                a_r_ei = r_ei
            if contig.name in full_size_map:
                f_r_st, f_r_ei = full_size_map[contig.name]
            else:
                f_r_st = r_st - contig.alignment.q_st
                f_r_ei = r_ei + (len(contig.seq) - contig.alignment.q_ei)
        else:
            f_r_st = 0
            f_r_ei = len(contig.seq)
            a_r_st = f_r_st
            a_r_ei = f_r_ei
        return (a_r_st, a_r_ei, f_r_st, f_r_ei)

    def get_tracks(repeatset: Set[str], group_ref: str, contig_name: str) -> Iterable[Track]:
        parts_names = final_children_mapping[contig_name]
        parts = [contig_map[name] for name in parts_names]
        parts = list(sorted(parts, key=lambda part: part.alignment.r_st if isinstance(part, AlignedContig) else -1))
        for prev_part, part, next_part in sliding_window(parts):
            if part.name in repeatset:
                continue

            if part.name in bad_contigs:
                continue

            if not isinstance(part, AlignedContig):
                continue

            if part.group_ref != group_ref:
                continue

            repeatset.add(part.name)
            indexes = name_map[part.name]
            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(part)

            if a_r_st < f_r_st:
                yield Track(min(a_r_st, f_r_st) + position_offset, max(a_r_st, f_r_st) + position_offset, color="yellow")

            if a_r_ei > f_r_ei:
                yield Track(min(a_r_ei, f_r_ei) + position_offset, max(a_r_ei, f_r_ei) + position_offset, color="yellow")

            yield Track(f_r_st + position_offset, f_r_ei + position_offset, label=f"{indexes}")

    def get_arrows(repeatset: Set[str], group_ref: str, contig_name: str, labels: bool) -> Iterable[Arrow]:
        parts = final_children_mapping[contig_name]
        for part_name in parts:
            part = contig_map[part_name]

            if part.name in repeatset:
                continue

            if part.name in bad_contigs:
                continue

            if not isinstance(part, AlignedContig):
                continue

            if part.group_ref != group_ref:
                continue

            repeatset.add(part.name)
            indexes = name_map[part.name] if labels else None
            height = 20 if labels else 1
            elevation = 1 if labels else -20
            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(part)
            yield Arrow(a_r_st + position_offset, a_r_ei + position_offset,
                        elevation=elevation,
                        h=height,
                        label=indexes)

    def get_all_arrows(group_ref: str, labels: bool) -> Iterable[Arrow]:
        repeatset: Set[str] = set()
        for parent_name in sorted_roots:
            yield from get_arrows(repeatset, group_ref, parent_name, labels)

    ################
    # Drawing part #
    ################

    landmarks_path = (Path(__file__).parent.parent / "data" /
                      "landmark_references.yaml")
    landmark_groups = yaml.safe_load(landmarks_path.read_text())
    projects = ProjectConfig.loadDefault()
    figure = Figure()
    for group_ref in group_refs:
        matching_groups = [group for group in landmark_groups if group['coordinates'] == group_ref]
        if matching_groups:
            reference_set = matching_groups[0]
        elif "HIV1" in group_ref:
            matching_groups = [group for group in landmark_groups if group['coordinates'] == "HIV1-B-FR-K03455-seed"]
            reference_set = matching_groups[0]
        else:
            reference_set = None

        #############
        # Landmarks #
        #############

        if reference_set:

            # Filling out missing ends.
            prev_landmark = None
            for landmark in sorted(reference_set['landmarks'], key=itemgetter('start')):
                landmark.setdefault('frame', 0)
                if prev_landmark and 'end' not in prev_landmark:
                    prev_landmark['end'] = landmark['start'] - 1
                prev_landmark = landmark

            # Computing the stretching factor.
            landmark_max = 0
            for landmark in reference_set['landmarks']:
                landmark_max = max(landmark_max, landmark['end'])

            stretch_c = group_refs[group_ref] / landmark_max

            # Drawing the landmarks.
            for frame, frame_landmarks in groupby(reference_set['landmarks'],
                                                  itemgetter('frame')):
                subtracks = []
                for landmark in frame_landmarks:
                    landmark_colour = landmark.get('colour')
                    if landmark_colour is None:
                        continue
                    subtracks.append(Track(landmark['start'] * stretch_c + position_offset,
                                           landmark['end'] * stretch_c + position_offset,
                                           label=landmark['name'],
                                           color=landmark_colour))
                figure.add(Multitrack(subtracks))

        # Drawing the reference sequence.
        r_st = 0
        r_ei = group_refs[group_ref]
        figure.add(Track(r_st + position_offset, r_ei + position_offset, label=f"{group_ref}"))

        ##########
        # Arrows #
        ##########

        ref_arrows = list(get_all_arrows(group_ref, labels=True))
        if ref_arrows:
            figure.add(ArrowGroup(ref_arrows))

        ###########
        # Contigs #
        ###########

        repeatset1: Set[str] = set()
        repeatset2: Set[str] = set()
        for parent_name in sorted_roots:
            arrows = list(get_arrows(repeatset1, group_ref, parent_name, labels=False))
            if arrows:
                figure.add(ArrowGroup(arrows))
            parts = list(get_tracks(repeatset2, group_ref, parent_name))
            if parts:
                figure.add(Multitrack(parts))

        #############
        # Discarded #
        #############

        if discarded:
            label = LeftLabel(text=f"discards:", x=0, font_size=12)
            pos = position_offset / 2
            figure.add(Track(pos, pos, h=40, label=label))
            for parent_name in sorted_roots:
                contigs = final_children_mapping.get(parent_name, [])
                for contig_name in contigs:
                    if contig_name not in discarded:
                        continue

                    contig = contig_map[contig_name]
                    (r_st, r_ei, f_r_st, f_r_ei) = get_contig_coordinates(contig)
                    name = name_map.get(contig_name, contig_name)
                    if isinstance(contig, AlignedContig) and contig.name not in unaligned:
                        colour = 'lightgrey'
                        figure.add(Arrow(r_st + position_offset, r_ei + position_offset, elevation=-20, h=1))
                    else:
                        colour = "yellow"
                    figure.add(Track(f_r_st + position_offset, f_r_ei + position_offset, label=name, color=colour))

        #############
        # Anomalies #
        #############

        if anomaly:
            label = LeftLabel(text=f"anomaly:", x=0, font_size=12)
            pos = position_offset / 2
            figure.add(Track(pos, pos, h=40, label=label))
            for parent_name in sorted_roots:
                contigs = final_children_mapping.get(parent_name, [])
                for contig_name in contigs:
                    if contig_name not in anomaly:
                        continue

                    contig = contig_map[contig_name]
                    (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(contig)
                    if isinstance(contig, AlignedContig):
                        colour = "lightgray"
                        if contig.strand == "reverse":
                            figure.add(Arrow(a_r_ei + position_offset, a_r_st + position_offset, elevation=-20, h=1))
                        else:
                            figure.add(Arrow(a_r_st + position_offset, a_r_ei + position_offset, elevation=-20, h=1))
                    else:
                        colour = "yellow"

                    name = name_map.get(contig_name, contig_name)
                    figure.add(Track(a_r_st + position_offset, a_r_ei + position_offset, color=colour, label=name))

    ###########
    # Unknown #
    ###########

    if unknown:
        label = LeftLabel(text=f"unknown:", x=0, font_size=12)
        pos = position_offset / 2
        figure.add(Track(pos, pos, h=40, label=label))
        for parent_name in sorted_roots:
            contigs = final_children_mapping.get(parent_name, [])
            for contig_name in contigs:
                if contig_name not in unknown:
                    continue

                contig = contig_map[contig_name]
                r_st = 0
                r_ei = len(contig.seq)
                colour = "yellow"
                name = name_map.get(contig_name, contig_name)
                figure.add(Track(r_st + position_offset, r_ei + position_offset, color=colour, label=name))

    if not figure.elements:
        figure.add(Track(1, max_position, label='No contigs found.', color='none'))
    return figure


def map_references(contig_ref_name: str,
                   coordinates_name: str,
                   projects: ProjectConfig) -> typing.Mapping[int, int]:
    ref_seq = projects.getReference(contig_ref_name)
    try:
        coordinates_seq = projects.getGenotypeReference(coordinates_name)
    except KeyError:
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
                 blast_rows,
                 use_concordance=False):
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
            field_names = (pos_field, 'coverage', 'dels')
            for field_name in field_names:
                field_text = contig_row[field_name]
                field_value = None if field_text == '' else int(field_text)
                contig_row[field_name] = field_value
            if use_concordance:
                field_text = contig_row['concordance']
                field_value = None if field_text == '' else 100 * float(field_text)
                contig_row['concordance'] = field_value
        start = contig_rows[0][pos_field]
        end = contig_rows[-1][pos_field]
        new_final_pos = -1
        while end is None:
            # this can happen if the match ends with an insertion. Backtrack to the last position that was not None
            new_final_pos -= 1
            end = contig_rows[new_final_pos][pos_field]
        coverage = [0] * (end - start + 1)
        concordance = [0] * (end - start + 1)
        pos = 0
        for contig_row in contig_rows:
            link = contig_row.get('link')
            pos = contig_row[pos_field]
            if pos is None:
                insertion_size += 1
            else:
                if insertion_size:
                    insertion_ranges.append((pos, pos+insertion_size-1))
                    insertion_size = 0
                if contig_row['coverage'] is not None:
                    coverage[pos - start] = (contig_row['coverage']) - contig_row['dels']
                if use_concordance and contig_row['concordance'] is not None:
                    concordance[pos - start] = contig_row['concordance']
                if link != 'D':
                    contig_pos = int(contig_row['query_nuc_pos'])
                    while event_positions and event_positions[-1] <= contig_pos:
                        event_pos = event_positions.pop()
                        for blast_num in blast_starts[event_pos]:
                            blast_ranges[blast_num-1][0] = pos
                        for blast_num in blast_ends[event_pos]:
                            blast_ranges[blast_num-1][1] = pos
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
                blast_ranges[blast_num - 1][0] = end
            for blast_num in blast_ends[event_pos]:
                blast_ranges[blast_num - 1][1] = end

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
            if not use_concordance:
                f.add(ShadedCoverage(start + position_offset,
                                     end + position_offset,
                                     coverage),
                      gap=-4)
            elif max(concordance) > 0:
                f.add(ConcordanceLine(start + position_offset,
                                      end + position_offset,
                                      concordance),
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


def summarize_figure(figure: Figure, is_concordance=False):
    """ Summarize the contents of a figure to text.

    Useful for testing.
    """
    figure.show()  # Test that all the display math works.
    if is_concordance:
        coverage_or_concordance = 'Concordance '
    else:
        coverage_or_concordance = 'Coverage '

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
                summary.write(coverage_or_concordance)
                summary.write(', '.join(map(str, ys)))
                continue
            coverage_groups = getattr(span, 'coverage_groups', None)
            if coverage_groups is not None:
                summary.write(coverage_or_concordance)
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
    parser.add_argument('--concordance',
                        help='Make concordance plot instead of coverage',
                        action='store_true')
    args = parser.parse_args()

    plot_genome_coverage(args.genome_coverage_csv,
                         args.blast,
                         args.genome_coverage_svg,
                         args.concordance)
    print('Wrote', args.genome_coverage_svg)


if __name__ == '__main__':
    main()
