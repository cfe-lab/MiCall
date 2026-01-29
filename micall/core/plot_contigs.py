import typing
from typing import Dict, Tuple, List, Set, Iterable, NoReturn, Sequence, Union
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictReader
from io import StringIO
from itertools import groupby
from math import log10, copysign, floor
from operator import itemgetter, attrgetter
from pathlib import Path
import logging

import yaml
from aligntools import CigarHit
from genetracks import Figure, Track, Multitrack, Coverage
# noinspection PyPep8Naming
import drawsvg as draw
from genetracks.elements import Element, Label
from matplotlib import cm, colors
from matplotlib.colors import Normalize

from micall.core.project_config import ProjectConfig
from micall.utils.alignment_wrapper import align_nucs
from micall.utils.contig_stitcher_contigs import Contig, GenotypedContig, AlignedContig
from micall.utils.contig_stitcher_context import ReferencefullStitcherContext
import micall.utils.referencefull_contig_stitcher_events as events
from micall.data.landmark_reader import LandmarkReader


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


def plot_stitcher_coverage(logs: Iterable[events.EventType], genome_coverage: Path):
    with ReferencefullStitcherContext.stage():
        f = build_stitcher_figure(logs)
        f.show(w=970).save_svg(genome_coverage, context=draw.Context(invert_y=True))
        return f


def build_stitcher_figure(logs: Iterable[events.EventType]) -> Figure:
    complete_contig_map: Dict[int, GenotypedContig] = {}
    name_map: Dict[int, str] = {}
    complete_parent_graph: Dict[int, List[int]] = {}
    alive_set: Set[int] = set()
    morphism_graph: Dict[int, List[int]] = {}
    reduced_parent_graph: Dict[int, List[int]] = {}
    transitive_parent_graph: Dict[int, List[int]] = {}
    discarded: List[int] = []
    unknown: List[int] = []
    anomaly: List[int] = []
    anomaly_data_map: Dict[int, Union[events.ZeroHits, events.StrandConflict]] = {}
    unaligned_map: Dict[int, List[CigarHit]] = {}
    overlaps_list: List[int] = []
    overlap_leftparent_map: Dict[int, int] = {}
    overlap_rightparent_map: Dict[int, int] = {}
    overlap_lefttake_map: Dict[int, int] = {}
    overlap_righttake_map: Dict[int, int] = {}
    overlap_left_sibling: Dict[int, int] = {}
    overlap_right_sibling: Dict[int, int] = {}
    combine_left_edge: Dict[int, int] = {}
    combine_right_edge: Dict[int, int] = {}
    children_join_points: List[int] = []
    query_position_map: Dict[int, Tuple[int, int]] = {}
    lstrip_map: Dict[int, int] = {}
    rstrip_map: Dict[int, int] = {}

    def remove_intermediate_edges(graph):
        tr_cl = transitive_closure(graph)
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                if all(other not in tr_cl.get(child, []) for other in children):
                    lst.append(child)
            ret[parent] = lst
        return ret

    def remove_transitive_edges(graph):
        tr_cl = transitive_closure(graph)
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                is_transitive = any(child in tr_cl.get(other_node, [])
                                    for other_node in children
                                    if other_node != child)
                if not is_transitive:
                    lst.append(child)
            ret[parent] = lst
        return ret

    def remove_duplicate_edges(graph):
        ret = {}
        for parent, children in graph.items():
            lst = []
            for child in children:
                if child not in lst:
                    lst.append(child)
            ret[parent] = lst
        return ret

    def get_transitive_children(recur, lst, graph, current):
        for child in graph.get(current, []):
            if child not in recur:
                recur.add(child)
                lst.append(child)
                get_transitive_children(recur, lst, graph, child)

    def transitive_closure(graph):
        ret = {}
        for parent in graph:
            children = []
            get_transitive_children(set(), children, graph, parent)
            ret[parent] = children
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

    def record_contig(contig: GenotypedContig, parents: Iterable[GenotypedContig]):
        complete_contig_map[contig.id] = contig
        if [contig.id] != [parent.id for parent in parents]:
            for parent in parents:
                complete_contig_map[parent.id] = parent
                if contig.id not in complete_parent_graph:
                    complete_parent_graph[contig.id] = []

                complete_parent_graph[contig.id].append(parent.id)

    def record_alive(contig: Contig):
        alive_set.add(contig.id)

    def record_bad_contig(contig: GenotypedContig, lst: List[int]):
        complete_contig_map[contig.id] = contig
        if contig.id not in lst:
            lst.append(contig.id)

    def record_lstrip(result: AlignedContig, original: AlignedContig):
        lstrip_map[result.id] = original.id

    def record_rstrip(result: AlignedContig, original: AlignedContig):
        rstrip_map[result.id] = original.id

    def hit_to_insertions(contig: GenotypedContig, hit: CigarHit):
        yield CigarHit.from_default_alignment(q_st=0, q_ei=hit.q_st - 1, r_st=hit.r_st, r_ei=hit.r_st - 1)
        yield from hit.insertions()
        yield CigarHit.from_default_alignment(q_st=hit.q_ei + 1, q_ei=len(contig.seq) - 1,
                                              r_st=hit.r_ei + 1, r_ei=hit.r_ei)

    def hits_to_insertions(contig: GenotypedContig, hits: Iterable[CigarHit]):
        for hit in hits:
            yield from hit_to_insertions(contig, hit)

    def record_initial_hit(contig: GenotypedContig, hits: Sequence[CigarHit]):
        insertions = [gap for gap in hits_to_insertions(contig, hits)]
        unaligned_map[contig.id] = insertions

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
        elif isinstance(event, events.ConnectedHit):
            record_contig(event.part, [event.contig])
            record_alive(event.part)
        elif isinstance(event, events.NoRef):
            record_bad_contig(event.contig, unknown)
            record_alive(event.contig)
        elif isinstance(event, events.ZeroHits):
            record_bad_contig(event.contig, anomaly)
            anomaly_data_map[event.contig.id] = event
            record_alive(event.contig)
        elif isinstance(event, events.StrandConflict):
            record_bad_contig(event.contig, anomaly)
            anomaly_data_map[event.contig.id] = event
            record_alive(event.contig)
        elif isinstance(event, events.ReverseComplement):
            record_contig(event.result, [event.contig])
            record_alive(event.result)
        elif isinstance(event, events.HitNumber):
            record_initial_hit(event.contig, event.connected)
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
            pass
        elif isinstance(event, events.Overlap):
            overlaps_list.append(event.left_overlap.id)
            overlaps_list.append(event.right_overlap.id)
            overlap_leftparent_map[event.left_remainder.id] = event.left.id
            overlap_rightparent_map[event.right_remainder.id] = event.right.id
            overlap_lefttake_map[event.left_remainder.id] = event.left_take.id
            overlap_righttake_map[event.right_remainder.id] = event.right_take.id
            overlap_left_sibling[event.left_remainder.id] = event.right_remainder.id
            overlap_right_sibling[event.right_remainder.id] = event.left_remainder.id
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
        elif isinstance(event, events.NoOverlap):
            record_alive(event.contig)
        elif isinstance(event, events.Cut):
            record_contig(event.left, [event.original])
            record_contig(event.right, [event.original])
        elif isinstance(event, events.Combine):
            record_alive(event.result)
            record_contig(event.result, event.contigs)
            if event.contigs:
                combine_left_edge[event.result.id] = event.contigs[0].id
                combine_right_edge[event.result.id] = event.contigs[-1].id
        elif isinstance(event, (events.IgnoreGap, events.InitialHit, events.IgnoreCoverage)):
            pass
        else:
            _x: NoReturn = event
            raise RuntimeError(f"Unrecognized action or event: {event}")

    notransitive_parent_graph = remove_transitive_edges(complete_parent_graph)
    nodup_parent_graph = remove_duplicate_edges(notransitive_parent_graph)

    # Close alive set by parents
    def extend_alive(contig_id):
        alive_set.add(contig_id)
        for parent_id in nodup_parent_graph.get(contig_id, []):
            extend_alive(parent_id)

    for contig_id in alive_set.copy():
        extend_alive(contig_id)

    parent_graph: Dict[int, List[int]] = {}
    for contig_id in nodup_parent_graph:
        if contig_id in alive_set:
            parent_graph[contig_id] = nodup_parent_graph[contig_id]

    contig_map: Dict[int, GenotypedContig] = {k: v for k, v in complete_contig_map.items() if k in alive_set}
    bad_contigs = anomaly + discarded + unknown
    group_refs = {contig.group_ref: len(contig.ref_seq) for contig in contig_map.values() if contig.ref_seq}
    children_graph = inverse_graph(parent_graph)
    transitive_parent_graph = transitive_closure(parent_graph)
    transitive_children_graph = transitive_closure(children_graph)
    reduced_parent_graph = remove_intermediate_edges(transitive_parent_graph)
    sorted_roots = list(sorted(parent_id for
                               parent_id in contig_map
                               if parent_id not in parent_graph))
    sorted_sinks = list(sorted(child_id for
                               child_id in contig_map
                               if child_id not in children_graph))

    lstrip_set = set(lstrip_map.keys())
    rstrip_set = set(rstrip_map.keys())

    for contig_id, parents in parent_graph.items():
        if len(parents) == 1:
            morphism_graph[parents[0]] = [contig_id]

    transitive_morphism_graph = transitive_closure(morphism_graph)
    reduced_morphism_graph = remove_intermediate_edges(transitive_morphism_graph)
    eqv_morphism_graph = reflexive_closure(symmetric_closure(transitive_morphism_graph))

    for contig_id, parents in parent_graph.items():
        if len(parents) > 1:
            children_join_points.append(contig_id)

    def set_query_position(contig_id: int) -> None:
        contig = contig_map[contig_id]
        children_ids = children_graph.get(contig.id, [])

        def copy_from_parent(contig: AlignedContig, parent_id: int) -> None:
            if parent_id in query_position_map:
                (original_q_st, original_q_ei) = query_position_map[parent_id]
                (current_q_st, current_q_ei) = (contig.alignment.q_st, contig.alignment.q_ei)
                current_query_len = abs(current_q_st - current_q_ei)

                if contig_id in lstrip_map:
                    query_position_map[contig.id] = (original_q_ei - current_query_len, original_q_ei)
                elif contig_id in rstrip_map:
                    query_position_map[contig.id] = (original_q_st, original_q_st + current_query_len)
                else:
                    query_position_map[contig_id] = query_position_map[parent_id]

        if contig_id not in query_position_map:
            if isinstance(contig, AlignedContig):
                regular_parents_ids = parent_graph.get(contig_id, [])
                regular_parents_ids = [name for name in regular_parents_ids if name in query_position_map]
                strip_parents_ids = lstrip_map.get(contig_id, None) or rstrip_map.get(contig_id, None)
                parents_ids = (strip_parents_ids and [strip_parents_ids]) or regular_parents_ids
                if parents_ids:
                    for parent_id in parents_ids:
                        copy_from_parent(contig, parent_id)
                else:
                    query_position_map[contig_id] = (contig.alignment.q_st, contig.alignment.q_ei)

        for child_id in children_ids:
            set_query_position(child_id)

    for contig_id in sorted_roots:
        set_query_position(contig_id)

    def copy_takes_one_side(edge_table, overlap_xtake_map, overlap_xparent_map, overlap_xsibling, xstrip_set):
        for parent in edge_table:
            child_remainder = edge_table[parent]
            for child_remainder_morph in eqv_morphism_graph.get(child_remainder, [child_remainder]):
                for parent_morph in eqv_morphism_graph.get(parent, [parent]):
                    if child_remainder_morph in xstrip_set:
                        xstrip_set.add(parent_morph)
                    if parent_morph in xstrip_set:
                        xstrip_set.add(child_remainder_morph)

                    if child_remainder_morph in overlap_xtake_map:
                        continue
                    for parent_remainder in overlap_xparent_map:
                        if overlap_xparent_map[parent_remainder] == parent_morph:
                            overlap_xtake_map[child_remainder_morph] = overlap_xtake_map[parent_remainder]
                            overlap_xsibling[child_remainder_morph] = overlap_xsibling[parent_remainder]
                            yield True

    # Closing `takes` by parents
    while list(copy_takes_one_side(combine_right_edge, overlap_lefttake_map,
                                   overlap_leftparent_map, overlap_left_sibling, rstrip_set)):
        pass
    while list(copy_takes_one_side(combine_left_edge, overlap_righttake_map,
                                   overlap_rightparent_map, overlap_right_sibling, lstrip_set)):
        pass

    final_nodes: List[int] = []
    final_parts: Dict[int, bool] = {}
    final_children_mapping: Dict[int, List[int]] = {}

    def add_join_parents(join_id):
        if join_id in children_join_points:
            for contig_id in parent_graph.get(join_id, [join_id]):
                add_join_parents(contig_id)
        else:
            final_nodes.append(join_id)

    for join_id in children_join_points + sorted_sinks:
        add_join_parents(join_id)

    def is_ancestor(contig_id, other_ids):
        for other in other_ids:
            if other == contig_id:
                continue

            if contig_id in transitive_children_graph.get(other, []):
                return True
        return False

    for contig_id in final_nodes[:]:
        if is_ancestor(contig_id, final_nodes):
            final_nodes.remove(contig_id)

    for contig_id in final_nodes:
        if any(contig_id in eqv_morphism_graph.get(bad, []) for bad in bad_contigs):
            continue

        if any(contig_id in eqv_morphism_graph.get(temp_id, [temp_id]) for temp_id in overlaps_list):
            continue

        final_parts[contig_id] = True

    for contig_id in bad_contigs:
        final_parts[contig_id] = True

    for parent_id in sorted_roots:
        children = []
        for final_contig in final_parts:
            if final_contig == parent_id or \
               parent_id in reduced_parent_graph.get(final_contig, [final_contig]):
                children.append(final_contig)

        final_children_mapping[parent_id] = children

    aligned_size_map: Dict[int, Tuple[int, int]] = {}
    full_size_map: Dict[int, Tuple[int, int]] = {}

    def get_neighbours(part, lookup):
        for clone in eqv_morphism_graph.get(part.id, [part.id]):
            maybe_id = lookup.get(clone, None)
            if maybe_id is not None:
                yield contig_map[maybe_id]

    def get_final_version(contig):
        [name] = reduced_morphism_graph.get(contig.id, [contig.id])
        return contig_map[name]

    def get_neighbour(part, lookup):
        if not part:
            return None
        lst = list(get_neighbours(part, lookup))
        ret = max(map(get_final_version, lst), key=lambda contig: contig.alignment.ref_length, default=None)
        return ret

    def get_contig_coordinates(contig: GenotypedContig) -> Tuple[int, int, int, int]:
        if isinstance(contig, AlignedContig) and contig.alignment.ref_length > 0:
            r_st = contig.alignment.r_st
            r_ei = contig.alignment.r_ei
            if contig.id in aligned_size_map:
                a_r_st, a_r_ei = aligned_size_map[contig.id]
            else:
                a_r_st = r_st
                a_r_ei = r_ei
            if contig.id in full_size_map:
                f_r_st, f_r_ei = full_size_map[contig.id]
            else:
                f_r_st = r_st - contig.alignment.q_st
                f_r_ei = r_ei + (len(contig.seq) - contig.alignment.q_ei)
        else:
            f_r_st = 0
            f_r_ei = len(contig.seq)
            a_r_st = f_r_st
            a_r_ei = f_r_ei
        return (a_r_st, a_r_ei, f_r_st, f_r_ei)

    for parent_id in sorted_roots:
        parts_ids = final_children_mapping[parent_id]
        for part_id in parts_ids:
            part = contig_map[part_id]
            if not isinstance(part, AlignedContig):
                continue

            prev_part = get_neighbour(part, overlap_righttake_map)
            next_part = get_neighbour(part, overlap_lefttake_map)

            if prev_part is not None:
                r_st = prev_part.alignment.r_st
            elif part_id in lstrip_set:
                r_st = part.alignment.r_st
            else:
                start_delta = -1 * part.alignment.q_st
                r_st = part.alignment.r_st + start_delta

            if next_part is not None:
                r_ei = next_part.alignment.r_ei
            elif part_id in rstrip_set:
                r_ei = part.alignment.r_ei
            else:
                end_delta = len(part.seq) - 1 - part.alignment.q_ei
                r_ei = part.alignment.r_ei + end_delta

            aligned_size_map[part.id] = (r_st, r_ei)

            sibling_left_id = ([overlap_left_sibling[name]
                                for name in eqv_morphism_graph.get(part.id, [part.id])
                                if name in overlap_left_sibling] or [0])[0]
            sibling_left = sibling_left_id and contig_map[sibling_left_id]
            sibling_right_id = ([overlap_right_sibling[name]
                                 for name in eqv_morphism_graph.get(part.id, [part.id])
                                 if name in overlap_right_sibling] or [0])[0]
            sibling_right = sibling_right_id and contig_map[sibling_right_id]
            prev_part = get_neighbour(sibling_right, overlap_lefttake_map)
            next_part = get_neighbour(sibling_left, overlap_righttake_map)

            if prev_part is not None:
                r_st = prev_part.alignment.r_st
            else:
                r_st = part.alignment.r_st

            if next_part is not None:
                r_ei = next_part.alignment.r_ei
            else:
                r_ei = part.alignment.r_ei

            full_size_map[part.id] = (r_st, r_ei)

    def carve_gap(gap: CigarHit, aligned_parts: Iterable[AlignedContig]):
        for contig in aligned_parts:
            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(contig)
            other_coords = query_position_map.get(contig.id, (-1, -2))

            other_q_st = min(other_coords) - max(0, abs(f_r_st - a_r_st))
            other_q_ei = max(other_coords) + max(0, abs(a_r_ei - f_r_ei))

            if gap.q_st <= other_q_st and gap.q_ei >= other_q_st:
                q_st = gap.q_st
                q_ei = other_q_st - 1
            elif gap.q_ei >= other_q_ei and gap.q_ei <= other_q_ei:
                q_st = other_q_ei + 1
                q_ei = gap.q_ei
            elif gap.q_st >= other_q_st and gap.q_ei <= other_q_ei:
                return None
            else:
                continue

            if q_st >= other_q_st and q_ei <= other_q_ei:
                return None

            if q_st > q_ei:
                return None

            gap = CigarHit.from_default_alignment(q_st=q_st, q_ei=q_ei, r_st=gap.r_st, r_ei=gap.r_ei)

        if gap.query_length > 0:
            return gap

    def collect_gaps(root: int, children_ids: List[int]):
        all_children = [contig_map[name] for name in children_ids]
        children = [child for child in all_children if isinstance(child, AlignedContig)]
        for name in unaligned_map:
            if reduced_parent_graph.get(name, [name]) == [root]:
                for gap in unaligned_map[name]:
                    carved = carve_gap(gap, children)
                    if carved is not None:
                        yield carved

    carved_unaligned_parts: Dict[int, List[int]] = {}
    for root in sorted_roots:
        existing: Set[Tuple[int, int]] = set()
        children = final_children_mapping[root]
        for gap in collect_gaps(root, children):
            coords = (gap.q_st, gap.q_ei)
            if coords not in existing:
                existing.add(coords)
                if root not in carved_unaligned_parts:
                    carved_unaligned_parts[root] = []
                fake = Contig(name=None, seq="", reads_count=None)
                carved_unaligned_parts[root].append(fake.id)
                query_position_map[fake.id] = coords

    merged_unaligned_parts: Dict[int, List[int]] = {}
    for root in sorted_roots:
        children = final_children_mapping[root]
        unaligned_children = carved_unaligned_parts.get(root, [])
        todo = children + unaligned_children
        todo = list(sorted(todo, key=lambda name: query_position_map.get(name, (-1, -1))))
        current_group = []
        for child_id in todo + [None]:
            if child_id in unaligned_children:
                coords = query_position_map[child_id]
                current_group.append(coords)
            elif current_group:
                coords = (min(q_st for q_st, q_ei in current_group),
                          max(q_ei for q_st, q_ei in current_group))
                if root not in merged_unaligned_parts:
                    merged_unaligned_parts[root] = []
                fake = Contig(name=None, seq="", reads_count=None)
                query_position_map[fake.id] = coords
                merged_unaligned_parts[root].append(fake.id)
                current_group = []

    name_map = {}
    for i, root in enumerate(sorted_roots):
        children = final_children_mapping[root]
        unaligned_children = merged_unaligned_parts.get(root, [])

        name_map[root] = f"{i + 1}"

        todo_ids = children + unaligned_children
        todo_ids = list(sorted(todo_ids, key=lambda name: query_position_map.get(name, (-1, -1))))
        for k, child_id in enumerate(todo_ids):
            if len(todo_ids) > 1:
                name_map[child_id] = f"{i + 1}.{k + 1}"
            else:
                name_map[child_id] = f"{i + 1}"

        for bad_id in bad_contigs:
            if bad_id not in children:
                if bad_id in transitive_parent_graph \
                   and root in transitive_parent_graph[bad_id]:
                    k += 1
                    name_map[bad_id] = f"{i + 1}.{k + 1}"

    for contig_id, name in name_map.items():
        if contig_id in complete_contig_map:
            contig = complete_contig_map[contig_id]
            logger.debug(f"Contig name {contig.unique_name} is displayed as {name!r}.")

    def get_tracks(parts: Iterable[GenotypedContig]) -> Iterable[Track]:
        for part in parts:
            name = name_map[part.id]
            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(part)

            if a_r_st < f_r_st:
                yield Track(min(a_r_st, f_r_st) + position_offset,
                            max(a_r_st, f_r_st) + position_offset, color="yellow")

            if a_r_ei > f_r_ei:
                yield Track(min(a_r_ei, f_r_ei) + position_offset,
                            max(a_r_ei, f_r_ei) + position_offset, color="yellow")

            if isinstance(part, AlignedContig):
                colour = 'lightgrey'
            else:
                colour = "yellow"

            yield Track(f_r_st + position_offset, f_r_ei + position_offset, label=f"{name}", color=colour)

    def get_anomaly_tracks(part: GenotypedContig, hits: Sequence[CigarHit], name: str) -> Iterable[Track]:
        """Generate tracks for anomaly contigs with alignment hits, showing mapped (grey) and unmapped (yellow) regions."""
        # Calculate the span of all hits
        min_r_st = min(hit.r_st for hit in hits)
        max_r_ei = max(hit.r_ei for hit in hits)

        # Calculate the aligned region on reference (all hits)
        # And the full region including unmapped parts
        min_q_st = min(hit.q_st for hit in hits)
        max_q_ei = max(hit.q_ei for hit in hits)

        # Full contig region (including unmapped parts at start and end)
        f_r_st = min_r_st - min_q_st
        f_r_ei = max_r_ei + (len(part.seq) - 1 - max_q_ei)

        # Aligned region
        a_r_st = min_r_st
        a_r_ei = max_r_ei

        # Yellow track for unmapped region at start
        if a_r_st > f_r_st:
            yield Track(f_r_st + position_offset, a_r_st + position_offset, color="yellow")

        # Yellow track for unmapped region at end
        if f_r_ei > a_r_ei:
            yield Track(a_r_ei + position_offset, f_r_ei + position_offset, color="yellow")

        # Grey track for the aligned/mapped region
        yield Track(a_r_st + position_offset, a_r_ei + position_offset, label=f"{name}", color="lightgrey")

    def get_anomaly_arrows(hits: Sequence[CigarHit], strands: Sequence[typing.Literal["forward", "reverse"]],
                          name: str) -> Iterable[Arrow]:
        """Generate arrows for anomaly contigs showing alignment hits with strand information."""
        arrows = []

        # Group hits by strand to assign elevations
        forward_indices = [i for i, strand in enumerate(strands) if strand == "forward"]
        reverse_indices = [i for i, strand in enumerate(strands) if strand == "reverse"]

        # Create arrows for forward strand hits (elevation = 0)
        for idx_pos, i in enumerate(forward_indices):
            hit = hits[i]
            r_st = hit.r_st + position_offset
            r_ei = hit.r_ei + position_offset
            arrows.append(Arrow(r_st, r_ei, h=3, elevation=-4, label=None))

        # Create arrows for reverse strand hits (elevation = 1 for separation)
        for idx_pos, i in enumerate(reverse_indices):
            hit = hits[i]
            r_st = hit.r_st + position_offset
            r_ei = hit.r_ei + position_offset
            # Reverse the arrow direction by swapping start/end
            arrows.append(Arrow(r_ei, r_st, h=3, elevation=-4, label=None))

        return arrows

    def get_arrows(parts: Iterable[GenotypedContig], labels: bool) -> Iterable[Arrow]:
        for part in parts:
            name = name_map[part.id] if labels else None
            height = 20 if labels else 1
            elevation = 1 if labels else -20
            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(part)

            if isinstance(part, AlignedContig) and part.strand == "reverse":
                tmp = a_r_st
                a_r_st = a_r_ei
                a_r_ei = tmp

            yield Arrow(a_r_st + position_offset, a_r_ei + position_offset,
                        elevation=elevation,
                        h=height,
                        label=name)

    def make_ray() -> Element:
        screen_size = (max_position - min_position) + position_offset / 2
        single_size = 0.02 * screen_size

        def generate_beams():
            for i in range(floor(screen_size / single_size) + 1):
                if i % 2 == 0:
                    yield Track(i * single_size + min_position + position_offset / 2,
                                (i + 1) * single_size + min_position + position_offset / 2,
                                h=0.1, color="green")

        return Multitrack(list(generate_beams()))

    def add_section(title: str) -> None:
        label = LeftLabel(text=title, x=0, font_size=12)
        pos = position_offset / 2
        figure.add(Arrow(pos, pos, h=0))
        figure.add(make_ray())
        figure.add(Arrow(pos, pos, h=0))
        figure.add(Track(pos, pos, label=label, h=0))

    min_position = 0
    max_position = max(group_refs.values(), default=1)
    for contig_id in final_parts:
        contig = contig_map[contig_id]
        if isinstance(contig, AlignedContig):
            positions = get_contig_coordinates(contig)
            max_position = max(max_position, max(positions))
            min_position = min(min_position, min(positions))
        else:
            max_position = max(max_position, len(contig.seq))

    position_offset = -1 * min_position + 0.05 * (max_position - min_position)

    ################
    # Drawing part #
    ################

    landmark_reader = LandmarkReader.load()
    figure = Figure()
    for group_ref in group_refs:
        try:
            if group_ref is not None:
                landmarks = landmark_reader.get_landmarks(group_ref)
            else:
                landmarks = None
        except ValueError:
            landmarks = None

        #############
        # Landmarks #
        #############

        if landmarks:
            # Filling out missing ends.
            prev_landmark = None
            for landmark in sorted(landmarks, key=itemgetter('start')):
                landmark.setdefault('frame', 0)
                if prev_landmark and 'end' not in prev_landmark:
                    prev_landmark['end'] = landmark['start'] - 1
                prev_landmark = landmark

            # Computing the stretching factor.
            landmark_max = 0
            for landmark in landmarks:
                landmark_max = max(landmark_max, landmark['end'])

            stretch_c = group_refs[group_ref] / landmark_max

            # Drawing the landmarks.
            for frame, frame_landmarks in groupby(landmarks, itemgetter('frame')):
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

        #############
        # Reference #
        #############

        r_st = 0
        r_ei = group_refs[group_ref]
        reference_tracks = []
        reference_min = r_st + position_offset
        reference_max = r_ei + position_offset
        reference_tracks.append(Track(r_st + position_offset, r_ei + position_offset, color="red"))

        for contig_id in final_parts:
            contig = contig_map[contig_id]
            if contig.group_ref != group_ref:
                continue

            if not isinstance(contig, AlignedContig):
                continue

            if contig_id in bad_contigs:
                continue

            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(contig)
            reference_tracks.append(Track(a_r_st + position_offset, a_r_ei + position_offset, color="yellow"))
            reference_min = min(a_r_st + position_offset, reference_min)
            reference_max = max(a_r_ei + position_offset, reference_max)

        for contig_id in final_parts:
            contig = contig_map[contig_id]
            if contig.group_ref != group_ref:
                continue

            if not isinstance(contig, AlignedContig):
                continue

            if contig_id in bad_contigs:
                continue

            (a_r_st, a_r_ei, f_r_st, f_r_ei) = get_contig_coordinates(contig)
            reference_tracks.append(Track(f_r_st + position_offset, f_r_ei + position_offset, color="lightgray"))
            reference_min = min(f_r_st + position_offset, reference_min)
            reference_max = max(f_r_ei + position_offset, reference_max)

        figure.add(Multitrack(reference_tracks))
        midpoint = round((reference_max - reference_min) / 2 + reference_min)
        figure.add(Track(midpoint, midpoint, label=group_ref, color="transparent", h=-11.5))

        ##########
        # Arrows #
        ##########

        ref_arrows: List[Arrow] = []
        for root in sorted_roots:
            parts_ids = final_children_mapping[root]
            parts_ids = [name for name in parts_ids if name not in bad_contigs]
            parts = [contig_map[name] for name in parts_ids]
            parts = [part for part in parts if part.group_ref == group_ref]
            ref_arrows.extend(get_arrows(parts, labels=True))

        if ref_arrows:
            figure.add(ArrowGroup(ref_arrows))

        ###########
        # Contigs #
        ###########

        for root in sorted_roots:
            parts_ids = final_children_mapping[root]
            parts_ids = [name for name in parts_ids if name not in bad_contigs]
            parts = [contig_map[name] for name in parts_ids]
            parts = [part for part in parts if part.group_ref == group_ref]
            if parts:
                figure.add(ArrowGroup(list(get_arrows(parts, labels=False))))
                figure.add(Multitrack(list(get_tracks(parts))))

        #############
        # Discarded #
        #############

        def get_group_discards(group_ref):
            for root in sorted_roots:
                if contig_map[root].group_ref != group_ref:
                    continue

                parts_ids = final_children_mapping[root]
                parts_ids = [id for id in parts_ids if id in discarded]
                unaligned_parts = merged_unaligned_parts.get(root, [])
                for id in sorted(parts_ids + unaligned_parts,
                                 key=lambda x: name_map[x.id] if isinstance(x, Contig) else name_map[x]):
                    if id in unaligned_parts:
                        (q_st, q_ei) = query_position_map[id]
                        label = name_map[id]
                        yield Track(position_offset, position_offset + abs(q_ei - q_st),
                                    label=label, color="yellow")
                    else:
                        part = contig_map[id]
                        yield Multitrack(list(get_tracks([part])))

        disc = list(get_group_discards(group_ref))
        if disc:
            add_section("discards:")
            for element in disc:
                figure.add(element)

        #############
        # Anomalies #
        #############

        def get_group_anomalies(group_ref):
            for root in sorted_roots:
                parts_ids = final_children_mapping[root]
                parts_ids = [name for name in parts_ids if name in anomaly]
                parts = [contig_map[name] for name in parts_ids]
                parts = [part for part in parts if part.group_ref == group_ref]
                for part in parts:
                    # Check if we have alignment data to display
                    has_alignment_data = False
                    if part.id in anomaly_data_map:
                        anomaly_event = anomaly_data_map[part.id]
                        if isinstance(anomaly_event, events.StrandConflict):
                            has_alignment_data = True
                            hits = anomaly_event.hits
                            strands = anomaly_event.strands
                            name = name_map.get(part.id, "")

                            # Generate arrows first (displayed on top)
                            arrows = get_anomaly_arrows(hits, strands, name)
                            if arrows:
                                yield ArrowGroup(arrows)

                            # Generate tracks with proper coloring (grey for mapped, yellow for unmapped)
                            yield Multitrack(list(get_anomaly_tracks(part, hits, name)))

                    # For anomalies without alignment data, display using default tracks
                    if not has_alignment_data:
                        yield Multitrack(list(get_tracks([part])))

        anom = list(get_group_anomalies(group_ref))
        if anom:
            add_section("anomaly:")
            for element in anom:
                figure.add(element)

    ###########
    # Unknown #
    ###########

    if unknown:
        add_section("unknown:")
        for parent_id in sorted_roots:
            parts_ids = final_children_mapping[parent_id]
            parts_ids = [name for name in parts_ids if name in unknown]
            parts = [contig_map[name] for name in parts_ids]
            for part in parts:
                figure.add(Multitrack(list(get_tracks([part]))))

    if not figure.elements:
        figure.add(Track(0, max_position, label='.', color='none'))
        figure.add(Track(0, max_position * 3 / 2, label='No contigs found.', color='none', h=-10))
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
