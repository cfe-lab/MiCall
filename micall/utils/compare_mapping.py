from base64 import standard_b64encode
from math import copysign
from turtle import Turtle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import defaultdict, Counter
from csv import DictReader

from drawSvg import Line, Lines, Drawing, Text, Circle

from micall.core.plot_contigs import build_coverage_figure
from micall.core.remap import SAM_FLAG_IS_UNMAPPED, SAM_FLAG_IS_FIRST_SEGMENT


def parse_args():
    parser = ArgumentParser(
        description='Compare two sets of mapped reads to see where they mapped.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('remap1_csv',
                        type=FileType(),
                        help='mapped reads from one run')
    parser.add_argument('remap2_csv',
                        type=FileType(),
                        help='mapped reads from the other run')
    parser.add_argument('genome_coverage1_csv',
                        type=FileType(),
                        help='genome coverage levels from one run')
    parser.add_argument('genome_coverage2_csv',
                        type=FileType(),
                        help='genome coverage levels from the other run')
    parser.add_argument('blast2_csv',
                        type=FileType(),
                        help='blast results from the other run')
    parser.add_argument('compare_mapping_svg',
                        help='diagram file to write')
    return parser.parse_args()


def main():
    args = parse_args()
    source_bins = load_read_bins(args.remap1_csv)
    target_bins = load_read_bins(args.remap2_csv)
    all_keys = set(source_bins)
    all_keys.update(target_bins)
    moves = defaultdict(Counter)  # {source_bin: {target_bin: count}}
    for read_key in all_keys:
        source_bin = source_bins.get(read_key)
        target_bin = target_bins.get(read_key)
        bin_moves = moves[source_bin]
        bin_moves[target_bin] += 1

    old_unmapped = moves.pop(None, None)
    sorted_moves = sorted(moves.items())
    if old_unmapped:
        sorted_moves.insert(0, (None, old_unmapped))

    diagram_path = args.compare_mapping_svg
    f = build_coverage_figure(args.genome_coverage1_csv)
    del f.elements[6:]
    ref_y, ref = f.elements[5]
    f2 = build_coverage_figure(args.genome_coverage2_csv)
    coverage3_y, coverage3 = f2.elements[4]
    contig3_y, contig3 = f2.elements[5]
    coverage1_y, coverage1 = f2.elements[6]
    contig1_y, contig1 = f2.elements[7]
    dashes_y, dashes = f2.elements[10]
    del f2.elements[4:]
    # contig_y, contig = f2.elements[6]

    # f.h += 50
    f.w = max(f.w, f2.w)
    f.add(coverage3, gap=-4)
    f.add(contig3, gap=30)
    contig3_y = f.elements[-1][0]
    coverage1.a = 0
    f.add(coverage1, gap=-4)
    contig1_shift = contig1.tracks[0].a
    for track in contig1.tracks:
        if track.a >= contig1_shift:
            track.a -= contig1_shift
            track.b -= contig1_shift
    f.add(contig1)
    contig1_y = f.elements[-1][0]

    if __name__ != '__live_coding__':
        drawing_width = 970
    else:
        # noinspection PyProtectedMember
        turtle_screen = Turtle._screen
        drawing_width = turtle_screen.cv.cget('width') - 10
        diagram_path = None

    drawing = f.show(w=drawing_width)
    seed1_y = f.h - f.elements[5][0] - 10
    seed2_y = f.h - f.elements[7][0] + 25
    ref_y = f.h - ref_y
    contig_y = f.h - contig1_y + 25
    contig_x = 0
    x_scale = drawing_width / f.w
    blast_display = BlastDisplay(drawing, x_scale, ref_y, contig_x, contig_y)
    blast_rows = list(DictReader(args.blast2_csv))
    blast_rows.sort(key=lambda match: int(match['score']), reverse=True)
    best_ref = None
    matched_positions = set()
    for row in blast_rows:
        if '003' not in row['contig_name']:
            continue
        if best_ref is None:
            best_ref = row['ref_name']
        elif row['ref_name'] != best_ref:
            continue
        start = int(row['start'])
        end = int(row['end'])
        new_positions = set(range(start, end))
        collisions = matched_positions & new_positions
        collision_fraction = len(collisions) / len(new_positions)
        if collision_fraction < 0.1:
            matched_positions |= new_positions
            ref_start = int(row['ref_start'])
            ref_end = int(row['ref_end'])
            blast_display.add_match(start, end, ref_start, ref_end)

    blast_display.draw()

    for source_bin, bin_moves in sorted_moves:
        total = sum(bin_moves.values())
        for i, (target_bin, count) in enumerate(bin_moves.most_common()):
            fraction = count / total
            if fraction < 0.1:
                break
            if i == 0:
                print(source_bin, end=':')
            print(f'\t{target_bin}({fraction})')
            source_shift = 0
            target_shift = contig_x*x_scale
            if source_bin is None:
                source_x = target_bin
                source_y = seed2_y - 10
                source_shift = target_shift
            else:
                source_x = source_bin
                source_y = seed1_y
            if target_bin is None:
                target_x = source_bin
                target_y = seed1_y + 10
                target_shift = source_shift
            else:
                target_x = target_bin
                target_y = seed2_y
            source_x *= 100*x_scale
            target_x *= 100*x_scale
            source_x += source_shift
            target_x += target_shift
            drawing.append(Line(source_x, source_y,
                                target_x, target_y,
                                stroke='black',
                                stroke_width=2*fraction,
                                stroke_opacity=0.25))

    if diagram_path is not None:
        drawing.saveSvg(diagram_path)
    else:
        display_image(drawing)


class BlastDisplay:
    def __init__(self,
                 drawing: Drawing,
                 x_scale: float,
                 ref_y: int,
                 contig_x: int,
                 contig_y: int):
        self.drawing = drawing
        self.x_scale = x_scale
        self.ref_y = ref_y
        self.contig_x = contig_x
        self.contig_y = contig_y
        self.matches = []  # [(start, end, ref_start, ref_end)]

    def add_match(self, start, end, ref_start, ref_end):
        self.matches.append((start, end, ref_start, ref_end))

    def draw(self, use_arrows: bool = True):
        if use_arrows:
            self.matches.sort()
        for i, (start, end, ref_start, ref_end) in enumerate(self.matches, 1):
            left1 = ref_start * self.x_scale
            left2 = (self.contig_x + start) * self.x_scale
            right1 = ref_end * self.x_scale
            right2 = (self.contig_x + end) * self.x_scale
            if use_arrows:
                self.draw_arrows(left1, right1, left2, right2, str(i))
            else:
                self.draw_region(left1, right1, left2, right2)

    def draw_region(self, left1, right1, left2, right2):
        top = self.ref_y
        bottom = self.contig_y
        self.drawing.append(Lines(left1, top,
                                  left2, bottom,
                                  right2, bottom,
                                  right1, top,
                                  left1, top,
                                  stroke='black',
                                  fill='ivory'))

    def draw_arrows(self, left1, right1, left2, right2, label):
        self.draw_arrow(left1, right1, self.ref_y - 15, label)
        self.draw_arrow(left2, right2, self.contig_y + 15, label)

    def draw_arrow(self, start, end, y, label):
        r = 10
        arrow_size = 7
        direction = copysign(1, end-start)
        centre = (start + end - direction*arrow_size)/2
        centre_start = centre - direction*r
        centre_end = centre + direction*r
        self.drawing.append(Line(start, y,
                                 centre_start, y,
                                 stroke='black'))
        self.drawing.append(Line(centre_end, y,
                                 end, y,
                                 stroke='black'))
        self.drawing.append(Circle(centre, y, r, fill='ivory', stroke='black'))
        self.drawing.append(Text(label,
                                 15,
                                 centre, y,
                                 text_anchor='middle',
                                 dy="0.3em"))
        arrow_end = end
        arrow_start = arrow_end - arrow_size*direction
        self.drawing.append(Lines(arrow_end, y,
                                  arrow_start, y + arrow_size/2,
                                  arrow_start, y - arrow_size/2,
                                  fill='black'))


def display_image(drawing):
    png_data = drawing.rasterize()
    encoded = standard_b64encode(png_data.pngData)
    image_text = str(encoded.decode('UTF-8'))
    # noinspection PyUnresolvedReferences
    Turtle.display_image(0, 0, image=image_text)


def load_read_bins(remap_csv):
    reader = DictReader(remap_csv)
    bins = {}  # {qname: bin_num}
    for i, row in enumerate(reader):
        if i == 10 and __name__ == '__live_coding__':
            break
        qname = row['qname']
        flag = int(row['flag'])
        if flag & SAM_FLAG_IS_UNMAPPED == 0:
            read_num = 1 if flag & SAM_FLAG_IS_FIRST_SEGMENT else 2
            pos = int(row['pos'])
            bin_num = pos // 100
            bins[(qname, read_num)] = bin_num
    return bins


if __name__ == '__main__':
    main()
