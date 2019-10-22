from base64 import standard_b64encode
from turtle import Turtle
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import defaultdict, Counter
from csv import DictReader

from drawSvg import Line, Lines, Drawing

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
    f2 = build_coverage_figure(args.genome_coverage2_csv)
    # del f2.elements[7:]
    contig_y, contig = f2.elements[6]

    f.h += 300
    f.w = max(f.w, f2.w)
    f.elements.extend(((y + 120 + (100 if i > 1 else 0), track)
                       for i, (y, track) in enumerate(f2.elements[4:])))

    if __name__ != '__live_coding__':
        drawing_width = 970
    else:
        # noinspection PyProtectedMember
        turtle_screen = Turtle._screen
        drawing_width = turtle_screen.cv.cget('width')
        diagram_path = None

    drawing = f.show(w=drawing_width)
    seed1_y = f.h - f.elements[5][0] - 10
    seed2_y = f.h - f.elements[7][0] + 25
    ref_y = f.h - f.elements[7][0]
    contig_y = f.h - f.elements[8][0] + 10
    contig_x = contig.tracks[0].a
    x_scale = drawing_width / f.w
    blast_display = BlastDisplay(drawing, x_scale, ref_y, contig_x, contig_y)
    blast_rows = list(DictReader(args.blast2_csv))
    blast_rows.reverse()
    best_score = 0
    best_ref = None
    for row in blast_rows:
        score = int(row['score'])
        if score > best_score:
            best_score = score
            best_ref = row['ref_name']
    for row in blast_rows:
        if row['ref_name'] == best_ref:
            start = int(row['start'])
            end = int(row['end'])
            ref_start = int(row['ref_start'])
            ref_end = int(row['ref_end'])
            blast_display.add_match(start, end, ref_start, ref_end)

    for source_bin, bin_moves in sorted_moves:
        total = sum(bin_moves.values())
        for i, (target_bin, count) in enumerate(bin_moves.most_common()):
            fraction = count / total
            if fraction < 0.1:
                break
            if i == 0:
                print(source_bin, end=':')
            print(f'\t{target_bin}({fraction})')
            if source_bin is None:
                source_x = target_bin
                source_y = seed2_y - 10
            else:
                source_x = source_bin
                source_y = seed1_y
            if target_bin is None:
                target_x = source_bin
                target_y = seed1_y + 10
            else:
                target_x = target_bin
                target_y = seed2_y
            source_x *= 100*x_scale
            target_x *= 100*x_scale
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

    def add_match(self, start, end, ref_start, ref_end):
        left1 = ref_start * self.x_scale
        top = self.ref_y
        bottom = self.contig_y
        left2 = (self.contig_x + start) * self.x_scale
        right1 = ref_end * self.x_scale
        right2 = (self.contig_x + end) * self.x_scale
        self.drawing.append(Lines(left1, top,
                                  left2, bottom,
                                  right2, bottom,
                                  right1, top,
                                  left1, top,
                                  stroke='black',
                                  fill='ivory'))


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


main()
