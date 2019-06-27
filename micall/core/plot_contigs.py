from csv import DictReader
from operator import itemgetter

from genetracks import Figure, Track, Alignment, Multitrack, Label


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
    for reading_frame in [first,]:
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


def plot_contigs(contigs_csv, contigs_svg_path):
    rows = list(DictReader(contigs_csv))
    for row in rows:
        row['match'] = float(row['match'])
        for field in ('start', 'end', 'ref_start', 'ref_end'):
            field_value = row[field]
            row[field] = int(field_value) if field_value else None
    rows.sort(key=itemgetter('genotype', 'ref_start'))
    f = Figure()
    is_empty = True
    genotype = None
    for row in rows:
        new_genotype = row['genotype']
        if genotype != new_genotype:
            if new_genotype.startswith('HIV'):
                hiv_region_tracks(f)
            elif new_genotype.startswith('HCV'):
                hcv_region_tracks(f)
            else:
                continue
            genotype = new_genotype
        f.add(Track(row['ref_start'], row['ref_end'], label=genotype))
        is_empty = False
    if is_empty:
        f.add(Track(4000, 4000, label='No contigs assembled.'))
    f.show(w=970).saveSvg(contigs_svg_path)