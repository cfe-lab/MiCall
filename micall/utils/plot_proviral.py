from genetracks import Figure, Track, Alignment, Multitrack, Label
import sys
import pysam

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
        f.add(Multitrack([Track(l, r, label=Label(0, text, offset=1),
            color=color) for l, r, text, color in reading_frame]), gap=0)

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
                          for l, r, text, color in reading_frame]), gap=0)


def get_alignments(fn, references, label=""):
    bf = pysam.AlignmentFile(fn)

    contigs = {}
    ref_names = {}
    for aln in bf:
        print(aln.reference_name)
        if not aln.reference_name:
            continue
        if aln.reference_name not in references:
            continue
        if aln.cigartuples is None:
            continue

        contig = aln.query_name
        if contig not in contigs:
            contigs[contig] = []
            ref_names[contig] = aln.reference_name

        for start, end in aln.get_blocks():
            print(start, end)
            contigs[contig].append(Track(start, end))
    bf.close()

    tracks = []
    for contig in sorted(contigs):
        tracks.append(Multitrack([Multitrack(contigs[contig], join=True),
            Track(4000,4000, label="{}-{}: {}".format(label, contig,
                ref_names[contig]))], join=False))
    return tracks

if __name__ == "__main__":
    ts = get_alignments(sys.argv[1], ["K03455.1"], label="WG")

    f = Figure()
    hiv_region_tracks(f)

    for t in ts:
        f.add(t)

    f.show(w=970).saveSvg("/tmp/alignment.svg")
    
