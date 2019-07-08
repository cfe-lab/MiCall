from genetracks import Figure, Track, Alignment, Multitrack, Label, Coverage
import sys
import pysam

subtyping = open("/tmp/subtyping.csv", 'w')

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


def get_alignments(fn, references, label="", coverage=None):
    bf = pysam.AlignmentFile(fn)

    contigs = {}
    ref_names = {}

    empty_tracks = []
    seen = set()
    for aln in bf:
        if aln.query_name not in coverage:
            coverage[aln.query_name] = [0 for _ in range(0, 10000)]
        print(aln.reference_name)
        if not aln.reference_name:
            if aln.query_name in seen:
                continue
            seen.add(aln.query_name)

            empty_tracks.append("Non-HCV: {} {}, length: {}, max depth: {}".format(label,
                aln.query_name, aln.query_length,
                max(coverage[aln.query_name])))
            continue
        if "HCV" not in aln.reference_name:
            if aln.query_name in seen:
                continue
            seen.add(aln.query_name)
            empty_tracks.append("Non-HCV: {} {} ({}), length: {}, max depth: {}".format(label,
                aln.query_name, aln.reference_name, aln.query_length,
                max(coverage[aln.query_name])))
 
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
    coverage_tracks = []
    for contig in sorted(contigs):
        coverage_tracks.append(Coverage(contigs[contig][0].a,
            contigs[contig][0].b, coverage[contig]))
        contigs[contig].append(Track(4000,4000, label="{} {}: {}".format(label,
            contig, ref_names[contig])))
        tracks.append(Multitrack(contigs[contig], join=False))
        print("{},{}".format(contig, ref_names[contig]), file=subtyping)
    return tracks, coverage_tracks, empty_tracks

def get_coverage(fn, length=10000):
    cs = {}
    bf = pysam.AlignmentFile(fn)

    for pc in bf.pileup():
        if pc.reference_name not in cs:
            cs[pc.reference_name] = [0 for _ in range(0, length)]
        if pc.pos < 1 or pc.pos >= length:
            continue
        cs[pc.reference_name][pc.pos] = pc.n
    return cs

if __name__ == "__main__":
    wg_cov = get_coverage(sys.argv[3])
    mid_cov = get_coverage(sys.argv[4])
    label = sys.argv[5]
    
    ts, cts, es1 = get_alignments(sys.argv[1], ["NC_004102.1"],
            label="{}-WG".format(label), coverage=wg_cov)

    f = Figure()
    hcv_region_tracks(f)

    for t, c in zip(ts, cts):
        f.add(c, gap=-5)
        f.add(t)

    tx, ctx, es2 = get_alignments(sys.argv[2], ["NC_004102.1"],
            label="{}-MIDI".format(label), coverage=mid_cov)

    for t, c in zip(tx, ctx):
        f.add(c, gap=-5)
        f.add(t)

    for e in es1 + es2:
        f.add(Track(4000,4000, label=e))

    d = f.show(w=970)
    d.saveSvg("/tmp/alignment.svg")
    d.savePng("/tmp/alignment.png")
    
subtyping.close()
