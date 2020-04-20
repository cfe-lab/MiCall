import os
from csv import DictReader
from operator import attrgetter
from pathlib import Path

from micall.utils.alignment_wrapper import align_nucs

try:
    # noinspection PyPackageRequirements
    from mappy import Aligner
except ImportError:
    Aligner = None

from micall.utils.fetch_sequences import fetch_by_accession

# Source: https://github.com/PoonLab/sam2conseq/wiki
MAPPING = """\
SRR11177792	RNA-Seq	2020-02-25	MT072688	sam2conseq reports	
SRR10903401	RNA-Seq	2020-01-18	MN988669.1	sam2conseq reports	
SRR10903402	RNA-Seq	2020-01-18	MN988668.1	sam2conseq reports	
SRR10971381	RNA-Seq	2020-01-27	MN908947	published sequence	
SRR11241254	RNA-Seq	2020-03-04	MT163716.1 or EPI_ISL_413025
SRR11241255	RNA-Seq	2020-03-04	MT163717.1	3 nt differences (
SRR11247075	RNA-Seq	2020-03-05	MT163721.1		
SRR11247076	RNA-Seq	2020-03-05	MT163720.1		
SRR11247077	RNA-Seq	2020-03-05	MT163719.1		
SRR11247078	RNA-Seq	2020-03-05	MT163718.1		
SRR11278090	RNA-Seq	2020-03-09	EPI_ISL_413601		
SRR11278091	RNA-Seq	2020-03-09	EPI_ISL_413563	Exact match except 3' end has extra "?"	Y
SRR11278092	RNA-Seq	2020-03-09	EPI_ISL_413562		
SRR11278164	RNA-Seq	2020-03-09	EPI_ISL_413653		
SRR11278165	RNA-Seq	2020-03-09	EPI_ISL_413652		
SRR11278166	RNA-Seq	2020-03-09	EPI_ISL_413651		
SRR11278167	RNA-Seq	2020-03-09	EPI_ISL_413650		
SRR11278168	RNA-Seq	2020-03-09	EPI_ISL_413649		
SRR11140744	WGS	2020-02-21	EPI_ISL_408670		
SRR11140746	WGS	2020-02-21	EPI_ISL_408670	Published has extra 3 A's at 3' end	
SRR11140748	WGS	2020-02-21	EPI_ISL_408670		
SRR11140750	WGS	2020-02-21	EPI_ISL_408670		
SRR11092057	RNA-Seq	2020-02-15	MN996528.1		
SRR11092058	RNA-Seq	2020-02-15	MN996527.1		
SRR11092064	RNA-Seq	2020-02-15	MN996531.1		
SRR11085733	RNA-Seq	2020-02-13	MN611525	non-human, failure to map	
SRR11085736	RNA-Seq	2020-02-13	MN611522	non-human, failure to map	
SRR11085737	RNA-Seq	2020-02-13	MN611521	non-human, failure to map	
SRR11085738	RNA-Seq	2020-02-13	MN611520	non-human, failure to map	
SRR11085740	RNA-Seq	2020-02-13	MN611518	non-human, failure to map	
SRR11085741	RNA-Seq	2020-02-13	MN611517	non-human, failure to map	
SRR11085797	RNA-Seq	2020-02-13	MN996532.1		
SRR11314339	RNA-Seq	2020-03-17	MT192765		
SRR11092056	RNA-Seq	2020-02-15	MN996530	Low coverage, sample is mostly human DNA (known issue)	Y
"""


def main():
    accessions = {}  # {reads_accession: conseq_accession}
    for line in MAPPING.splitlines():
        fields = line.split()
        reads_accession = fields[0]
        conseq_accession = fields[3]
        accessions[reads_accession] = conseq_accession

    source_path = Path(__file__)
    working_path = source_path.parent.parent / 'tests' / 'working'
    conseq_path = working_path / 'corona' / 'sars_results20200417' / 'conseq_all.csv'
    with conseq_path.open() as f:
        for row in DictReader(f):
            cutoff = row['consensus-percent-cutoff']
            if cutoff != 'MAX':
                continue
            sample_name = row['sample']
            reads_accession = sample_name.split('-')[0]
            conseq_accession = accessions[reads_accession]
            try:
                print('Fetching', conseq_accession)
                published_seq = fetch_by_accession(conseq_accession)
                print('Fetched', conseq_accession)
            except ValueError as ex:
                print(ex)
                continue
            if Aligner is not None:
                compare_using_mappy(row, conseq_accession, published_seq)
            else:
                compare_using_gotoh(row, conseq_accession, published_seq)
            print()


def compare_using_gotoh(row, conseq_accession, published_seq):
    conseq = row['sequence']
    sample_name = row['sample']
    print(sample_name, conseq_accession)
    print(len(conseq), len(published_seq))
    aln_pub_seq, aln_conseq, score = align_nucs(published_seq,
                                                conseq)
    screen_width = int(os.environ.get('COLUMNS', 100))
    print(aln_pub_seq[:screen_width])
    print(aln_conseq[:screen_width])
    mismatch_count = add_count = missing_count = 0
    for i, (theirs, ours) in enumerate(zip(aln_pub_seq, aln_conseq)):
        if theirs == ours:
            pass
        elif ours == '-':
            missing_count += 1
        elif theirs == '-':
            add_count += 1
        else:
            mismatch_count += 1
            if mismatch_count < 100:
                print(f'{i}: {theirs} => {ours}')
    print(f'{mismatch_count} mismatches, {missing_count} missing, and '
          f'{add_count} added out of {len(published_seq)}.')


def compare_using_mappy(row, conseq_accession, published_seq):
    conseq = row['sequence']
    sample_name = row['sample']
    print(sample_name, conseq_accession)
    print(len(conseq), len(published_seq))
    conseq = 'x' * int(row['offset']) + conseq
    # noinspection PyCallingNonCallable
    aligner = Aligner(seq=published_seq)
    prev_ref_pos = prev_query_pos = 0
    for hit in sorted(aligner.map(conseq), key=attrgetter('q_st')):  # traverse alignments
        print(published_seq[prev_ref_pos:hit.r_st])
        print(conseq[prev_query_pos:hit.q_st])
        print(published_seq[hit.r_st:hit.r_en])
        print(conseq[hit.q_st:hit.q_en])
        prev_ref_pos = hit.r_en
        prev_query_pos = hit.q_en
        print(hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.cigar_str, hit.is_primary)
        print(hit.blen, hit.mlen, hit.NM, hit.mlen / (hit.mlen + hit.NM))


main()
