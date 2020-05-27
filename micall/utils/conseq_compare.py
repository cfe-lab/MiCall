import csv
import os
import sys
from csv import DictReader
from operator import attrgetter
from pathlib import Path
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.utils.alignment_wrapper import align_nucs

try:
    # noinspection PyPackageRequirements
    from mappy import Aligner
except ImportError:
    Aligner = None

from micall.utils.fetch_sequences import fetch_by_accession


# Dan stuff
import sys
from micall.core.project_config import ProjectConfig

REFERENCE = ProjectConfig.loadDefault()
REFERENCE = REFERENCE.getReference('SARS-CoV-2-seed')

def load_coverage(csv):
    result = {}
    with open(csv) as csvfile:
        reader = DictReader(csvfile)
        for row in reader:
            result[int(row['query_nuc_pos'])] = int(row['coverage'])
    return result

BATCH = 'batch_01'

ROOT = (
    Path('/wow')
    / BATCH
    / 'Results'
    / 'scratch'
)

DATA = (
    Path(os.path.realpath(__file__)).parent.parent
    / 'data'
)


def load_yaml(yaml_file):
    with open(yaml_file) as f:
        return yaml.load(f, Loader=yaml.SafeLoader)
# end Dan stuff

# Source: https://github.com/PoonLab/sam2conseq/wiki
MAPPING = """\
SRR11177792	RNA-Seq	2020-02-25	MT072688	sam2conseq reports
	RNA-Seq	2020-01-18	MN988669.1	sam2conseq reports
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
SRR11578349_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_427024
SRR11578341_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426901
SRR11578342_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426900
SRR11578343_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426899
SRR11578344_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426899
SRR11578345_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426656
SRR11578346_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_426898
SRR11578347_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_427026
SRR11578348_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_427025
SRR11578349_1.fastq	RNA-Seq	2020-05-04	EPI_ISL_427024
SRR11593354_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414507
SRR11593355_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414574
SRR11593356_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414509
SRR11593357_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414508
SRR11593358_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414506
SRR11593359_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414505
SRR11593360_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414504
SRR11593361_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414499
SRR11593362_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414498
SRR11593364_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_414497
SRR11593365_1.fastq	RNA-Seq	2020-05-08	EPI_ISL_413488
"""


def main():
    accessions = {}  # {reads_accession: conseq_accession}
    for line in MAPPING.splitlines():
        fields = line.split()
        reads_accession = fields[0]
        conseq_accession = fields[3]
        accessions[reads_accession] = conseq_accession

    summaries = []
    source_path = Path(__file__)
    working_path = source_path.parent.parent.parent / 'micall' / 'tests' / 'working'
    conseq_path = working_path / f'{BATCH}_conseq_all.csv'
    accessions_path = working_path / 'fetched_accessions.fasta'
    accessions_file = accessions_path.open('a')
    accessions_index = None
    with conseq_path.open() as f:
        for row in DictReader(f):
            cutoff = row['consensus-percent-cutoff']
            if cutoff != 'MAX':
                continue
            region = row['region']
            if region != '':
                continue
            sample_name = row['sample']
            reads_accession = sample_name.split('-')[0]
            conseq_accession = accessions[reads_accession]
            if accessions_index is None:
                accessions_index = SeqIO.index(str(accessions_path), 'fasta')

            try:
                published_seq = str(accessions_index[conseq_accession].seq)
            except KeyError:
                try:
                    print('Fetching', conseq_accession)
                    published_seq = fetch_by_accession(conseq_accession)
                    print('Fetched', conseq_accession)
                except ValueError as ex:
                    print(ex)
                    published_seq = ''
                SeqIO.write([SeqRecord(Seq(published_seq), id=conseq_accession)],
                            accessions_file,
                            'fasta')
                accessions_file.flush()
                accessions_index = None
            if published_seq == '':
                print(f'Skipping conseq {conseq_accession}.')
                summary = f'{sample_name}|{conseq_accession}|Failed to fetch comparison sequence.'
                summaries.append(summary)
                continue
            if Aligner is not None:
                compare_using_mappy(row, conseq_accession, published_seq)
            else:
                summary = compare_using_gotoh(row, conseq_accession, published_seq)
                summaries.append(summary)
            print()
    print('Run|Compared to|Differences')
    print('---|-----------|-----------')
    print(os.linesep.join(summaries))


def compare_using_gotoh(row, conseq_accession, published_seq):

    # try:
    #     nlocs = yaml.load(open(DATA / f'{conseq_accession}.yaml'), Loader=yaml.SafeLoader)
    # except Exception as e:
    #     print(e)
    #     nlocs = []
    conseq = row['sequence'].replace('x', '')
    sample_name = row['sample']
    n_counts = load_yaml(DATA / 'n_counts.yaml')
    header = f'{conseq_accession} => {sample_name}'
    print(header)
    print(len(conseq), len(published_seq))
    aln_pub_seq, aln_conseq, score = align_nucs(published_seq,
                                                conseq)
    screen_width = int(os.environ.get('COLUMNS', 100))
    print(aln_pub_seq[:screen_width])
    print(aln_conseq[:screen_width])
    print('...', aln_pub_seq[-screen_width+4:])
    print('...', aln_conseq[-screen_width+4:])
    mismatch_count = add_count = missing_count = matching_count = 0

    mycsv = ROOT / sample_name / 'genome_coverage.csv'
    conseq_csv = ROOT / sample_name / 'conseq_all.csv'
    offset = None
    with open(conseq_csv) as csvfile:
        reader = DictReader(csvfile)
        for row in reader:
            offset = int(row['seed-offset'])
            break
    coverages = load_coverage(mycsv)
    # matchfile = open(ROOT / sample_name / 'mymatches.txt', 'w')
    # mismatchfile = open(ROOT / sample_name / 'mymismatches.txt', 'w')
    # addedfile = open(ROOT / sample_name / 'myadded.txt', 'w')
    datafile = open(ROOT / sample_name / 'data.csv', 'w')
    header = [
        'sample',
        'type',
        'pos',
        'micall_base',
        'gisaid_base',
        'ref_base',
        # 'is_nloc',
        'coverage'
        # 'mean_quality'
    ]
    writer = csv.DictWriter(datafile, header)
    writer.writeheader()
    for i, (theirs, ours) in enumerate(zip(aln_pub_seq, aln_conseq)):
        pos = i + 1 + offset
        # is_nloc = True if pos in nlocs else False
        coverage = '0'
        refbase = '?'
        try:
            refbase = REFERENCE[pos - 1]
        except IndexError:
            pass
        try:
            coverage = str(coverages[pos])
        except KeyError:
            pass
        pos = str(pos)
        if theirs == ours:
            # try:
            #     matchfile.write(str(coverages[i]) + '\n')
            # except KeyError:
            #     matchfile.write('0' + '\n')
            matching_count += 1
        elif ours == '-':
            my_type = 'deletion'
            writer.writerow({
                'sample': sample_name,
                'type': my_type,
                'pos': pos,
                'micall_base': ours,
                'gisaid_base': theirs,
                'ref_base': refbase,
                # 'is_nloc': is_nloc,
                'coverage': coverage
            })
            missing_count += 1
        elif theirs == '-':
            my_type = 'addition'
            writer.writerow({
                'sample': sample_name,
                'type': my_type,
                'pos': pos,
                'micall_base': ours,
                'gisaid_base': theirs,
                'ref_base': refbase,
                # 'is_nloc': is_nloc,
                'coverage': coverage
            })
            add_count += 1
            # try:
            #     addedfile.write(str(coverages[i]) + '\n')
            # except KeyError:
            #     addedfile.write('0' + '\n')
        else:
            my_type = 'mismatch'
            writer.writerow({
                'sample': sample_name,
                'type': my_type,
                'pos': pos,
                'micall_base': ours,
                'gisaid_base': theirs,
                'ref_base': refbase,
                # 'is_nloc': is_nloc,
                'coverage': coverage
            })
            # try:
            #     mismatchfile.write(str(coverages[i]) + '\n')
            # except KeyError:
            #     mismatchfile.write('0' + '\n')
            mismatch_count += 1
            if mismatch_count < 100:
                print(f'{i}: {theirs} => {ours}')
    summary = (f'{mismatch_count} mismatches, {missing_count} missing, and '
               f'{add_count} added out of {len(published_seq)}.')
    # matchfile.close()
    # mismatchfile.close()
    # addedfile.close()
    print(summary)

    with open(ROOT / sample_name / 'data_extended.csv', 'w') as o:
        header = [
            'sample',
            'conseq_len',
            'ref_len',
            'matches',
            'mismatches',
            'missing',
            'added',
            'num_n_removed',
            'concordance'
        ]
        csv_writer = csv.DictWriter(o, header)
        csv_writer.writeheader()

        # Compute data
        conseq_len = len(conseq)
        ref_len = len(published_seq)
        disagree = sum([mismatch_count, missing_count, add_count])
        removed_n = n_counts[conseq_accession]
        concordance = matching_count / ref_len * 100
        csv_writer.writerow({
            'sample': sample_name,
            'conseq_len': conseq_len,
            'ref_len': ref_len,
            'matches': matching_count,
            'mismatches': mismatch_count,
            'missing': missing_count,
            'added': add_count,
            'num_n_removed': n_counts[conseq_accession],
            'concordance': round(concordance, 2)
        })

    return f'{sample_name}|{conseq_accession}|{summary}'


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
