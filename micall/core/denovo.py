import argparse
import logging
import os
import typing
from collections import Counter
from csv import DictWriter, DictReader
from datetime import datetime
from glob import glob
from io import StringIO
from itertools import groupby
from operator import itemgetter
from shutil import rmtree
from subprocess import run, PIPE, CalledProcessError, STDOUT
from tempfile import mkdtemp

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.core.project_config import ProjectConfig

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from iva.assembly import Assembly
from pyfastaq.tasks import deinterleave
from micall.core.remap import remap, map_to_contigs

HAPLOFLOW = "haploflow"
IVA = "iva"
RAGTAG = "/home/charlotte/Documents/Git/MiCall/.venv/bin/RagTag/ragtag.py"
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs.fasta')
logger = logging.getLogger(__name__)


def write_contig_refs(contigs_fasta_path,
                      contigs_csv,
                      merged_contigs_csv=None,
                      blast_csv=None):
    """ Run BLAST search to identify contig sequences.

    :param str contigs_fasta_path: path to file to read contig sequences from
        and append merged contigs to
    :param contigs_csv: open file to write assembled contigs to
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    """
    writer = DictWriter(contigs_csv,
                        ['ref', 'match', 'group_ref', 'contig'],
                        lineterminator=os.linesep)
    writer.writeheader()
    with open(contigs_fasta_path, 'a') as contigs_fasta:
        if merged_contigs_csv is not None:
            contig_reader = DictReader(merged_contigs_csv)
            for i, row in enumerate(contig_reader, 1):
                contig_name = f'merged-contig-{i}'
                contigs_fasta.write(f">{contig_name}\n{row['contig']}\n")
    group_refs = {}
    genotypes = genotype(contigs_fasta_path,
                         blast_csv=blast_csv,
                         group_refs=group_refs)
    genotype_count = 0
    for i, record in enumerate(SeqIO.parse(contigs_fasta_path, "fasta")):
        (ref_name, match_fraction) = genotypes.get(record.name, ('unknown', 0))
        seq = record.seq
        if match_fraction < 0:
            seq = seq.reverse_complement()
            match_fraction *= -1
        writer.writerow(dict(ref=ref_name,
                             match=match_fraction,
                             group_ref=group_refs.get(ref_name),
                             contig=seq))
        genotype_count += 1
    return genotype_count


def genotype(fasta, db=DEFAULT_DATABASE, blast_csv=None, group_refs=None):
    """ Use Blastn to search for the genotype of a set of reference sequences.

    :param str fasta: file path of the FASTA file containing the query
        sequences
    :param str db: file path of the database to search for matches
    :param blast_csv: open file to write the blast matches to, or None
    :param dict group_refs: {contig_ref: group_ref} or None. The dictionary
        will get filled in with the mapping from each contig's reference name
        to the best matched reference for the whole seed group.
    :return: {query_name: (ref_name, matched_fraction)} where query_name is a
        sequence header from the query sequences FASTA file, ref_name is the
        name of the best match from the database, and matched_fraction is the
        fraction of the query that aligned against the reference (matches and
        mismatches).
    """
    contig_nums = {}  # {contig_name: contig_num}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_name = line[1:-1]
                contig_nums[contig_name] = len(contig_nums) + 1
    blast_columns = ['qaccver',
                     'saccver',
                     'pident',
                     'score',
                     'qcovhsp',
                     'qstart',
                     'qend',
                     'sstart',
                     'send']
    cline = NcbiblastnCommandline(query=fasta,
                                  db=db,
                                  outfmt=f'"10 {" ".join(blast_columns)}"',
                                  evalue=0.0001,
                                  gapopen=5,
                                  gapextend=2,
                                  penalty=-3,
                                  reward=1,
                                  max_target_seqs=5000)
    stdout, _ = cline()
    samples = {}  # {query_name: (subject_name, matched_fraction)}
    matches = sorted(DictReader(StringIO(stdout), blast_columns),
                     key=lambda row: (row['qaccver'], float(row['score'])))
    if not blast_csv:
        blast_writer = None
    else:
        blast_writer = DictWriter(blast_csv,
                                  ['contig_num',
                                   'ref_name',
                                   'score',
                                   'match',
                                   'pident',
                                   'start',
                                   'end',
                                   'ref_start',
                                   'ref_end'],
                                  lineterminator=os.linesep)
        blast_writer.writeheader()
    contig_top_matches = {match['qaccver']: match['saccver']
                          for match in matches}
    top_refs = set(contig_top_matches.values())
    projects = ProjectConfig.loadDefault()
    match_scores = Counter()
    for contig_name, contig_matches in groupby(matches, itemgetter('qaccver')):
        contig_top_ref = contig_top_matches[contig_name]
        contig_seed_group = projects.getSeedGroup(contig_top_ref)
        for match in contig_matches:
            ref_name = match['saccver']
            if ref_name not in top_refs:
                continue
            match_seed_group = projects.getSeedGroup(ref_name)
            if match_seed_group == contig_seed_group:
                match_scores[ref_name] += float(match['score'])

    if group_refs is not None:
        group_top_refs = {projects.getSeedGroup(ref_name): ref_name
                          for ref_name, count in reversed(match_scores.most_common())}
        for ref_name in contig_top_matches.values():
            group_refs[ref_name] = group_top_refs[projects.getSeedGroup(ref_name)]

    for match in matches:
        matched_fraction = float(match['qcovhsp']) / 100
        if int(match['send']) < int(match['sstart']):
            matched_fraction *= -1
        pident = round(float(match['pident']))
        contig_name = match['qaccver']
        samples[contig_name] = (match['saccver'], matched_fraction)
        if blast_writer:
            blast_writer.writerow(dict(contig_num=contig_nums[contig_name],
                                       ref_name=match['saccver'],
                                       score=match['score'],
                                       match=matched_fraction,
                                       pident=pident,
                                       start=match['qstart'],
                                       end=match['qend'],
                                       ref_start=match['sstart'],
                                       ref_end=match['send']))
    return samples


def separate_contigs(contigs_csv, ref_contigs_csv, noref_contigs_csv):
    """ Separate contigs into those that mapped to or did not map to a reference.
    :param contigs_csv:         file with contigs, open in read mode
    :param ref_contigs_csv:     file for contigs that mapped to a reference, open in write mode
    :param noref_contigs_csv:   file for contigs that did not map to a reference, open in write mode
    """
    threshold = 0.1
    # is a match threshold sufficient or do we need info from blast_csv as well?
    fieldnames = ['ref', 'match', 'group_ref', 'contig']
    ref_contig_writer = DictWriter(ref_contigs_csv, fieldnames)
    ref_contig_writer.writeheader()
    noref_contig_writer = DictWriter(noref_contigs_csv, fieldnames)
    noref_contig_writer.writeheader()
    contig_reader = DictReader(contigs_csv)
    num_total = 0
    num_match = 0
    for row in contig_reader:
        num_total += 1
        if float(row['match']) > threshold:
            ref_contig_writer.writerow(row)
            num_match += 1
        else:
            noref_contig_writer.writerow(row)
    return num_total - num_match


def separate_reads(remap_csv, ref_reads_file, noref_reads_file, unmapped1, unmapped2):
    """ Separate reads from remap.csv file into those that mapped to un unknown partial and the rest.

    :param remap_csv: remap output file created by map_to_contigs, open in read mode
    :param ref_reads_file: file to write potentially useful reads (that mapped to useful contigs or that did not map)
    :param noref_reads_file: file to write useless reads (that mapped to unknown contig)
    :param unmapped1: fasta file 1 of reads that did not map
    :param unmapped2: fasta file 2 of reads that did not map
    """
    fieldnames = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual']
    remap_reader = DictReader(remap_csv)
    for row in remap_reader:
        if row['rname'][-16:] == "-unknown-partial":
            file_to_write = noref_reads_file
        else:
            file_to_write = ref_reads_file
        file_to_write.write('@'+row['qname']+'\n')
        file_to_write.write(row['seq']+'\n')
        file_to_write.write('+\n')
        file_to_write.write(row['qual']+'\n')
    for line in unmapped1:
        ref_reads_file.write(line)
    for line in unmapped2:
        ref_reads_file.write(line)


def denovo(fastq1_path: str,
           fastq2_path: str,
           contigs_csv: typing.TextIO,
           work_dir: str = '.',
           merged_contigs_csv: typing.TextIO = None,
           blast_csv: typing.TextIO = None,
           haplo_args=None):
    """ Use de novo assembly to build contigs from reads.

    :param fastq1_path: FASTQ file name for read 1 reads
    :param fastq2_path: FASTQ file name for read 2 reads
    :param contigs_csv: open file to write assembled contigs to
    :param work_dir: path for writing temporary files
    :param merged_contigs_csv: open file to read contigs that were merged from
        amplicon reads
    :param blast_csv: open file to write BLAST search results for each contig
    """
    old_tmp_dirs = glob(os.path.join(work_dir, 'assembly_*'))
    for old_tmp_dir in old_tmp_dirs:
        rmtree(old_tmp_dir, ignore_errors=True)

    tmp_dir = mkdtemp(dir=work_dir, prefix='assembly_')
    start_time = datetime.now()
    start_dir = os.getcwd()
    joined_path = os.path.join(tmp_dir, 'joined.fastq')
    run(['merge-mates',
         fastq1_path,
         fastq2_path,
         '--interleave',
         '-o', joined_path],
        check=True)

    if haplo_args is None:
        haplo_args = {'long': 0,
                      'filter': 500,
                      'thres': -1,
                      'strict': 5,
                      'error': 0.02,
                      'kmer': 41,
                      'merge': False,
                      'scaffold': False,
                      'patch': False,
                      'ref': None,
                      'RP': False,
                      'IVA': False}
    if not haplo_args['IVA']:
        assembly_out_path = os.path.join(tmp_dir, 'haplo_out')
        contigs_fasta_path = os.path.join(assembly_out_path, 'contigs.fa')
        haplo_cmd = [HAPLOFLOW,
                      '--read-file', joined_path,
                      '--out', assembly_out_path,
                      '--k', str(haplo_args['kmer']),
                      '--error-rate', str(haplo_args['error']),
                      '--strict', str(haplo_args['strict']),
                      '--filter', str(haplo_args['filter']),
                      '--thres', str(haplo_args['thres']),
                      '--long', str(haplo_args['long'])]
        try:
            run(haplo_cmd, check=True, stdout=PIPE, stderr=STDOUT)
        except CalledProcessError as ex:
            output = ex.output and ex.output.decode('UTF8')
            if output != 'Failed to make first seed. Cannot continue\n':
                logger.warning('Haploflow failed to assemble.', exc_info=True)
                logger.warning(output)
            with open(contigs_fasta_path, 'a'):
                pass
    else:
        assembly_out_path = os.path.join(tmp_dir, 'iva_out')
        contigs_fasta_path = os.path.join(assembly_out_path, 'contigs.fasta')
        iva_args = [IVA, '--fr', joined_path, '-t', '2']
        if merged_contigs_csv is not None:
            seeds_fasta_path = os.path.join(tmp_dir, 'seeds.fasta')
            with open(seeds_fasta_path, 'w') as seeds_fasta:
                SeqIO.write((SeqRecord(Seq(row['contig']), f'seed-{i}', '', '')
                             for i, row in enumerate(DictReader(merged_contigs_csv))),
                            seeds_fasta,
                            'fasta')
                seeds_size = seeds_fasta.tell()
            if seeds_size > 0:
                iva_args.extend(['--contigs', seeds_fasta_path, '--make_new_seeds'])
        iva_args.append(assembly_out_path)
        try:
            run(iva_args, check=True, stdout=PIPE, stderr=STDOUT)
        except CalledProcessError as ex:
            output = ex.output and ex.output.decode('UTF8')
            if output != 'Failed to make first seed. Cannot continue\n':
                logger.warning('iva failed to assemble.', exc_info=True)
                logger.warning(output)
            with open(contigs_fasta_path, 'a'):
                pass

    if haplo_args['RP']:
        contigs_firstpass = os.path.join(assembly_out_path, "contigs_firstpass.csv")
        blast_firstpass = os.path.join(assembly_out_path, "blast_firstpass.csv")
        ref_contigs = os.path.join(assembly_out_path, "ref_contigs.csv")
        noref_contigs = os.path.join(assembly_out_path, "noref_contigs.csv")
        with open(contigs_firstpass, 'w') as contigs_firstpass_csv, \
                open(blast_firstpass, 'w') as blast_firstpass_csv:
            contig_count = write_contig_refs(contigs_fasta_path,
                                         contigs_firstpass_csv,
                                         blast_csv=blast_firstpass_csv)
        with open(contigs_firstpass, 'r') as contigs_firstpass_csv, \
                open(ref_contigs, 'w') as ref_contigs_csv, \
                open(noref_contigs, 'w') as noref_contigs_csv:
            num_noref = separate_contigs(contigs_firstpass_csv, ref_contigs_csv, noref_contigs_csv)
        print(f"Assembled {contig_count} contigs in the first pass, of which {num_noref} did not map to a reference.")
        unmapped1_path = os.path.join(assembly_out_path, 'firstpass_unmapped1.fastq')
        unmapped2_path = os.path.join(assembly_out_path, 'firstpass_unmapped2.fastq')
        remap_path = os.path.join(assembly_out_path, 'firstpass_remap.csv')
        if num_noref:
            with open(remap_path, 'w') as remap_csv, \
                    open(os.path.join(assembly_out_path, 'firstpass_remap_counts.csv'), 'w') as counts_csv, \
                    open(os.path.join(assembly_out_path, 'firstpass_remap_conseq.csv'), 'w') as conseq_csv, \
                    open(unmapped1_path, 'w') as unmapped1, \
                    open(unmapped2_path, 'w') as unmapped2, \
                    open(contigs_firstpass, 'r') as contigs_firstpass_csv:
                map_to_contigs(fastq1_path,
                               fastq2_path,
                               contigs_firstpass_csv,
                               remap_csv,
                               counts_csv,
                               conseq_csv,
                               unmapped1,
                               unmapped2,
                               assembly_out_path, )
            # we want to discard the reads that mapped to the contigs that did not blast to the refs
            ref_reads_path = os.path.join(assembly_out_path, 'ref_reads.fasta')
            noref_reads_path = os.path.join(assembly_out_path, 'noref_reads.fasta')
            with open(remap_path, 'r') as remap_csv, \
                    open(ref_reads_path, 'w') as ref_reads_file, \
                    open(noref_reads_path, 'w') as noref_reads_file, \
                    open(unmapped1_path, 'r') as unmapped1, \
                    open(unmapped2_path, 'r') as unmapped2:
                separate_reads(remap_csv, ref_reads_file, noref_reads_file, unmapped1, unmapped2)
            assembly_out_path = os.path.join(tmp_dir, 'haplo_secondpass_out')
            contigs_fasta_path = os.path.join(assembly_out_path, 'contigs.fa')
            haplo_cmd = [HAPLOFLOW,
                         '--read-file', ref_reads_path,
                         '--out', assembly_out_path,
                         '--k', str(haplo_args['kmer']),
                         '--error-rate', str(haplo_args['error']),
                         '--strict', str(haplo_args['strict']),
                         '--filter', str(haplo_args['filter']),
                         '--thres', str(haplo_args['thres']),
                         '--long', str(haplo_args['long'])]
            run(haplo_cmd, check=True, stdout=PIPE, stderr=STDOUT)

    if haplo_args['merge']:
        fh, abs_path = mkstemp()
        i = 0
        with fdopen(fh, 'w') as new_file:
            with open(contigs_fasta_path) as old_file:
                for line in old_file:
                    if line.startswith('>'):
                        new_file.write(f">contig{i}\n")
                        i += 1
                    else:
                        new_file.write(line)
        copymode(contigs_fasta_path, abs_path)
        remove(contigs_fasta_path)
        move(abs_path, contigs_fasta_path)
        print(f"Number of contigs before trimming and joining: {i}")

        haplo_assembly = Assembly(contigs_file=contigs_fasta_path)
        reads_prefix = os.path.join(tmp_dir, 'reads')
        reads_1 = reads_prefix + '_1.fa'
        reads_2 = reads_prefix + '_2.fa'
        deinterleave(joined_path, reads_1, reads_2, fasta_out=True)
        haplo_assembly._trim_strand_biased_ends(reads_prefix, tag_as_trimmed=True)
        haplo_assembly._remove_contained_contigs(list(haplo_assembly.contigs.keys()))
        haplo_assembly._merge_overlapping_contigs(list(haplo_assembly.contigs.keys()))
        contigs_fasta_path = os.path.join(assembly_out_path, 'contigs_merged.fasta')
        haplo_assembly.write_contigs_to_file(contigs_fasta_path)

    if haplo_args['scaffold']:
        scaffolding_path = os.path.join(assembly_out_path, 'scaffolding')
        scaffold_cmd = ['python3.8',
                        RAGTAG,
                        'scaffold',
                        haplo_args['ref'],
                        contigs_fasta_path,
                        '-o', scaffolding_path,
                        '--aligner', 'nucmer',
                        '--nucmer-params', '--maxmatch -l 30 -c 20']
        run(scaffold_cmd, check=True, stdout=PIPE, stderr=STDOUT)
        new_contigs_fasta_path = os.path.join(scaffolding_path, 'ragtag.scaffold.fasta')
        if os.path.getsize(new_contigs_fasta_path) > 0:
            print('Scaffolding was successful!')
            contigs_fasta_path = new_contigs_fasta_path
        else:
            print('Scaffolding was not successful')

    if haplo_args['patch']:
        patching_path = os.path.join(assembly_out_path, 'patching')
        patch_cmd = ['python3.8',
                     RAGTAG,
                     'patch',
                     contigs_fasta_path,
                     haplo_args['ref'],
                     '-o', patching_path,
                     '--nucmer-params', '--maxmatch -l 30 -c 20']
        run(patch_cmd, check=True, stdout=PIPE, stderr=STDOUT)
        new_contigs_fasta_path = os.path.join(patching_path, 'ragtag.patch.fasta')
        if os.path.getsize(new_contigs_fasta_path) > 0:
            print('Patching was successful!')
            contigs_fasta_path = new_contigs_fasta_path
        else:
            print('Patching was not successful')

    os.chdir(start_dir)
    duration = datetime.now() - start_time
    contig_count = write_contig_refs(contigs_fasta_path,
                                     contigs_csv,
                                     blast_csv=blast_csv)
    logger.info('Assembled %d contigs in %s (%ds) on %s.',
                contig_count,
                duration,
                duration.total_seconds(),
                fastq1_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1')
    parser.add_argument('fastq2')
    parser.add_argument('contigs', type=argparse.FileType('w'))

    args = parser.parse_args()
    denovo(args.fastq1, args.fastq2, args.contigs)
