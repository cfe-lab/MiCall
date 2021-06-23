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
from micall.core.remap import remap, map_to_contigs

IVA = "iva"
DEFAULT_DATABASE = os.path.join(os.path.dirname(__file__),
                                '..',
                                'blast_db',
                                'refs_justHIV.fasta')
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


def run_iva(tmp_dir: str,
            iva_out_path: str,
            reads1: str = None,
            reads2: str = None,
            interleaved: str = None,
            merged_contigs_csv: typing.TextIO = None,
            is_pessimistic: bool = False,
            max_num_contigs: int = 0):
    """ Run IVA with specified arguments.

    :param tmp_dir:             directory for temporary files
    :param iva_out_path:        path for iva ouput files (may not exist yet!)
    :param reads1:              path to forward reads file
    :param reads2:              path to reverse reads file
    :param interleaved:         path to interleaved file
    :param merged_contigs_csv   open file to read contigs that were merged from amplicon reads
    :param is_pessimistic       run IVA in pessimistic mode?
    """
    contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
    if reads1 and reads2:
        iva_args = [IVA, '-f', reads1, '-r', reads2, '-t', '2', '-vv']
    else:
        assert interleaved is not None
        iva_args = [IVA, '--fr', interleaved, '-t', '2', '-vv']
    if is_pessimistic: iva_args.append('--pessimistic')
    if max_num_contigs: iva_args.extend(['--max_contigs', str(max_num_contigs)])
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
    iva_args.append(iva_out_path)
    iva_logfile = os.path.join(tmp_dir, 'iva.log')
    with open(iva_logfile, 'a') as logfile:
        try:
            run(iva_args, check=True, stdout=logfile, stderr=STDOUT)
        except CalledProcessError as ex:
            output = ex.output and ex.output.decode('UTF8')
            if output != 'Failed to make first seed. Cannot continue\n':
                logger.warning('iva failed to assemble.', exc_info=True)
                logger.warning(output)
            with open(contigs_fasta_path, 'a'):
                pass


def separate_contigs(contigs_csv, blast_csv, ref_contigs_csv, noref_contigs_csv):
    """ Separate contigs into those that mapped to or did not map to a reference.

    :param contigs_csv:         file with contigs, open in read mode
    :param blast_csv:           file with blast results, open in read mode (unused so far)
    :param ref_contigs_csv:     file for contigs that mapped to a reference, open in write mode
    :param noref_contigs_csv:   file for contigs that did not map to a reference, open in write mode
    """
    threshold = 0.1
    # is a match threshold sufficient or do we need info from blast_csv as well?
    # we might want to keep all blast matches (regardless of how well they match)
    fieldnames = ['ref', 'match', 'group_ref', 'contig']
    ref_contig_writer = DictWriter(ref_contigs_csv, fieldnames)
    noref_contig_writer = DictWriter(noref_contigs_csv, fieldnames)
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
    return num_match, num_total - num_match


def pess_iva_iterations(tmp_dir, interleaved):
    num_contigs = 0
    num_iterations = 0
    iva_out_path = os.path.join(tmp_dir, 'pessiva_iteration0')
    fieldnames = ['ref', 'match', 'group_ref', 'contig']
    ref_contigs_path = os.path.join(tmp_dir, 'ref_contigs.csv')
    noref_contigs_path = os.path.join(tmp_dir, 'noref_contigs.csv')
    finalcontigs_path = os.path.join(tmp_dir, 'finalcontigs.fasta')
    with open(ref_contigs_path, 'w') as ref_contigs_csv, \
            open(noref_contigs_path, 'w') as noref_contigs_csv:
        ref_contig_writer = DictWriter(ref_contigs_csv, fieldnames)
        ref_contig_writer.writeheader()
        noref_contig_writer = DictWriter(noref_contigs_csv, fieldnames)
        noref_contig_writer.writeheader()
    while True:
        if num_iterations == 0:
            run_iva(tmp_dir, iva_out_path, interleaved=interleaved, is_pessimistic=True, max_num_contigs=1)
        else:
            run_iva(tmp_dir, iva_out_path, reads1=reads1, reads2=reads2, is_pessimistic=True, max_num_contigs=1)
        contigs_dir = os.path.join(iva_out_path, 'pessimistic_contigs.csv')
        blast_dir = os.path.join(iva_out_path, 'pessimistic_blast.csv')
        contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
        with open(contigs_dir, 'w') as pess_contigs_csv, \
                open(blast_dir, 'w') as pess_blast_csv:
            contig_count = write_contig_refs(contigs_fasta_path, pess_contigs_csv, blast_csv=pess_blast_csv)
        if contig_count == 0:
            rmtree(iva_out_path)
            logger.info('Pessimistic IVA finished after %d iterations with %d useful contigs.', num_iterations, num_contigs)
            break
        with open(ref_contigs_path, 'a') as ref_contigs_csv, \
                open(noref_contigs_path, 'a') as noref_contigs_csv, \
                open(contigs_dir, 'r') as pess_contigs_csv, \
                open(blast_dir, 'r') as pess_blast_csv:
            num_match, num_noref = separate_contigs(pess_contigs_csv, pess_blast_csv, ref_contigs_csv, noref_contigs_csv)
        if num_match:
            with open(contigs_fasta_path, 'r') as contigs_fasta, \
                    open(finalcontigs_path, 'a') as finalcontigs:
                finalcontigs.write(contigs_fasta.read())
        num_contigs += num_match # update number of useful contigs
        logger.info('Pessimistic IVA, iteration %d: Assembled %d useful contigs.', num_iterations, num_contigs)
        # we want to use IVA's filtered reads for the next iteration
        reads1 = os.path.join(iva_out_path, 'ivafiltered_1.fa')
        reads2 = os.path.join(iva_out_path, 'ivafiltered_2.fa')
        if num_iterations != 0:
            last_outpath = os.path.join(tmp_dir, f'pessiva_iteration{num_iterations-1}')
            rmtree(last_outpath) # clean up all temp files from the iteration before last
        num_iterations += 1
        iva_out_path = os.path.join(tmp_dir, f'pessiva_iteration{num_iterations}')
    last_outpath = os.path.join(tmp_dir, f'pessiva_iteration{num_iterations - 1}')
    rmtree(last_outpath)  # clean up all temp files from the iteration before last
    return finalcontigs_path


def denovo(fastq1_path: str,
           fastq2_path: str,
           contigs_csv: typing.TextIO,
           work_dir: str = '.',
           merged_contigs_csv: typing.TextIO = None,
           blast_csv: typing.TextIO = None):
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

    if True: # set to pessimistic option!
        logger.info('First pass over reads with pessimistic IVA...')
        contigs_fasta_path = pess_iva_iterations(tmp_dir, joined_path)

    #iva_out_path = os.path.join(tmp_dir, 'iva_out')
    #run_iva(tmp_dir, joined_path, iva_out_path, merged_contigs_csv, is_pessimistic=False)

    #contigs_fasta_path = os.path.join(iva_out_path, 'contigs.fasta')
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
