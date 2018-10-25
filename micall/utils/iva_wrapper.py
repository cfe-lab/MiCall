#!/mnt/data/don/git/MiCall/venv_micall/bin/python3
# Copyright (c) 2014-2016 Genome Research Ltd.
#
# This file is part of IVA.
#
# IVA is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

import argparse
import os
import sys
import multiprocessing
import pyfastaq
import iva


parser = argparse.ArgumentParser(
    usage = '%(prog)s [options] {-f reads_fwd -r reads_rev | --fr reads} <output directory>')

parser.add_argument('outdir', help='Name of output directory (must not already exist)', metavar='Output directory')


io_group = parser.add_argument_group('Input and output')
io_group.add_argument('-f', '--reads_fwd', action=iva.common.abspathAction, help='Name of forward reads fasta/q file. Must be used in conjunction with --reads_rev', metavar='filename[.gz]')
io_group.add_argument('-r', '--reads_rev', action=iva.common.abspathAction, help='Name of reverse reads fasta/q file. Must be used in conjunction with --reads_fwd', metavar='filename[.gz]')
io_group.add_argument('--fr', action=iva.common.abspathAction, dest='reads', help='Name of interleaved fasta/q file', metavar='filename[.gz]')
io_group.add_argument('--keep_files', action='store_true', help='Keep intermediate files (could be many!). Default is to delete all unnecessary files')
io_group.add_argument('--contigs', action=iva.common.abspathAction, help='Fasta file of contigs to be extended. Incompatible with --reference', metavar='filename[.gz]')
io_group.add_argument('--reference', action=iva.common.abspathAction, help='EXPERIMENTAL! This option is EXPERIMENTAL, not recommended, and has not been tested! Fasta file of reference genome, or parts thereof. IVA will try to assemble one contig per sequence in this file. Incompatible with --contigs', metavar='filename[.gz]')
io_group.add_argument('-v', '--verbose', action='count', help='Be verbose by printing messages to stdout. Use up to three times for increasing verbosity.', default=0)


mapping_group = parser.add_argument_group('SMALT mapping options')
mapping_group.add_argument('-k', '--smalt_k', type=int, help='kmer hash length in SMALT (the -k option in smalt index) [%(default)s]', default=19, metavar='INT')
mapping_group.add_argument('-s', '--smalt_s', type=int, help='kmer hash step size in SMALT (the -s option in smalt index) [%(default)s]', default=11, metavar='INT')
mapping_group.add_argument('-y', '--smalt_id', type=float, help='Minimum identity threshold for mapping to be reported (the -y option in smalt map) [%(default)s]', default=0.5, metavar='FLOAT')


contig_group = parser.add_argument_group('Contig options')
contig_group.add_argument('--ctg_first_trim', type=int, help='Number of bases to trim off the end of every contig before extending for the first time [%(default)s]', default=25, metavar='INT')
contig_group.add_argument('--ctg_iter_trim', type=int, help='During iterative extension, number of bases to trim off the end of a contig when extension fails (then try extending again) [%(default)s]', default=10, metavar='INT')
contig_group.add_argument('--ext_min_cov', type=int, help='Minimum kmer depth needed to use that kmer to extend a contig [%(default)s]', default=10, metavar='INT')
contig_group.add_argument('--ext_min_ratio', type=float, help='Sets N, where kmer for extension must be at least N times more abundant than next most common kmer [%(default)s]', default=4, metavar='FLOAT')
contig_group.add_argument('--ext_max_bases', type=int, help='Maximum number of bases to try to extend on each iteration [%(default)s]', default=100, metavar='INT')
contig_group.add_argument('--ext_min_clip', type=int, help='Set minimum number of bases soft clipped off a read for those bases to be used for extension [%(default)s]', default=3, metavar='INT')
contig_group.add_argument('--max_contigs', type=int, help='Maximum number of contigs allowed in the assembly. No more seeds generated if the cutoff is reached [%(default)s]', metavar='INT', default=50)


seed_group = parser.add_argument_group('Seed generation options')
seed_group.add_argument('--make_new_seeds', action='store_true', help='When no more contigs can be extended, generate a new seed. This is forced to be true when --contigs is not used')
seed_group.add_argument('--seed_start_length', type=int, help='When making a seed sequence, use the most common kmer of this length. Default is to use the minimum of (median read length, 95). Warning: it is not recommended to set this higher than 95', metavar='INT', default=None)
seed_group.add_argument('--seed_stop_length', type=int, help='Stop extending seed using perfect matches from reads when this length is reached. Future extensions are then made by treating the seed as a contig [0.9*max_insert]', default=0, metavar='INT')
seed_group.add_argument('--seed_min_kmer_cov', type=int, help='Minimum kmer coverage of initial seed [%(default)s]', default=25, metavar='INT')
seed_group.add_argument('--seed_max_kmer_cov', type=int, help='Maximum kmer coverage of initial seed [%(default)s]', default=1000000, metavar='INT')
seed_group.add_argument('--seed_ext_max_bases', type=int, help='Maximum number of bases to try to extend on each iteration [%(default)s]', default=50, metavar='INT')
seed_group.add_argument('--seed_overlap_length', type=int, help='Number of overlapping bases needed between read and seed to use that read to extend [seed_start_length]', metavar='INT')
seed_group.add_argument('--seed_ext_min_cov', type=int, help='Minimum kmer depth needed to use that kmer to extend a contig [%(default)s]', default=10, metavar='INT')
seed_group.add_argument('--seed_ext_min_ratio', type=float, help='Sets N, where kmer for extension must be at least N times more abundant than next most common kmer [%(default)s]', default=4, metavar='FLOAT')


trimming_group = parser.add_argument_group('Read trimming options')
trimming_group.add_argument('--trimmomatic', action=iva.common.abspathAction, help='Provide location of trimmomatic.jar file to enable read trimming. Required if --adapters used', metavar='FILENAME')
trimming_group.add_argument('--trimmo_qual', help='Trimmomatic options used to quality trim reads [%(default)s]', default='LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20', metavar='STRING')
trimming_group.add_argument('--adapters', action=iva.common.abspathAction, help='Fasta file of adapter sequences to be trimmed off reads. If used, must also use --trimmomatic. Default is file of adapters supplied with IVA', metavar='FILENAME')
trimming_group.add_argument('--min_trimmed_length', type=int, help='Minimum length of read after trimming [%(default)s]', default=50, metavar='INT')
trimming_group.add_argument('--pcr_primers', action=iva.common.abspathAction, help='FASTA file of primers. The first perfect match found to a sequence in the primers file will be trimmed off the start of each read. This is run after trimmomatic (if --trimmomatic used)', metavar='FILENAME')


other_group = parser.add_argument_group('Other options')
other_group.add_argument('-i', '--max_insert', type=int, help='Maximum insert size (includes read length). Reads with inferred insert size more than the maximum will not be used to extend contigs [%(default)s]', default=800, metavar='INT')
other_group.add_argument('-t', '--threads', type=int, help='Number of threads to use [%(default)s]', default=1, metavar='INT')
other_group.add_argument('--kmc_onethread', action='store_true', help='Force kmc to use one thread. By default the value of -t/--threads is used when running kmc')
other_group.add_argument('--strand_bias', type=float, help='Set strand bias cutoff of mapped reads when trimming contig ends, in the interval [0,0.5]. A value of x means that a base needs min(fwd_depth, rev_depth) / total_depth <= x. The only time this should be used is with libraries with overlapping reads (ie fragment length < 2*read length), and even then, it can make results worse. If used, try a low value like 0.1 first [%(default)s]', default=0, metavar='FLOAT in [0,0.5]')
other_group.add_argument('--test', action='store_true', help='Run using built in test data. All other options will be ignored, except the mandatory output directory, and --trimmomatic and --threads can be also be used')
other_group.add_argument('--version', action='version', version=iva.common.version)

options = parser.parse_args()

if options.test:
    print('Running iva in test mode...')
    this_script = os.path.abspath(__file__)
    tester = iva.test_data_runner.Tester(options.outdir, this_script, trimmo_jar=options.trimmomatic, threads=options.threads)
    tester.run()
    sys.exit()


if options.seed_stop_length == 0:
    options.seed_stop_length = int(0.9 * options.max_insert)

if not (0 <= options.strand_bias <= 0.5):
    print('Error! strand bias must in the interval [0, 0.5]. Cannot continue because it\'s', options.strand_bias, file=sys.stderr)
    sys.exit(1)


if options.adapters and not options.trimmomatic:
    print('Error! --adapters used, but not --trimmomatic. I need the trimmomatic jar file. Cannot contiue', file=sys.stderr)
    sys.exit(1)

if not (bool(options.reads) ^ bool(options.reads_fwd and options.reads_rev)) or (bool(options.reads_fwd) != bool(options.reads_rev)):
    print('Error! Must use options: -f/--reads_fwd and -r/--reads_rev together, or just use --fr on its own. Cannot continue', file=sys.stderr)
    sys.exit(1)

if options.contigs and options.reference:
    print('Error! Cannot use both of --contgs and --reference. Cannot continue', file=sys.stderr)
    sys.exit(1)

if options.reference:
    print('WARNING. The option --reference has been used. It is EXPERIMENTAL and it is probably better to not use it!', file=sys.stderr)

if os.path.exists(options.outdir):
    print('Error! Output directory', options.outdir, 'already exists. Cannot continue', file=sys.stderr)
    sys.exit(1)


if options.kmc_onethread:
    kmc_threads = 1
else:
    kmc_threads = options.threads


iva.external_progs.get_all_versions(iva.external_progs.assembly_progs)

try:
    os.mkdir(options.outdir)
except:
    print('Error making output directory', options.outdir)
    sys.exit(1)

os.chdir(options.outdir)

log_file = 'info.txt'
iva.external_progs.write_prog_info('iva', log_file)

reads_prefix = 'reads'
reads_1 = reads_prefix + '_1.fa'
reads_2 = reads_prefix + '_2.fa'
original_line_length = pyfastaq.sequences.Fasta.line_length
pyfastaq.sequences.Fasta.line_length = 0

if options.reads and not options.trimmomatic:
    pyfastaq.tasks.deinterleave(options.reads, reads_1, reads_2, fasta_out=True)
else:
    to_delete = []

    if options.reads:
        reads_for_trimming_1 = 'reads.untrimmed_1.fq'
        reads_for_trimming_2 = 'reads.untrimmed_2.fq'
        pyfastaq.tasks.deinterleave(options.reads, reads_for_trimming_1, reads_for_trimming_2)
        to_delete.append(reads_for_trimming_1)
        to_delete.append(reads_for_trimming_2)
    else:
        reads_for_trimming_1 = options.reads_fwd
        reads_for_trimming_2 = options.reads_rev

    if options.trimmomatic:
        trimmed_reads_prefix = 'reads.trimmed'
        if options.adapters is None:
            extractor = iva.egg_extract.Extractor(os.path.abspath(os.path.join(os.path.dirname(iva.__file__), os.pardir)))
            egg_adapters = os.path.join('iva', 'read_trim', 'adapters.fasta')
            options.adapters = 'adapters.fasta'
            extractor.copy_file(egg_adapters, options.adapters)

        assert os.path.exists(options.adapters)

        iva.read_trim.run_trimmomatic(reads_for_trimming_1, reads_for_trimming_2, trimmed_reads_prefix, options.trimmomatic, options.adapters, minlen=options.min_trimmed_length, verbose=options.verbose, threads=options.threads, qual_trim=options.trimmo_qual)
        fq_to_convert_to_fa_1 = trimmed_reads_prefix + '_1.fq'
        fq_to_convert_to_fa_2 = trimmed_reads_prefix + '_2.fq'
        to_delete.append(fq_to_convert_to_fa_1)
        to_delete.append(fq_to_convert_to_fa_2)
    else:
        fq_to_convert_to_fa_1 = reads_for_trimming_1
        fq_to_convert_to_fa_2 = reads_for_trimming_2

    p1 = multiprocessing.Process(target=pyfastaq.tasks.to_fasta, args=(fq_to_convert_to_fa_1, reads_1), kwargs={'line_length':0})
    p2 = multiprocessing.Process(target=pyfastaq.tasks.to_fasta, args=(fq_to_convert_to_fa_2, reads_2), kwargs={'line_length':0})
    p1.start()

    if options.threads == 1:
        p1.join()

    p2.start()
    p2.join()

    if options.threads > 1:
        p1.join()

    if p1.exitcode != 0 or p2.exitcode != 0:
        print('Error in input reads files. Cannot continue', file=sys.stderr)
        sys.exit(1)

    for fname in to_delete:
        os.unlink(fname)

if options.pcr_primers:
    tmp_reads_1 = 'reads_1.pcr_trim.fa'
    tmp_reads_2 = 'reads_2.pcr_trim.fa'
    pyfastaq.tasks.sequence_trim(reads_1, reads_2, tmp_reads_1, tmp_reads_2, options.pcr_primers, min_length=options.min_trimmed_length, check_revcomp=True)
    os.rename(tmp_reads_1, reads_1)
    os.rename(tmp_reads_2, reads_2)

pyfastaq.sequences.Fasta.line_length = original_line_length

if options.contigs:
    contigs = 'contigs_to_extend.fasta'
    pyfastaq.tasks.to_fasta(options.contigs, contigs, line_length=60, strip_after_first_whitespace=True)
elif options.reference:
    print('WARNING. The option --reference has been used. Trying to use reference file to generate starting contig.', file=sys.stderr)
    print(' ... if this throws errors, then try running without the --reference option', file=sys.stderr)
    reference = 'reference_in.fasta'
    pyfastaq.tasks.to_fasta(options.reference, reference, line_length=60, strip_after_first_whitespace=True)
    p = iva.seed_processor.SeedProcessor(
        reference,
        reads_1,
        reads_2,
        'seeds.fasta',
        index_k = options.smalt_k,
        index_s = options.smalt_s,
        threads = options.threads,
        kmc_threads = kmc_threads,
        max_insert = options.max_insert,
        minid = 0.9,
        seed_stop_length = options.seed_stop_length,
        extend_length = options.seed_ext_max_bases,
        overlap_length = options.seed_overlap_length,
        ext_min_cov = options.seed_ext_min_cov,
        ext_min_ratio = options.seed_ext_min_ratio,
        verbose = options.verbose,
        seed_length = options.seed_start_length,
        seed_min_count = options.seed_min_kmer_cov,
        seed_max_count = options.seed_max_kmer_cov
    )
    p.process()
    contigs = 'seeds.fasta'
else:
    contigs = None
    options.make_new_seed = True

assembly = iva.assembly.Assembly(
    contigs,
    verbose = options.verbose,
    clean = not options.keep_files,
    map_index_k = options.smalt_k,
    map_index_s = options.smalt_s,
    threads = options.threads,
    kmc_threads = kmc_threads,
    map_minid = options.smalt_id,
    contig_iter_trim = options.ctg_iter_trim,
    ext_min_cov = options.ext_min_cov,
    ext_min_ratio = options.ext_min_ratio,
    ext_bases = options.ext_max_bases,
    min_clip = options.ext_min_clip,
    max_contigs = options.max_contigs,
    make_new_seeds = options.make_new_seeds,
    seed_start_length = options.seed_start_length,
    seed_stop_length = options.seed_stop_length,
    seed_min_kmer_count = options.seed_min_kmer_cov,
    seed_max_kmer_count = options.seed_max_kmer_cov,
    seed_ext_max_bases = options.seed_ext_max_bases,
    seed_overlap_length = options.seed_overlap_length,
    seed_min_cov = options.seed_ext_min_cov,
    seed_min_ratio = options.seed_ext_min_ratio,
    max_insert = options.max_insert,
    strand_bias = options.strand_bias
)

if options.contigs:
    assembly.trim_contigs(options.ctg_first_trim)
elif not options.reference:
    seed_name = assembly.add_new_seed_contig(reads_1, reads_2)
    if seed_name is None:
        error_message = 'Failed to make first seed. Cannot continue'
        print(error_message, file=sys.stderr)
        with open(log_file, 'a') as f:
            print(error_message, file=f)
        f.close()
        sys.exit(1)

assembly.read_pair_extend(reads_prefix, 'iteration')

final_contigs = 'contigs.fasta'

if options.trimmomatic or options.pcr_primers:
    pre_trim_contigs = 'contigs.pre_trim.fasta'
    assembly.write_contigs_to_file(pre_trim_contigs, min_length=100, order_by_orfs=True, prefix='contig')
    iva.contig_trim.trim_primers_and_adapters(pre_trim_contigs, final_contigs, options.adapters, options.pcr_primers, min_length=100)
    if not options.keep_files:
        os.unlink(pre_trim_contigs)
        os.unlink(pre_trim_contigs + '.fai')
else:
    assembly.write_contigs_to_file(final_contigs, min_length=100, order_by_orfs=True, prefix='contig')

if not options.keep_files:
    os.unlink(reads_1)
    os.unlink(reads_2)
    if contigs is not None:
        os.unlink(contigs)

