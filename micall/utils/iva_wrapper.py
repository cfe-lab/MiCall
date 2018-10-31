import multiprocessing
import os
from datetime import datetime, timedelta

import pyfastaq
import iva

DEFAULT_TIMEOUT_SECONDS = 60.0


class IvaError(Exception):
    """ Base class for this module's exceptions. """
    pass


class InputError(IvaError):
    """ Raised when input files can't be read. """
    pass


class SeedingError(IvaError):
    """ Raised when assembly fails to make the first seed. """
    pass


class TimedAssembly(iva.assembly.Assembly):
    def __init__(self, *args, timeout_seconds=DEFAULT_TIMEOUT_SECONDS, **kwargs):
        super().__init__(*args, **kwargs)
        self.timeout_delta = timedelta(seconds=timeout_seconds)
        self.stop_time = None

    def read_pair_extend(self, reads_prefix, out_prefix):
        self.stop_time = datetime.now() + self.timeout_delta
        super().read_pair_extend(reads_prefix, out_prefix)

    def _worth_extending(self):
        if self.stop_time is not None and datetime.now() > self.stop_time:
            return False
        return super()._worth_extending()


def assemble(outdir,
             reads_fwd,
             reads_rev,
             max_insert=800,
             threads=1,
             smalt_k=19,
             smalt_s=11,
             smalt_id=0.5,
             ctg_iter_trim=10,
             ext_min_cov=10,
             ext_min_ratio=4.0,
             verbose=0,
             keep_files=False,
             ext_max_bases=100,
             ext_min_clip=3,
             max_contigs=50,
             seed_start_length=None,
             seed_min_kmer_cov=25,
             seed_max_kmer_cov=1_000_000,
             seed_ext_max_bases=50,
             seed_overlap_length=None,
             seed_ext_min_cov=10,
             seed_ext_min_ratio=4.0,
             strand_bias=0.0,
             timeout_seconds=DEFAULT_TIMEOUT_SECONDS):
    """ Run de novo assembly to build contigs from overlapping reads.

    :param str outdir: Name of output directory (must not already exist)
    :param str reads_fwd: Name of forward reads fasta/q file.
    :param str reads_rev: Name of reverse reads fasta/q file.
    :param int max_insert: Maximum insert size (includes read length). Reads
        with inferred insert size more than the maximum will not be used to
        extend contigs
    :param int threads: Number of threads to run in parallel
    :param int smalt_k: kmer hash length in SMALT (the -k option in smalt
        index)
    :param int smalt_s: kmer hash step size in SMALT (the -s option in smalt
        index)
    :param float smalt_id: Minimum identity threshold for mapping to be
        reported (the -y option in smalt map)
    :param int ctg_iter_trim: During iterative extension, number of bases to
        trim off the end of a contig when extension fails (then try extending
        again)
    :param int ext_min_cov: Minimum kmer depth needed to use that kmer to
        extend a contig
    :param float ext_min_ratio: Sets N, where kmer for extension must be at
        least N times more abundant than next most common kmer
    :param int ext_max_bases: Maximum number of bases to try to extend on each
        iteration
    :param ext_min_clip: Set minimum number of bases soft clipped off a read
        for those bases to be used for extension
    :param max_contigs: Maximum number of contigs allowed in the assembly. No
        more seeds generated if the cutoff is reached
    :param seed_start_length: When making a seed sequence, use the most common
        kmer of this length. Default is to use the minimum of (median read
        length, 95). Warning: it is not recommended to set this higher than 95
    :param int seed_min_kmer_cov: Minimum kmer coverage of initial seed
    :param int seed_max_kmer_cov: Maximum kmer coverage of initial seed
    :param int seed_ext_max_bases: Maximum number of bases to try to extend on
        each iteration
    :param int seed_overlap_length: Number of overlapping bases needed between
        read and seed to use that read to extend [seed_start_length]
    :param seed_ext_min_cov: Minimum kmer depth needed to use that kmer to
        extend a contig
    :param float seed_ext_min_ratio: Sets N, where kmer for extension must be
        at least N times more abundant than next most common kmer
    :param float strand_bias: Set strand bias cutoff of mapped reads when
        trimming contig ends, in the interval [0,0.5]. A value of x means that
        a base needs min(fwd_depth, rev_depth) / total_depth <= x. The only
        time this should be used is with libraries with overlapping reads
        (ie fragment length < 2*read length), and even then, it can make
        results worse. If used, try a low value like 0.1 first.
    :param int verbose: Be verbose by printing messages to stdout. Use up to
        3 for increasing verbosity.
    :param bool keep_files: Keep intermediate files (could be many!).
    :param float timeout_seconds: number of seconds to continue trying to
        assemble contigs. (Not a hard timeout.)
    """
    seed_stop_length = int(0.9 * max_insert)

    if os.path.exists(outdir):
        raise FileExistsError(f"Output directory {outdir} already exists.")

    kmc_threads = threads

    iva.external_progs.get_all_versions(iva.external_progs.assembly_progs)

    os.mkdir(outdir)
    os.chdir(outdir)

    log_file = 'info.txt'
    iva.external_progs.write_prog_info('iva', log_file)

    reads_prefix = 'reads'
    reads_1 = reads_prefix + '_1.fa'
    reads_2 = reads_prefix + '_2.fa'
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0

    reads_for_trimming_1 = reads_fwd
    reads_for_trimming_2 = reads_rev

    fq_to_convert_to_fa_1 = reads_for_trimming_1
    fq_to_convert_to_fa_2 = reads_for_trimming_2

    if threads <= 1:
        pyfastaq.tasks.to_fasta(fq_to_convert_to_fa_1, reads_1, line_length=0)
        pyfastaq.tasks.to_fasta(fq_to_convert_to_fa_2, reads_2, line_length=0)
    else:
        p1 = multiprocessing.Process(target=pyfastaq.tasks.to_fasta,
                                     args=(fq_to_convert_to_fa_1, reads_1),
                                     kwargs={'line_length': 0})
        p2 = multiprocessing.Process(target=pyfastaq.tasks.to_fasta,
                                     args=(fq_to_convert_to_fa_2, reads_2),
                                     kwargs={'line_length': 0})
        p1.start()
        p2.start()
        p1.join()
        p2.join()

        if p1.exitcode != 0 or p2.exitcode != 0:
            raise InputError('Error in input reads files. Cannot continue.')

    pyfastaq.sequences.Fasta.line_length = original_line_length

    contigs = None
    make_new_seeds = True

    assembly = TimedAssembly(
        contigs,
        verbose=verbose,
        clean=not keep_files,
        map_index_k=smalt_k,
        map_index_s=smalt_s,
        threads=threads,
        kmc_threads=kmc_threads,
        map_minid=smalt_id,
        contig_iter_trim=ctg_iter_trim,
        ext_min_cov=ext_min_cov,
        ext_min_ratio=ext_min_ratio,
        ext_bases=ext_max_bases,
        min_clip=ext_min_clip,
        max_contigs=max_contigs,
        make_new_seeds=make_new_seeds,
        seed_start_length=seed_start_length,
        seed_stop_length=seed_stop_length,
        seed_min_kmer_count=seed_min_kmer_cov,
        seed_max_kmer_count=seed_max_kmer_cov,
        seed_ext_max_bases=seed_ext_max_bases,
        seed_overlap_length=seed_overlap_length,
        seed_min_cov=seed_ext_min_cov,
        seed_min_ratio=seed_ext_min_ratio,
        max_insert=max_insert,
        strand_bias=strand_bias,
        timeout_seconds=timeout_seconds)

    seed_name = assembly.add_new_seed_contig(reads_1, reads_2)
    if seed_name is None:
        error_message = 'Failed to make first seed. Cannot continue'
        with open(log_file, 'a') as f:
            print(error_message, file=f)
        raise SeedingError(error_message)

    assembly.read_pair_extend(reads_prefix, 'iteration')

    final_contigs = 'contigs.fasta'

    assembly.write_contigs_to_file(final_contigs, min_length=100, order_by_orfs=True, prefix='contig')

    if not keep_files:
        os.unlink(reads_1)
        os.unlink(reads_2)
