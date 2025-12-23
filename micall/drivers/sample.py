""" Data needed to process a sample """
import logging
import os

import typing
from csv import DictReader
from pathlib import Path

from micall.core.aln2counts import aln2counts
from micall.core.amplicon_finder import write_merge_lengths_plot, merge_for_entropy
from micall.core.cascade_report import CascadeReport
from micall.core.coverage_plots import coverage_plot, concordance_plot
from micall.core.plot_contigs import plot_genome_coverage
from micall.core.prelim_map import prelim_map
from micall.core.project_config import ProjectConfig
from micall.core.remap import remap, map_to_contigs
from micall.core.sam2aln import sam2aln
from micall.core.trim_fastqs import trim
from micall.core.denovo import denovo
from micall.g2p.fastq_g2p import fastq_g2p, DEFAULT_MIN_COUNT, MIN_VALID, MIN_VALID_PERCENT
from micall.utils.driver_utils import makedirs
from micall.utils.fasta_to_csv import fasta_to_csv
from micall.utils.csv_to_fasta import csv_to_fasta, NoContigsInCSV
from micall.utils.referencefull_contig_stitcher import referencefull_contig_stitcher
from micall.utils.cat import cat as concatenate_files
from micall.utils.work_dir import WorkDir
from micall.utils.exact_coverage import calculate_exact_coverage_from_csv, write_coverage_csv
from contextlib import contextmanager

logger = logging.getLogger(__name__)


def prepend_prefix_to_basename(prefix: str, path: str):
    dir_name, base_name = os.path.split(path)
    return os.path.join(dir_name, prefix + base_name)


@contextmanager
def open_files(**files):
    """ Context manager that will open files and close them at the end.

    @param files: keyworded files to be opened. The keyword corresponds to the key of the output dictionary of opened
        file handles. The input arguments should be tuples of the path of the file to be opened, and the mode to open
        it: (path, mode).

    Yields: A dictionary of opened files, {keyword: file handle}
    """
    files_opened: typing.Dict[str, typing.Optional[typing.IO]] = {}
    try:
        for name, file_info in files.items():
            if file_info is None:
                files_opened[name] = None
            else:
                file = file_info[0]
                mode = file_info[1]
                files_opened[name] = open(file, mode)
            # Could add support for open file handles in the future by checking if file_info is an instance of TextIO
        yield files_opened
    finally:
        filenames = []
        for filename, file in files_opened.items():
            if file is None:
                pass
            else:
                try:
                    file.close()
                except (IOError, OSError):
                    filenames.append(filename)
                    pass
        if len(filenames) > 0:
            logger.error(f"The following files could not be closed: {filenames}")
            raise IOError


def exclude_extra_seeds(excluded_seeds: typing.Iterable[str],
                        project_code: typing.Optional[str] = None) -> typing.Sequence[str]:
    if project_code == 'HIVGHA':
        return tuple(excluded_seeds)
    projects = ProjectConfig.loadDefault()
    hivgha_seeds = projects.getProjectSeeds('HIVGHA')
    extra_exclusions = {seed
                        for seed in hivgha_seeds
                        # Exclude any circulating recombinant forms (CRF).
                        if seed.startswith('HIV1-CRF')}
    return sorted(extra_exclusions | set(excluded_seeds))


class Sample:
    def __init__(self,
                 basespace_id=None,
                 basespace_href=None,
                 rank=None,
                 debug_remap=False,
                 scratch_path=None,
                 skip: typing.Iterable[str] = (),
                 **paths):
        """ Record the details.

        :param str fastq1: path to first FASTQ file
        :param str fastq2: path to second FASTQ file
        :param str basespace_id: id from BaseSpace
        :param str basespace_href: web address within BaseSpace API
        :param str rank: rank within the run, like '2 of 10'
        :param bool debug_remap: write debug files during remap step
        :param str scratch_path: path where unknown outputs are written
        :param skip: (step_name,) of steps to skip
        :param paths: {output_name: file_path}
        """
        fastq1 = paths.get('fastq1')
        if 'fastq2' in paths:
            pass
        elif fastq1:
            if '_R1_' not in fastq1:
                raise ValueError(
                    "fastq2 not given, and fastq1 does not contain '_R1_'.")
            paths['fastq2'] = fastq1.replace('_R1_', '_R2_')
        if fastq1:
            self.name: typing.Optional[str] = '_'.join(os.path.basename(fastq1).split('_')[:2])
        else:
            self.name = None
        self.basespace_id = basespace_id
        self.basespace_href = basespace_href
        self.rank = rank
        self.debug_remap = debug_remap
        self.scratch_path = scratch_path
        self.skip = skip
        self.paths = paths
        self.project_code = None

    def __repr__(self):
        fastq1 = self.paths.get('fastq1')
        if fastq1 is None:
            return 'Sample()'
        return 'Sample(fastq1={!r})'.format(fastq1)

    def __str__(self):
        result = 'Sample'
        if self.name is not None:
            result += ' '
            result += self.name
        if self.rank is not None:
            result += ' ({})'.format(self.rank)
        return result

    def __getattr__(self, output_name):
        if output_name.startswith('__'):
            raise AttributeError(output_name)
        output_path = self.paths.get(output_name)
        if output_path is None:
            output_path = self.get_default_path(output_name)
        return output_path

    def get_default_path(self, output_name):
        if self.scratch_path is None:
            raise AttributeError(
                'Unknown output {} and no scratch path.'.format(output_name))
        for extension in ('csv', 'fastq', 'pdf', 'svg', 'png', 'fasta'):
            if output_name.endswith('_'+extension):
                file_name = output_name[:-(len(extension)+1)] + '.' + extension
                break
        else:
            file_name = output_name
        return os.path.join(self.scratch_path, file_name)

    def get_scratch_path(self):
        if self.scratch_path is not None:
            return self.scratch_path
        return os.path.dirname(self.cascade_csv)

    def process(self,
                pssm,
                excluded_seeds: typing.Iterable[str] = (),
                excluded_projects: typing.Iterable[str] = (),
                force_gzip=False,
                use_denovo=False):
        """ Process a single sample.

        :param pssm: the pssm library for running G2P analysis
        :param excluded_seeds: seeds to exclude from mapping
        :param excluded_projects: project codes to exclude from reporting
        :param bool force_gzip: treat FASTQ files as gzipped, even when they
            don't end in .gz
        :param bool use_denovo: True if de novo assembly should be used,
            instead of bowtie2 mapping against references.
        """
        logger.info('Processing %s (%r).', self, self.fastq1)
        scratch_path = self.get_scratch_path()
        makedirs(scratch_path)
        use_gzip = force_gzip or self.fastq1.endswith('.gz')

        sample_info = self.load_sample_info()
        excluded_seeds = exclude_extra_seeds(excluded_seeds,
                                             sample_info.get('project'))

        with open(self.read_summary_csv, 'w') as read_summary:
            trim((self.fastq1, self.fastq2),
                 self.bad_cycles_csv,
                 (self.trimmed1_fastq, self.trimmed2_fastq),
                 summary_file=read_summary,
                 use_gzip=use_gzip,
                 skip=self.skip,
                 project_code=sample_info.get('project'))

        if use_denovo:
            logger.info('Running merge_for_entropy on %s.', self)
            with open(self.read_entropy_csv, 'w') as read_entropy_csv:
                merge_for_entropy(self.trimmed1_fastq,
                                  self.trimmed2_fastq,
                                  read_entropy_csv,
                                  scratch_path)

            write_merge_lengths_plot(self.read_entropy_csv,
                                     self.merge_lengths_svg)

        logger.info('Running fastq_g2p on %s.', self)
        with open(self.trimmed1_fastq) as fastq1, \
                open(self.trimmed2_fastq) as fastq2, \
                open(self.g2p_csv, 'w') as g2p_csv, \
                open(self.g2p_summary_csv, 'w') as g2p_summary_csv, \
                open(self.g2p_unmapped1_fastq, 'w') as g2p_unmapped1, \
                open(self.g2p_unmapped2_fastq, 'w') as g2p_unmapped2, \
                open(self.g2p_aligned_csv, 'w') as g2p_aligned_csv, \
                open(self.merged_contigs_csv, 'w') as merged_contigs_csv:

            fastq_g2p(pssm=pssm,
                      fastq1=fastq1,
                      fastq2=fastq2,
                      g2p_csv=g2p_csv,
                      g2p_summary_csv=g2p_summary_csv,
                      unmapped1=g2p_unmapped1,
                      unmapped2=g2p_unmapped2,
                      aligned_csv=g2p_aligned_csv,
                      min_count=DEFAULT_MIN_COUNT,
                      min_valid=MIN_VALID,
                      min_valid_percent=MIN_VALID_PERCENT,
                      merged_contigs_csv=merged_contigs_csv)

        if use_denovo:
            self.run_denovo(excluded_seeds)
        else:
            self.run_mapping(excluded_seeds)

        if use_denovo:
            # Run exact coverage after remap_conseq.csv has been generated
            logger.info('Running exact_coverage on %s.', self)
            with open(self.remap_csv, 'r') as aligned_csv, \
                 open(self.remap_conseq_csv, 'r') as remap_conseq_file, \
                 open(self.exact_coverage_csv, 'w') as exact_coverage_csv:
                coverage, contigs = calculate_exact_coverage_from_csv(
                    aligned_csv,
                    remap_conseq_file,
                    overlap_size=70)
                write_coverage_csv(coverage, contigs, exact_coverage_csv)

        self.process_post_assembly(prefix="",
                                   use_denovo=use_denovo,
                                   excluded_projects=excluded_projects)

        if use_denovo:
            self.process_post_assembly(prefix="unstitched_",
                                       use_denovo=use_denovo,
                                       excluded_projects=excluded_projects)

        logger.info('Finished sample %s.', self)

    def process_post_assembly(self,
                              use_denovo: bool,
                              excluded_projects: typing.Iterable[str],
                              prefix: str,
                              ):

        def with_prefix(path):
            return prepend_prefix_to_basename(prefix, path)

        logger.info('Running sam2aln on %s.', self)
        with open(with_prefix(self.remap_csv)) as remap_csv, \
                open(with_prefix(self.aligned_csv), 'w') as aligned_csv, \
                open(with_prefix(self.conseq_ins_csv), 'w') as conseq_ins_csv, \
                open(with_prefix(self.failed_csv), 'w') as failed_csv, \
                open(with_prefix(self.clipping_csv), 'w') as clipping_csv:

            sam2aln(remap_csv,
                    aligned_csv,
                    conseq_ins_csv,
                    failed_csv,
                    clipping_csv=clipping_csv)

        logger.info('Running aln2counts on %s.', self)
        with open_files(aligned_csv=(with_prefix(self.aligned_csv), 'r'),

                        # Does not need a prefix because it is produced before the denovo/remap split.
                        g2p_aligned_csv=(self.g2p_aligned_csv, 'r'),

                        clipping_csv=(with_prefix(self.clipping_csv), 'r'),
                        nuc_csv=(with_prefix(self.nuc_csv), 'w'),
                        conseq_ins_csv=(with_prefix(self.conseq_ins_csv), 'r'),
                        remap_conseq_csv=(with_prefix(self.remap_conseq_csv), 'r'),
                        contigs_csv=(with_prefix(self.contigs_csv), 'r') if use_denovo else None,
                        exact_coverage_csv=(self.exact_coverage_csv, 'r') if use_denovo and prefix == "" else None,
                        nuc_detail_csv=(with_prefix(self.nuc_details_csv), 'w') if use_denovo else None,
                        amino_csv=(with_prefix(self.amino_csv), 'w'),
                        amino_detail_csv=(with_prefix(self.amino_details_csv), 'w') if use_denovo else None,
                        insertions_csv=(with_prefix(self.insertions_csv), 'w'),
                        conseq_csv=(with_prefix(self.conseq_csv), 'w'),
                        conseq_region_csv=(with_prefix(self.conseq_region_csv), 'w') if use_denovo else None,
                        failed_align_csv=(with_prefix(self.failed_align_csv), 'w'),
                        coverage_summary_csv=(with_prefix(self.coverage_summary_csv), 'w'),
                        genome_coverage_csv=(with_prefix(self.genome_coverage_csv), 'w'),
                        conseq_all_csv=(with_prefix(self.conseq_all_csv), 'w'),
                        conseq_stitched_csv=(with_prefix(self.conseq_stitched_csv), 'w') if use_denovo else None,
                        minimap_hits_csv=(with_prefix(self.minimap_hits_csv), 'w'),
                        alignments_csv=(with_prefix(self.alignments_csv), 'w'),
                        alignments_unmerged_csv=(with_prefix(self.alignments_unmerged_csv), 'w'),
                        alignments_intermediate_csv=(with_prefix(self.alignments_intermediate_csv), 'w'),
                        alignments_overall_csv=(with_prefix(self.alignments_overall_csv), 'w'),
                        concordance_csv=(with_prefix(self.concordance_csv), 'w'),
                        concordance_detailed_csv=(with_prefix(self.concordance_detailed_csv), 'w'),
                        concordance_seed_csv=(with_prefix(self.concordance_seed_csv), 'w')) as opened_files:

            aln2counts(opened_files['aligned_csv'],
                       opened_files['nuc_csv'],
                       opened_files['amino_csv'],
                       opened_files['insertions_csv'],
                       opened_files['conseq_csv'],
                       opened_files['failed_align_csv'],
                       coverage_summary_csv=opened_files['coverage_summary_csv'],
                       clipping_csv=opened_files['clipping_csv'],
                       conseq_ins_csv=opened_files['conseq_ins_csv'],
                       g2p_aligned_csv=opened_files['g2p_aligned_csv'],
                       remap_conseq_csv=opened_files['remap_conseq_csv'],
                       conseq_region_csv=opened_files['conseq_region_csv'],
                       amino_detail_csv=opened_files['amino_detail_csv'],
                       nuc_detail_csv=opened_files['nuc_detail_csv'],
                       genome_coverage_csv=opened_files['genome_coverage_csv'],
                       contigs_csv=opened_files['contigs_csv'],
                       exact_coverage_csv=opened_files['exact_coverage_csv'],
                       conseq_all_csv=opened_files['conseq_all_csv'],
                       conseq_stitched_csv=opened_files['conseq_stitched_csv'],
                       minimap_hits_csv=opened_files['minimap_hits_csv'],
                       alignments_csv=opened_files['alignments_csv'],
                       alignments_unmerged_csv=opened_files['alignments_unmerged_csv'],
                       alignments_intermediate_csv=opened_files['alignments_intermediate_csv'],
                       alignments_overall_csv=opened_files['alignments_overall_csv'],
                       concordance_csv=opened_files['concordance_csv'],
                       concordance_detailed_csv=opened_files['concordance_detailed_csv'],
                       concordance_seed_csv=opened_files['concordance_seed_csv'])

        logger.info('Running coverage_plots on %s.', self)
        os.makedirs(with_prefix(self.coverage_maps))
        with open(with_prefix(self.amino_csv)) as amino_csv, \
             open(with_prefix(self.coverage_scores_csv), 'w') as coverage_scores_csv:
            coverage_plot(amino_csv,
                          coverage_scores_csv,
                          coverage_maps_path=self.coverage_maps,
                          coverage_maps_prefix=self.name,
                          excluded_projects=excluded_projects)

        with open(with_prefix(self.genome_coverage_csv)) as genome_coverage_csv, \
             open(with_prefix(self.minimap_hits_csv)) as minimap_hits_csv:
            plot_genome_coverage(genome_coverage_csv,
                                 minimap_hits_csv if use_denovo else None,
                                 with_prefix(self.genome_coverage_svg))

        with open(with_prefix(self.genome_coverage_csv)) as genome_coverage_csv, \
             open(with_prefix(self.minimap_hits_csv)) as minimap_hits_csv:
            plot_genome_coverage(genome_coverage_csv,
                                 minimap_hits_csv if use_denovo else None,
                                 with_prefix(self.genome_concordance_svg),
                                 use_concordance=True)

        with open(with_prefix(self.concordance_detailed_csv)) as concordance_detailed_csv:
            concordance_plot(concordance_detailed_csv,
                             plot_path=with_prefix(self.coverage_maps),
                             concordance_prefix=self.name)

        logger.info('Running cascade_report on %s.', self)
        with open(self.g2p_summary_csv) as g2p_summary_csv, \
             open(with_prefix(self.remap_counts_csv)) as remap_counts_csv, \
             open(with_prefix(self.aligned_csv)) as aligned_csv, \
             open(with_prefix(self.cascade_csv), 'w') as cascade_csv:
            cascade_report = CascadeReport(cascade_csv)
            cascade_report.g2p_summary_csv = g2p_summary_csv
            cascade_report.remap_counts_csv = remap_counts_csv
            cascade_report.aligned_csv = aligned_csv
            cascade_report.generate()

    def load_sample_info(self) -> dict[str, str]:
        path = Path(self.sample_info_csv)
        if not path.exists():
            sample_info: dict[str, str] = {}
            if self.project_code is not None:
                sample_info['project'] = self.project_code
            return sample_info
        with path.open() as info_file:
            reader = DictReader(info_file)
            return next(reader)

    def run_mapping(self, excluded_seeds):
        logger.info('Running prelim_map on %s.', self)
        scratch_path = self.get_scratch_path()
        with WorkDir.using(Path(scratch_path)):
            prelim_map(Path(self.trimmed1_fastq),
                       Path(self.trimmed2_fastq),
                       Path(self.prelim_csv),
                       excluded_seeds=excluded_seeds)
        logger.info('Running remap on %s.', self)
        if self.debug_remap:
            debug_file_prefix = os.path.join(scratch_path, 'debug')
        else:
            debug_file_prefix = None
        with WorkDir.using(Path(scratch_path)):
            remap(Path(self.trimmed1_fastq),
                  Path(self.trimmed2_fastq),
                  Path(self.prelim_csv),
                  Path(self.remap_csv),
                  Path(self.remap_counts_csv),
                  Path(self.remap_conseq_csv),
                  Path(self.unmapped1_fastq),
                  Path(self.unmapped2_fastq),
                  debug_file_prefix=debug_file_prefix)

    def run_denovo(self, excluded_seeds):
        logger.info('Running de novo assembly on %s.', self)
        scratch_path = self.get_scratch_path()

        # Set work_dir via dynamic scoping
        with WorkDir.using(Path(self.scratch_path)):
            denovo(Path(self.trimmed1_fastq),
                   Path(self.trimmed2_fastq),
                   Path(self.unstitched_contigs_fasta),
                   Path(self.merged_contigs_csv))

        with open(self.merged_contigs_csv, 'r') as merged_contigs_csv:
            output = Path(self.merged_contigs_fasta)
            try:
                csv_to_fasta(merged_contigs_csv, output)
            except NoContigsInCSV:
                output.touch()

        concatenate_files(inputs=[self.unstitched_contigs_fasta,
                                  self.merged_contigs_fasta],
                          output=self.combined_contigs_fasta)

        with open(self.unstitched_contigs_csv, 'w') as unstitched_contigs_csv, \
             open(self.blast_csv, 'w') as blast_csv:
            fasta_to_csv(Path(self.combined_contigs_fasta),
                         unstitched_contigs_csv,
                         blast_csv=blast_csv,
                         )

        if self.debug_remap:
            debug_file_prefix = os.path.join(scratch_path, 'debug')
        else:
            debug_file_prefix = None

        def with_prefix(path):
            return path and prepend_prefix_to_basename("unstitched_", path)

        logger.info('Running remap on %s.', self)

        with open(self.unstitched_contigs_csv) as contigs_csv, \
                open(with_prefix(self.remap_csv), 'w') as remap_csv, \
                open(with_prefix(self.remap_counts_csv), 'w') as counts_csv, \
                open(with_prefix(self.remap_conseq_csv), 'w') as remap_conseq_csv, \
                open(with_prefix(self.unmapped1_fastq), 'w') as unmapped1, \
                open(with_prefix(self.unmapped2_fastq), 'w') as unmapped2:

            map_to_contigs(self.trimmed1_fastq,
                           self.trimmed2_fastq,
                           contigs_csv,
                           remap_csv,
                           counts_csv,
                           remap_conseq_csv,
                           unmapped1,
                           unmapped2,
                           scratch_path,
                           debug_file_prefix=with_prefix(debug_file_prefix),
                           excluded_seeds=excluded_seeds)

        with open(self.unstitched_contigs_csv, 'r') as unstitched_contigs_csv, \
             open(with_prefix(self.remap_counts_csv)) as unstitched_counts_csv, \
             open(self.contigs_csv, 'w') as contigs_csv:
            referencefull_contig_stitcher(unstitched_contigs_csv, contigs_csv, self.stitcher_plot_svg, unstitched_counts_csv)

        with open(self.contigs_csv) as contigs_csv, \
                open(self.remap_csv, 'w') as remap_csv, \
                open(self.remap_counts_csv, 'w') as counts_csv, \
                open(self.remap_conseq_csv, 'w') as remap_conseq_csv, \
                open(self.unmapped1_fastq, 'w') as unmapped1, \
                open(self.unmapped2_fastq, 'w') as unmapped2:

            map_to_contigs(self.trimmed1_fastq,
                           self.trimmed2_fastq,
                           contigs_csv,
                           remap_csv,
                           counts_csv,
                           remap_conseq_csv,
                           unmapped1,
                           unmapped2,
                           scratch_path,
                           debug_file_prefix=debug_file_prefix,
                           excluded_seeds=excluded_seeds)
