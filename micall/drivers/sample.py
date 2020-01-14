""" Data needed to process a sample """
import logging
import os

from micall.core.aln2counts import aln2counts
from micall.core.amplicon_finder import write_merge_lengths_plot, merge_for_entropy
from micall.core.cascade_report import CascadeReport
from micall.core.coverage_plots import coverage_plot
from micall.core.plot_contigs import plot_genome_coverage
from micall.core.prelim_map import prelim_map
from micall.core.remap import remap, map_to_contigs
from micall.core.sam2aln import sam2aln
from micall.core.trim_fastqs import trim
from micall.core.denovo import denovo
from micall.g2p.fastq_g2p import fastq_g2p, DEFAULT_MIN_COUNT, MIN_VALID, MIN_VALID_PERCENT

logger = logging.getLogger(__name__)


class Sample:
    def __init__(self,
                 basespace_id=None,
                 basespace_href=None,
                 rank=None,
                 debug_remap=False,
                 scratch_path=None,
                 **paths):
        """ Record the details.

        :param str fastq1: path to first FASTQ file
        :param str fastq2: path to second FASTQ file
        :param str basespace_id: id from BaseSpace
        :param str basespace_href: web address within BaseSpace API
        :param str rank: rank within the run, like '2 of 10'
        :param bool debug_remap: write debug files during remap step
        :param str scratch_path: path where unknown outputs are written
        :param paths: {output_name: file_path}
        """
        fastq1 = paths.get('fastq1')
        if 'fastq2' in paths:
            pass
        elif 'fastq1' in paths:
            if '_R1_' not in fastq1:
                raise ValueError(
                    "fastq2 not given, and fastq1 does not contain '_R1_'.")
            paths['fastq2'] = fastq1.replace('_R1_', '_R2_')
        if fastq1:
            self.name = '_'.join(os.path.basename(fastq1).split('_')[:2])
        else:
            self.name = None
        self.basespace_id = basespace_id
        self.basespace_href = basespace_href
        self.rank = rank
        self.debug_remap = debug_remap
        self.scratch_path = scratch_path
        self.paths = paths

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
        for extension in ('csv', 'fastq', 'pdf', 'svg', 'png'):
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
                excluded_seeds=(),
                excluded_projects=(),
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
        os.mkdir(scratch_path)
        use_gzip = force_gzip or self.fastq1.endswith('.gz')

        with open(self.read_summary_csv, 'w') as read_summary:
            trim((self.fastq1, self.fastq2),
                 self.bad_cycles_csv,
                 (self.trimmed1_fastq, self.trimmed2_fastq),
                 summary_file=read_summary,
                 use_gzip=use_gzip)

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

        logger.info('Running sam2aln on %s.', self)
        with open(self.remap_csv) as remap_csv, \
                open(self.aligned_csv, 'w') as aligned_csv, \
                open(self.conseq_ins_csv, 'w') as conseq_ins_csv, \
                open(self.failed_csv, 'w') as failed_csv, \
                open(self.clipping_csv, 'w') as clipping_csv:

            sam2aln(remap_csv,
                    aligned_csv,
                    conseq_ins_csv,
                    failed_csv,
                    clipping_csv=clipping_csv)

        logger.info('Running aln2counts on %s.', self)
        if use_denovo:
            contigs_path = self.contigs_csv
        else:
            contigs_path = os.devnull
        with open(self.aligned_csv) as aligned_csv, \
                open(self.g2p_aligned_csv) as g2p_aligned_csv, \
                open(self.clipping_csv) as clipping_csv, \
                open(self.conseq_ins_csv) as conseq_ins_csv, \
                open(self.remap_conseq_csv) as remap_conseq_csv, \
                open(contigs_path) as contigs_csv, \
                open(self.nuc_csv, 'w') as nuc_csv, \
                open(self.nuc_detail_csv, 'w') as nuc_detail_csv, \
                open(self.amino_csv, 'w') as amino_csv, \
                open(self.amino_detail_csv, 'w') as amino_detail_csv, \
                open(self.coord_ins_csv, 'w') as coord_ins_csv, \
                open(self.conseq_csv, 'w') as conseq_csv, \
                open(self.conseq_region_csv, 'w') as conseq_region_csv, \
                open(self.failed_align_csv, 'w') as failed_align_csv, \
                open(self.coverage_summary_csv, 'w') as coverage_summary_csv, \
                open(self.genome_coverage_csv, 'w') as genome_coverage_csv:

            if not use_denovo:
                for f in (amino_detail_csv, nuc_detail_csv):
                    f.close()
                    os.remove(f.name)
                amino_detail_csv = nuc_detail_csv = None
            aln2counts(aligned_csv,
                       nuc_csv,
                       amino_csv,
                       coord_ins_csv,
                       conseq_csv,
                       failed_align_csv,
                       coverage_summary_csv=coverage_summary_csv,
                       clipping_csv=clipping_csv,
                       conseq_ins_csv=conseq_ins_csv,
                       g2p_aligned_csv=g2p_aligned_csv,
                       remap_conseq_csv=remap_conseq_csv,
                       conseq_region_csv=conseq_region_csv,
                       amino_detail_csv=amino_detail_csv,
                       nuc_detail_csv=nuc_detail_csv,
                       genome_coverage_csv=genome_coverage_csv,
                       contigs_csv=contigs_csv)

        logger.info('Running coverage_plots on %s.', self)
        os.makedirs(self.coverage_maps)
        with open(self.amino_csv) as amino_csv, \
                open(self.coverage_scores_csv, 'w') as coverage_scores_csv:
            coverage_plot(amino_csv,
                          coverage_scores_csv,
                          coverage_maps_path=self.coverage_maps,
                          coverage_maps_prefix=self.name,
                          excluded_projects=excluded_projects)

        if use_denovo:
            blast_path = self.blast_csv
        else:
            blast_path = os.devnull
        with open(self.genome_coverage_csv) as genome_coverage_csv, \
                open(blast_path) as blast_csv:
            if not use_denovo:
                blast_csv = None
            plot_genome_coverage(genome_coverage_csv,
                                 blast_csv,
                                 self.genome_coverage_svg)

        logger.info('Running cascade_report on %s.', self)
        with open(self.g2p_summary_csv) as g2p_summary_csv, \
                open(self.remap_counts_csv) as remap_counts_csv, \
                open(self.aligned_csv) as aligned_csv, \
                open(self.cascade_csv, 'w') as cascade_csv:
            cascade_report = CascadeReport(cascade_csv)
            cascade_report.g2p_summary_csv = g2p_summary_csv
            cascade_report.remap_counts_csv = remap_counts_csv
            cascade_report.aligned_csv = aligned_csv
            cascade_report.generate()
        logger.info('Finished sample %s.', self)

    def run_mapping(self, excluded_seeds):
        logger.info('Running prelim_map on %s.', self)
        scratch_path = self.get_scratch_path()
        with open(self.prelim_csv, 'w') as prelim_csv:
            prelim_map(self.trimmed1_fastq,
                       self.trimmed2_fastq,
                       prelim_csv,
                       work_path=scratch_path,
                       excluded_seeds=excluded_seeds)
        logger.info('Running remap on %s.', self)
        if self.debug_remap:
            debug_file_prefix = os.path.join(scratch_path, 'debug')
        else:
            debug_file_prefix = None
        with open(self.prelim_csv) as prelim_csv, \
                open(self.remap_csv, 'w') as remap_csv, \
                open(self.remap_counts_csv, 'w') as counts_csv, \
                open(self.remap_conseq_csv, 'w') as conseq_csv, \
                open(self.unmapped1_fastq, 'w') as unmapped1, \
                open(self.unmapped2_fastq, 'w') as unmapped2:

            remap(self.trimmed1_fastq,
                  self.trimmed2_fastq,
                  prelim_csv,
                  remap_csv,
                  counts_csv,
                  conseq_csv,
                  unmapped1,
                  unmapped2,
                  scratch_path,
                  debug_file_prefix=debug_file_prefix)

    def run_denovo(self, excluded_seeds):
        logger.info('Running de novo assembly on %s.', self)
        scratch_path = self.get_scratch_path()
        with open(self.merged_contigs_csv) as merged_contigs_csv, \
                open(self.contigs_csv, 'w') as contigs_csv, \
                open(self.blast_csv, 'w') as blast_csv:
            denovo(self.trimmed1_fastq,
                   self.trimmed2_fastq,
                   contigs_csv,
                   self.scratch_path,
                   merged_contigs_csv,
                   blast_csv=blast_csv)
        logger.info('Running remap on %s.', self)
        if self.debug_remap:
            debug_file_prefix = os.path.join(scratch_path, 'debug')
        else:
            debug_file_prefix = None
        with open(self.contigs_csv) as contigs_csv, \
                open(self.remap_csv, 'w') as remap_csv, \
                open(self.remap_counts_csv, 'w') as counts_csv, \
                open(self.remap_conseq_csv, 'w') as conseq_csv, \
                open(self.unmapped1_fastq, 'w') as unmapped1, \
                open(self.unmapped2_fastq, 'w') as unmapped2:

            map_to_contigs(self.trimmed1_fastq,
                           self.trimmed2_fastq,
                           contigs_csv,
                           remap_csv,
                           counts_csv,
                           conseq_csv,
                           unmapped1,
                           unmapped2,
                           scratch_path,
                           debug_file_prefix=debug_file_prefix,
                           excluded_seeds=excluded_seeds)
