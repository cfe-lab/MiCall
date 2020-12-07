import logging

import os

from micall.resistance.genreport import gen_report
from micall.resistance.resistance import report_resistance

logger = logging.getLogger(__name__)


class EmptySample:
    def __getattr__(self, item):
        return None


class SampleGroup:
    def __init__(self, main_sample, midi_sample=None):
        self.main_sample = main_sample
        self.midi_sample = midi_sample

    def __repr__(self):
        samples = ', '.join(repr(sample)
                            for sample in self)
        return 'SampleGroup({})'.format(samples)

    def __iter__(self):
        yield self.main_sample
        if self.midi_sample is not None:
            yield self.midi_sample

    def process_resistance(self, run_info):
        midi_sample = self.midi_sample or self.main_sample
        logger.info('Running resistance on %s.', self.main_sample)
        with open(self.main_sample.amino_csv) as amino_csv, \
                open(midi_sample.amino_csv) as midi_amino_csv, \
                open(self.main_sample.nuc_csv) as nuc_csv, \
                open(self.main_sample.resistance_csv, 'w') as resistance_csv, \
                open(self.main_sample.mutations_csv, 'w') as mutations_csv, \
                open(self.main_sample.nuc_mutations_csv, 'w') as nuc_mutations_csv, \
                open(self.main_sample.resistance_fail_csv, 'w') as fail_csv, \
                open(self.main_sample.resistance_consensus_csv, 'w') as \
                resistance_consensus_csv:
            report_resistance(amino_csv,
                              midi_amino_csv,
                              nuc_csv,
                              resistance_csv,
                              mutations_csv,
                              nuc_mutations_csv,
                              fail_csv,
                              run_info.reports,
                              resistance_consensus_csv=resistance_consensus_csv)

        logger.info('Running resistance report on %s.', self.main_sample)
        source_path = os.path.dirname(__file__)
        version_filename = os.path.join(source_path,
                                        '..',
                                        '..',
                                        'version.txt')
        if not os.path.exists(version_filename):
            git_version = 'v0-dev'
        else:
            with open(version_filename) as version_file:
                git_version = version_file.read().strip()
        with open(self.main_sample.resistance_csv) as resistance_csv, \
                open(self.main_sample.mutations_csv) as mutations_csv, \
                open(self.main_sample.resistance_pdf, 'wb') as report_pdf:
            gen_report(resistance_csv,
                       mutations_csv,
                       report_pdf,
                       self.main_sample.name,
                       git_version=git_version)
        if not os.stat(self.main_sample.resistance_pdf).st_size:
            os.remove(self.main_sample.resistance_pdf)
        logger.info('Finished resistance of %s.', self.main_sample)
