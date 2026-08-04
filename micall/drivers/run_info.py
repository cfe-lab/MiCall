from collections import namedtuple
from xml.etree import ElementTree
from pathlib import Path

import os

ReadSizes = namedtuple('ReadSizes', 'read1 read2 index1 index2')


class RunInfo:
    def __init__(self,
                 sample_groups,
                 reports=None,
                 interop_path=None,
                 scratch_path=None,
                 output_path=None,
                 read_sizes=None,
                 href_app_session=None,
                 is_denovo=False) -> None:
        """ Initialize the run details.

        :param list[Sample] sample_groups: list of sample details
        :param list[str] reports: resistance reports to generate, or None for all
        :param str interop_path: folder to read Interop data from
        :param str scratch_path: folder that holds all the temporary files
        :param str output_path: folder to write all the output files in
        :param ReadSizes read_sizes: details from run parameters
        :param str href_app_session: identifier for BaseSpace
        :param bool is_denovo: should the reads be assembled, instead of
            remapped?
        """
        self.sample_groups = sample_groups
        self.reports = reports
        self.interop_path = interop_path
        self.scratch_path = scratch_path
        self.quality_csv = scratch_path and os.path.join(scratch_path,
                                                         'quality.csv')
        self.bad_cycles_csv = scratch_path and os.path.join(scratch_path,
                                                            'bad_cycles.csv')
        self.output_path = output_path
        self.bad_tiles_csv = output_path and os.path.join(output_path,
                                                          'bad_tiles.csv')
        self.read_sizes = read_sizes
        self.href_app_session=href_app_session
        self.is_denovo = is_denovo

    def get_all_samples(self):
        for sample_group in self.sample_groups:
            yield from sample_group


def parse_read_sizes(run_info_path: Path) -> ReadSizes:
    run_info = ElementTree.parse(run_info_path).getroot()
    read1 = run_info.find('.//Read[@Number="1"][@IsIndexedRead="N"]')
    read2 = run_info.find('.//Read[@IsIndexedRead="N"][last()]')
    index1 = run_info.find('.//Read[@Number="2"][@IsIndexedRead="Y"]')
    index2 = run_info.find('.//Read[@Number="3"][@IsIndexedRead="Y"]')

    if read1 is None or read2 is None or index1 is None:
        raise ValueError(f"Failed to parse {run_info_path}")

    read_sizes = ReadSizes(read1=int(read1.attrib['NumCycles']),
                           read2=int(read2.attrib['NumCycles']),
                           index1=int(index1.attrib['NumCycles']),
                           index2=(int(index2.attrib['NumCycles'])
                                   if index2 is not None
                                   else 0))
    return read_sizes
