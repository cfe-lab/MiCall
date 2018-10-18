from collections import namedtuple

from micall.utils.sample_sheet_parser import sample_sheet_parser

SampleGroup = namedtuple('SampleGroup', 'enum names')


def find_groups(file_names, sample_sheet_path, included_projects=None):
    """ Group HCV samples with their MIDI partners.

    :param list[str] file_names: a list of FASTQ file names without paths
    :param sample_sheet_path: path to the SampleSheet.csv file
    :param included_projects: project codes to include, or None to include
        all
    """
    with open(sample_sheet_path) as sample_sheet_file:
        run_info = sample_sheet_parser(sample_sheet_file)

    midi_files = {row['sample']: row['filename']
                  for row in run_info['DataSplit']
                  if row['project'] == 'MidHCV'}
    wide_names = {row['filename']: row['sample']
                  for row in run_info['DataSplit']
                  if (row['project'] != 'MidHCV' and
                      (included_projects is None or
                       row['project'] in included_projects))}
    trimmed_names = {'_'.join(file_name.split('_')[:2]): file_name
                     for file_name in file_names}
    for trimmed_name, file_name in sorted(trimmed_names.items()):
        sample_name = wide_names.get(trimmed_name)
        if sample_name is None:
            # Project was not included.
            continue
        midi_trimmed = midi_files.pop(sample_name + 'MIDI', None)
        midi_name = trimmed_names.get(midi_trimmed)
        yield SampleGroup(sample_name, (file_name, midi_name))

    for sample_name, trimmed_name in midi_files.items():
        file_name = trimmed_names.get(trimmed_name)
        if file_name is not None:
            yield SampleGroup(sample_name, (file_name, None))
