from argparse import ArgumentParser
import errno
import fnmatch
import gzip
import json
import os
import shutil
from tempfile import NamedTemporaryFile

from micall.core.prelim_map import prelim_map


def parse_args():
    parser = ArgumentParser(description='Map FASTQ files to references.')
    parser.add_argument('data_path',
                        default='/data',
                        nargs='?',
                        help='data folder filled in by BaseSpace')
    return parser.parse_args()


def parse_json(json_file):
    """ Load JSON from an open file, and pull out the arguments for this run.

    :param json_file: an open file that contains JSON in the BaseSpace
    AppSession format.
    :return: an object with an attribute for each argument
    """
    class Args(object):
        pass
    args = Args()

    raw_args = json.load(json_file)
    arg_map = {item['Name']: item
               for item in raw_args['Properties']['Items']}

    args.name = arg_map['Input.app-session-name']['Content']
    args.run_id = arg_map['Input.run-id']['Content']['Id']
    args.samples = arg_map['Input.sample-ids']['Items']
    args.project_id = arg_map['Input.project-id']['Content']['Id']

    return args


def unzip(filename, scratch_path):
    """ Extract a FASTQ file, and return the new file name.

    :param str filename: the FASTQ file's name, possibly compressed
    :return: the uncompressed FASTQ file name
    """
    if not filename.endswith('.gz'):
        return filename

    gz_name = os.path.basename(filename)  # drop path
    fastq_name = os.path.splitext(gz_name)[0]  # drop .gz
    basename = os.path.splitext(fastq_name)[0]  # drop .fastq

    with gzip.open(filename, 'rb') as zip_src, NamedTemporaryFile(
            mode='w',
            prefix=basename,
            suffix='.fastq',
            dir=scratch_path,
            delete=False) as dest:
        shutil.copyfileobj(zip_src, dest)
        return dest.name


def process_sample(sample_info, project_id, data_path):
    scratch_path = os.path.join(data_path, 'scratch')
    makedirs(scratch_path)
    sample_id = sample_info['Id']
    sample_name = sample_info['Name']
    sample_dir = os.path.join(data_path,
                              'input/samples',
                              sample_id,
                              'Data/Intensities/BaseCalls')
    if not os.path.exists(sample_dir):
        sample_dir = os.path.join(data_path,
                                  'input/samples',
                                  sample_id)
    sample_path = None
    for root, _dirs, files in os.walk(sample_dir):
        sample_paths = fnmatch.filter(files, '*_R1_*')
        if sample_paths:
            sample_path = os.path.join(root, sample_paths[0])
            break
    if sample_path is None:
        raise RuntimeError('No R1 file found for sample id {}.'.format(sample_id))
    sample_path2 = sample_path.replace('_R1_', '_R2_')
    if not os.path.exists(sample_path2):
        raise RuntimeError('R2 file missing for sample id {}: {!r}.'.format(
            sample_id,
            sample_path2))
    print('Processing sample {}.'.format(sample_path))
    sample_path = unzip(sample_path, scratch_path)
    sample_path2 = unzip(sample_path2, scratch_path)

    sample_out_path = os.path.join(data_path,
                                   'output/appresults',
                                   project_id,
                                   'out')
    makedirs(sample_out_path)

    prelim_csv_path = os.path.join(sample_out_path, sample_name + '_prelim.csv')
    print('Running prelim_map.')
    with open(prelim_csv_path, 'wb') as prelim_csv:
        prelim_map(sample_path,
                   sample_path2,
                   prelim_csv)
    print('Done.')


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(path):
            raise


def main():
    args = parse_args()
    json_path = os.path.join(args.data_path, 'input/AppSession.json')
    with open(json_path, 'rU') as json_file:
        json = parse_json(json_file)

    for sample_info in json.samples:
        process_sample(sample_info, json.project_id, args.data_path)

if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    from cStringIO import StringIO
    json_file = StringIO("""\
{"Properties": {"Items": [{"Name": "Input.app-session-name",
                           "Content": "MiCall 04/05/2016 3:14:23"},
                          {"Name": "Input.sample-ids",
                           "Items": [{"Id": "11111",
                                      "Name": "Example-Sample1"},
                                     {"Id": "22222",
                                      "Name": "Example-Sample2"}]},
                          {"Name": "Input.project-id",
                           "Content": {"Id": "33333",
                                       "Name": "Example-Project"}},
                          {"Name": "Input.run-id",
                           "Content": {"Id": "44444",
                                       "Name": "160115_Project"}}]}}
""")
    parse_json(json_file)
