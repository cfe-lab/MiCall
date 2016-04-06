import os
import fnmatch
import json

from micall.core.prelim_map import prelim_map
from argparse import ArgumentParser
import gzip
import shutil
import errno


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
    args.sample_id = arg_map['Input.sample-id']['Content']['Id']

    return args


def unzip(filename):
    """ Extract a FASTQ file, and return the new file name.

    :param str filename: the FASTQ file's name, possibly compressed
    :return: the uncompressed FASTQ file name
    """
    if not filename.endswith('.gz'):
        return filename
    dest = os.path.splitext(filename)[0]

    # neither file type is present
    with gzip.open(filename, 'rb') as zip_src, open(dest, 'w') as fastq_dest:
        shutil.copyfileobj(zip_src, fastq_dest)
    return dest


def process_sample(sample_id, data_path):
    sampleDir = os.path.join(data_path,
                             'input/samples',
                             sample_id,
                             'Data/Intensities/BaseCalls')
    if not os.path.exists(sampleDir):
        sampleDir = os.path.join(data_path,
                                 'input/samples',
                                 sample_id)
    sample_path = None
    for root, _dirs, files in os.walk(sampleDir):
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
    sample_path = unzip(sample_path)
    sample_path2 = unzip(sample_path2)

    sample_out_path = os.path.join(data_path, 'output/appresults/out')
    try:
        os.makedirs(sample_out_path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(sample_out_path):
            raise

    prelim_csv_path = os.path.join(sample_out_path, 'prelim.csv')
    print('Running prelim_map.')
    with open(prelim_csv_path, 'wb') as prelim_csv:
        prelim_map(sample_path,
                   sample_path2,
                   prelim_csv)
    print('Done.')


def main():
    args = parse_args()
    json_path = os.path.join(args.data_path, 'input/AppSession.json')
    with open(json_path, 'rU') as json_file:
        json = parse_json(json_file)

    process_sample(json.sample_id, args.data_path)

if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    from cStringIO import StringIO
    json_file = StringIO("""\
{"Properties": {"Items": [{"Name": "Input.app-session-name",
                           "Content": "MiCall 04/05/2016 3:14:23"},
                          {"Name": "Input.sample-id",
                           "Content": {"Id": "32896881",
                                       "Name": "Se1-lib2-70x"}}]}}
""")
    parse_json(json_file)
