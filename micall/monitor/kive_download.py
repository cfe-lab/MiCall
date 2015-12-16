# To run this as a script, use python -m micall.monitor.kive_download

from argparse import ArgumentParser
import logging
import os
import shutil
import tarfile

from kiveapi import KiveAPI
from kiveapi.runstatus import RunStatus
from requests.adapters import HTTPAdapter

from micall.settings import kive_server_url, kive_user, kive_password, home
from micall.core.miseq_logging import init_logging


def kive_login(server_url, user, password):
    kive = KiveAPI(server_url)
    kive.mount('https://', HTTPAdapter(max_retries=20))
    kive.login(user, password)
    return kive


def download_results(kive_runs, results_folder, run_folder):
    """ Retrieve pipeline output files from Kive

    @param kive_runs: a list of sample names and RunStatus objects, in tuples
        [(sample_name, run_status)]
    @param results_folder: the path where the results should be copied
    @param run_folder: the local folder that will hold working files.
    """
    tar_path = os.path.join(run_folder, 'coverage_maps.tar')
    untar_path = os.path.join(run_folder, 'untar')
    coverage_source_path = os.path.join(untar_path, 'coverage_maps')
    coverage_dest_path = os.path.join(results_folder, 'coverage_maps')
    if not os.path.isdir(untar_path):
        os.mkdir(untar_path)
    if not os.path.isdir(coverage_source_path):
        os.mkdir(coverage_source_path)
    os.mkdir(coverage_dest_path)

    for i, (sample_name, kive_run) in enumerate(kive_runs):
        outputs = kive_run.get_results()
        output_names = ['remap_counts',
                        'conseq',
                        'conseq_ins',
                        'failed_read',
                        'nuc',
                        'amino',
                        'coord_ins',
                        'failed_align',
                        'nuc_variants',
                        'g2p',
                        'coverage_scores',
                        'coverage_maps_tar']
        for output_name in output_names:
            dataset = outputs.get(output_name, None)
            if not dataset:
                continue
            if not output_name.endswith('_tar'):
                filename = os.path.join(results_folder, output_name + '.csv')
                with open(filename, 'a') as result_file:
                    for j, line in enumerate(dataset.readlines()):
                        if i == 0 and j == 0:
                            result_file.write('sample,' + line)
                        elif j != 0:
                            result_file.write(sample_name + ',' + line)
            else:
                with open(tar_path, 'wb') as f:
                    dataset.download(f)
                with tarfile.open(tar_path) as tar:
                    tar.extractall(untar_path)
                for image_filename in os.listdir(coverage_source_path):
                    source = os.path.join(coverage_source_path, image_filename)
                    destination = os.path.join(coverage_dest_path, sample_name + '.' + image_filename)
                    shutil.move(source, destination)

    os.rmdir(coverage_source_path)
    os.rmdir(untar_path)
    os.remove(tar_path)


def find_old_runs(kive, **kwargs):
    params = {}
    param_count = 0
    for key, val in kwargs.iteritems():
        if val is not None:
            params['filters[{}][key]'.format(param_count)] = key
            params['filters[{}][val]'.format(param_count)] = val
            param_count += 1
    response = kive.get('/api/runs/status/', params=params)
    response.raise_for_status()
    json = response.json()
    runs = []
    for entry in json:
        status = RunStatus(entry, kive)
        status.json = entry
        sample_name = entry['display_name']
        status_response = kive.get(entry['run_status'])
        status_response.raise_for_status()
        sample_filename = status_response.json()['inputs']['2']['dataset_name']
        sample_name = '_'.join(sample_filename.split('_')[:2])
        runs.append((sample_name, status))
    runs.sort()
    return runs


def main():
    logger = init_logging(os.path.join(home, 'kive_download.log'),
                          file_log_level=logging.INFO,
                          console_log_level=logging.INFO)
    parser = ArgumentParser(description='Download runs from Kive.')
    parser.add_argument('--startafter',
                        '-a',
                        help='Old runs start after "DD Mon YYYY HH:MM".')
    parser.add_argument('--startbefore',
                        '-b',
                        help='Old runs start before "DD Mon YYYY HH:MM".')
    parser.add_argument('--workfolder',
                        '-w',
                        help='Work folder to download temporary files to.')
    parser.add_argument('--resultfolder',
                        '-r',
                        help='Result folder to copy result files to.')

    args = parser.parse_args()
    logger.info('Starting.')
    kive = kive_login(kive_server_url,
                      kive_user,
                      kive_password)
    runs = find_old_runs(kive,
                         startafter=args.startafter,
                         startbefore=args.startbefore)
    unfinished_count = 0
    for sample_name, run in runs:
        progress = run.json.get('run_progress')
        if progress:
            start_time = progress['start']
            end_time = progress['end']
        else:
            start_time = end_time = None
        if end_time is None:
            unfinished_count += 1
        print(run.json['display_name'])
        print('  ' + sample_name)
        print('  {} - {}'.format(start_time, end_time))
    if args.workfolder or args.resultfolder:
        if not args.workfolder:
            parser.error('argument --workfolder is required with --resultfolder')
        if not args.resultfolder:
            parser.error('argument --resultfolder is required with --workfolder')
        if not os.path.isdir(args.workfolder):
            os.makedirs(args.workfolder)
        if not os.path.isdir(args.resultfolder):
            os.makedirs(args.resultfolder)
        download_results(runs, args.resultfolder, args.workfolder)
    logger.info('%d runs found (%d unfinished).', len(runs), unfinished_count)

if __name__ == '__main__':
    main()
