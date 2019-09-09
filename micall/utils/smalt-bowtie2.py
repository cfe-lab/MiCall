#!/usr/bin/env python
# This is a wrapper script to replace smalt mapping features with the equivalent
# bowtie2 features. We built this to work around inconsistent mapping when IVA
# calls smalt.
import os
from argparse import ArgumentParser
from subprocess import run, PIPE
from sys import argv
from traceback import print_exc


def delegate(smalt_args):
    smalt_args[0] = 'smalt-original'
    run(smalt_args)


def get_version_parser():
    return ArgumentParser()


def version(parsed_args):
    assert parsed_args is not None
    result = run(['bowtie2', '--version'], stdout=PIPE)
    version_line = result.stdout.decode('utf8').splitlines()[0]
    version_fields = version_line.split()
    print('Version:', 'bowtie2', version_fields[-1])


def get_index_parser():
    parser = ArgumentParser()
    parser.add_argument('-k', dest='wordlen', type=int)
    parser.add_argument('-s', dest='skipstep', type=int)
    parser.add_argument('index_path')
    parser.add_argument('ref_path')
    return parser


def index(parsed_args):
    expected_args = ((9, 1),
                     (19, 11))
    assert (parsed_args.wordlen,
            parsed_args.skipstep) in expected_args, parsed_args
    run(['bowtie2-build', parsed_args.ref_path, parsed_args.index_path])
    for extension in ('.smi', '.sma'):
        with open(parsed_args.index_path + extension, 'w'):
            pass


def get_map_reads_parser():
    parser = ArgumentParser()
    parser.add_argument('-i', dest='insert_max')
    parser.add_argument('-y', dest='min_identical')
    parser.add_argument('-n', dest='nthreads')
    parser.add_argument('-O', dest='reorder', action='store_true')
    parser.add_argument('index_path')
    parser.add_argument('read_path1')
    parser.add_argument('read_path2')
    return parser


def map_reads(parsed_args):
    assert parsed_args.min_identical in ('0.5', '0.9')
    bowtie_args = ['bowtie2',
                   '--maxins', parsed_args.insert_max,
                   '-x', parsed_args.index_path,
                   '-1', parsed_args.read_path1,
                   '-2', parsed_args.read_path2]
    if parsed_args.min_identical == '0.5':
        bowtie_args.append('--local')
    if parsed_args.nthreads is not None:
        bowtie_args.append('--threads')
        bowtie_args.append(parsed_args.nthreads)
    if parsed_args.reorder:
        bowtie_args.append('--reorder')
    if parsed_args.read_path1.endswith('.fa'):
        bowtie_args.append('-f')
    run(bowtie_args)


def main():
    smalt_args = argv[:]
    use_original = True
    if use_original:
        delegate(smalt_args)
        return
    smalt_args.pop(0)
    command = smalt_args.pop(0)
    parser_factories = dict(version=get_version_parser,
                            index=get_index_parser,
                            map=get_map_reads_parser)
    command_functions = dict(version=version,
                             index=index,
                             map=map_reads)
    parser_factory = parser_factories.get(command)
    if parser_factory:
        parser = parser_factory()
        parsed_args, extras = parser.parse_known_args(smalt_args)
    else:
        parsed_args = extras = None
    index_path = getattr(parsed_args, 'index_path', './index.txt')
    log_path = os.path.join(os.path.dirname(index_path),
                            'smalt.log')
    if parsed_args is None or extras:
        with open(log_path, 'a') as f:
            f.write(f'{argv!r}\n')
            if parsed_args is None:
                message = f'Unknown command: {command}'
            else:
                message = f'Extra arguments: {extras!r}'
            f.write(message + '\n')
        exit(message)
    command_function = command_functions[command]
    try:
        command_function(parsed_args)
    except Exception:
        with open(log_path, 'a') as f:
            print_exc(file=f)
        raise


main()
