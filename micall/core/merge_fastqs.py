import argparse
import logging
import shutil
import sys
import os
import tempfile
import gzip
import csv
from typing import TextIO, List, Optional, Dict, IO, Set, Tuple

from micall.core.trim_fastqs import TrimSteps, trim


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger(__name__)


class AmbiguousBadCycles(ValueError): pass
class BadTrimplanFieldNames(ValueError): pass


def concatenate_files(input_files: List[str], output_file: str) -> None:
    with open(output_file, 'wb') as dst:
        for input_file in input_files:
            with open(input_file, 'rb') as src:
                shutil.copyfileobj(src, dst)


def parse_inputs_and_merge_fastqs(trim_file: TextIO, mergeplan_file: TextIO, zip_file: Optional[TextIO]) -> None:
    mergeplan: Dict[str, List[str]] = {}
    trimplan: Set[Tuple[str, ...]] = set()
    zipped: Set[str] = set()
    bad_cycles: Dict[str, str] = {}

    mergeplan_reader = csv.DictReader(mergeplan_file)
    trim_reader = csv.DictReader(trim_file)
    zip_reader = csv.DictReader(zip_file) if zip_file is not None else None

    for row in mergeplan_reader:
        input_path = row["input"]
        output_path = row["output"]

        if output_path not in mergeplan:
            mergeplan[output_path] = []
        mergeplan[output_path].append(input_path)

    no_badcycles_fields = list(sorted(
        (field for field in trim_reader.fieldnames or [] if field != "bad_cycles"),
        key=lambda field: field.lower()))
    expected_no_badcycles_fields = [f"r{i + 1}" for i in range(len(no_badcycles_fields))]

    if [field.lower() for field in no_badcycles_fields] \
       != expected_no_badcycles_fields:
        raise BadTrimplanFieldNames(f"Bad field names: {no_badcycles_fields}, expected {expected_no_badcycles_fields}")

    for row in trim_reader:
        input_paths = tuple(row[field] for field in no_badcycles_fields)
        trimplan.add(input_paths)
        bad_cycles_path = row.get("bad_cycles", "")
        if bad_cycles_path:
            for input_path in input_paths:
                bad_cycles[input_path] = bad_cycles_path

    if zip_reader is not None:
        for row in zip_reader:
            zipped.add(row["file"])

    return merge_fastqs(trimplan, mergeplan, zipped, bad_cycles)


def compress_file(input_path: str, output_path: str) -> None:
    with open(input_path, 'rb') as input_file, \
         open(output_path, 'wb') as output_file:
        with gzip.GzipFile(fileobj=output_file, mode='wb') as gzip_file:
            shutil.copyfileobj(input_file, gzip_file)


def uncompress_file(input_path: str, output_path: str) -> None:
    with open(input_path, 'rb') as compressed_file, \
         open(output_path, 'w+b') as ret:
        with gzip.GzipFile(fileobj=compressed_file, mode='rb') as gzip_file:
            shutil.copyfileobj(gzip_file, ret)


def get_transitive(graph: Dict[str, str], key: str) -> str:
    if key in graph:
        return get_transitive(graph, graph[key])
    else:
        return key


def merge_fastqs(trimplan: Set[Tuple[str, ...]],
                 mergeplan: Dict[str, List[str]],
                 zipped: Set[str] = set(),
                 bad_cycles: Dict[str, str] = {},
                 ) -> None:

    inputs = [value for values in mergeplan.values() for value in values]
    outputs = list(mergeplan.keys())
    name_mapping: Dict[str, str] = {}
    temporary: List[IO] = []

    for output_path in outputs:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

    for input_path in inputs:
        if input_path in zipped:
            logger.info('Uncompressing %s.', input_path)
            uncompressed = tempfile.NamedTemporaryFile(mode="w+b")
            uncompressed_path = uncompressed.name
            uncompress_file(input_path, uncompressed_path)
            temporary.append(uncompressed)
            name_mapping[input_path] = uncompressed_path

    for to_trim in trimplan:
        assert len(to_trim) > 0

        trim_outputs: List[str] = []
        trim_inputs: List[str] = []

        all_bad_cycles_paths = set(bad_cycles[path] for path in to_trim if path in bad_cycles)
        if len(all_bad_cycles_paths) == 0:
            bad_cycles_path = None
        elif len(all_bad_cycles_paths) == 1:
            bad_cycles_path = next(iter(all_bad_cycles_paths))
        else:
            raise AmbiguousBadCycles(f"Ambiguous bad_cycles for {to_trim}: {all_bad_cycles_paths}.")

        for path in to_trim:
            path = get_transitive(name_mapping, path)
            tmp = tempfile.NamedTemporaryFile()
            trim_inputs.append(path)
            trim_outputs.append(tmp.name)
            temporary.append(tmp)
            name_mapping[path] = tmp.name

        logger.info('Trimming samples %s.', ','.join(map(repr, to_trim)))
        trim(trim_inputs, bad_cycles_path, trim_outputs, use_gzip=False)

    for output_path in mergeplan:
        input_files = mergeplan[output_path]
        logger.info('Merging results %s to %s.', ','.join(map(repr, input_files)), output_path)
        input_files = [get_transitive(name_mapping, path) for path in input_files]
        tmp = tempfile.NamedTemporaryFile()
        temporary.append(tmp)
        name_mapping[output_path] = tmp.name
        output_path = tmp.name
        concatenate_files(input_files, output_path)

    for output_path in mergeplan:
        concatenated = name_mapping[output_path]
        if output_path in zipped:
            logger.info('Compressing output %s.', output_path)
            compress_file(concatenated, output_path)
        else:
            os.rename(concatenated, output_path)

    for toclose in temporary:
        try: toclose.close()
        except: pass

    logger.info('Done.')


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(description="Combine and filter the FASTQ files from multiple samples into single output files.")
    parser.add_argument("trimplan", type=argparse.FileType('rt'), help="A CSV file containing the lists of files to be trimmed.")
    parser.add_argument("mergeplan", type=argparse.FileType('rt'), help="A CSV file containing merge plan.")
    parser.add_argument("--zipfile", type=argparse.FileType('rt'), help="A CSV file containing a list of files that are compressed or need to be compressed.")
    args = parser.parse_args(argv)

    parse_inputs_and_merge_fastqs(args.trimplan, args.mergeplan, args.zipfile)
    return 0


if __name__ == '__main__': exit(main(sys.argv[1:]))
