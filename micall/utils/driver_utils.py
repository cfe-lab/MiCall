import re
import os
import errno
import shutil
import logging
from argparse import ArgumentDefaultsHelpFormatter

logger = logging.getLogger(__name__)


class MiCallFormatter(ArgumentDefaultsHelpFormatter):
    def _fill_text(self, text, width, indent):
        import textwrap
        lines = text.splitlines(keepends=True)
        wrapped = ''
        paragraph = ''
        is_done = False
        while not is_done:
            if not lines:
                is_done = True
                is_blank = False
                extra_indent = next_line = ''
            else:
                next_line = lines.pop(0)
                indent_match = re.match(r'\s+', next_line)
                extra_indent = indent_match.group(0) if indent_match else ''
                is_blank = extra_indent == next_line
                if is_blank:
                    extra_indent = ''
            if (is_done or is_blank or extra_indent) and paragraph:
                paragraph = self._whitespace_matcher.sub(' ', paragraph).strip()
                wrapped += textwrap.fill(paragraph,
                                         width,
                                         initial_indent=indent,
                                         subsequent_indent=indent)
                wrapped += os.linesep
                paragraph = ''
            if extra_indent:
                wrapped += textwrap.fill(next_line,
                                         width,
                                         initial_indent=indent,
                                         subsequent_indent=indent+extra_indent)
                wrapped += os.linesep
            elif not is_done and not is_blank:
                paragraph += next_line
        return wrapped


class MiCallArgs:
    """
    Wrapper that performs some path translation on our inputs and outputs.

    Optionally map parameters to their "MIDI" equivalent, which is useful
    because the HCV sample processing mode has two of every input argument; one
    for the main sample and one for the MIDI sample.
    """
    PREFIX_WITH_RUN_FOLDER = [
        "fastq1",
        "fastq2",
        "bad_cycles_csv",
        "results_folder",
    ]
    MIDI_ARGUMENT_PREFIX = "midi_"

    def __init__(self, args, map_midi=False):
        self.original_args = vars(args)
        self.map_midi = map_midi

    def __getattr__(self, arg_name):
        if arg_name.startswith('__'):
            raise AttributeError(arg_name)

        abs_run_folder = os.path.abspath(self.original_args.get("run_folder"))
        if arg_name == "run_folder":
            return abs_run_folder

        mapped_arg_name = self.MIDI_ARGUMENT_PREFIX + arg_name if self.map_midi else arg_name
        resolved_path = self.original_args.get(mapped_arg_name)
        if resolved_path is None:
            return None
        if not os.path.isabs(resolved_path):
            results = self.original_args["results_folder"]
            if not os.path.isabs(results):
                results = os.path.join(abs_run_folder, results)
            io_prefix = (abs_run_folder
                         if arg_name in self.PREFIX_WITH_RUN_FOLDER
                         else results)
            resolved_path = os.path.join(io_prefix, resolved_path)
        return resolved_path


def safe_file_move(src, dst):
    """
    Helper that attempts to move a file from src to dst.

    Because os.rename may fail on certain platforms, we fall back to
    the boneheaded way of copying-and-deleting.
    :param src:
    :param dst:
    :return:
    """
    try:
        os.rename(src, dst)
    except OSError as e:
        if e.errno == errno.EXDEV:
            logger.debug(
                "Failed to rename %s to %s; copying and deleting the original.",
                src,
                dst,
            )
            shutil.copy(src, dst)
            os.remove(src)
        else:
            raise


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(path):
            raise
