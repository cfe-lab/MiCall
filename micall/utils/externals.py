import subprocess
import os
import re
from subprocess import CalledProcessError
from pathlib import Path
import logging
from typing import Optional, List, Any, Iterator
from functools import cached_property
from abc import ABC, abstractmethod, abstractproperty
import shutil
import importlib.resources as resources
import contextlib


class ExternalResource(ABC):
    @abstractproperty
    def identifier(sef) -> str: ...

    @abstractmethod
    @contextlib.contextmanager
    def path(self) -> Iterator[Path]: ...


class CommandWrapper(ExternalResource):
    """ Wraps an external tool, and builds the command lines for it. """

    def __init__(self) -> None:
        self._logger: Optional[logging.Logger] = None
        self._validate_version()

    @abstractproperty
    def executable_name(self) -> str: ...

    @abstractproperty
    def expected_version(self) -> Optional[str]: ...

    @property
    def identifier(self) -> str:
        return self.executable_name

    @cached_property
    def executable_path(self) -> Path:
        path = shutil.which(self.executable_name)
        if path is not None:
            return Path(path)

        raise RuntimeError(f"Cannot find {self.identifier!r} executable."
                           "Make sure you have installed the package.")

    @contextlib.contextmanager
    def path(self) -> Iterator[Path]:
        yield self.executable_path

    def _validate_version(self) -> None:
        if self.expected_version is None:
            pass
        elif self.expected_version != self.version:
            message = '{} version incompatibility: expected {}, found {}'.format(
                self.identifier,
                self.expected_version,
                self.version)
            raise RuntimeError(message)

    @abstractmethod
    def get_version(self) -> str: ...

    @cached_property
    def version(self) -> str:
        return self.get_version()

    @property
    def logger(self) -> logging.Logger:
        """ Raise an exception if no logger is set for this command. """
        if self._logger is None:
            raise RuntimeError(f'logger not set for command {self.identifier!r}')
        return self._logger

    def set_logger(self, logger: logging.Logger) -> None:
        self._logger = logger

    def build_args(self, args: List[str]) -> List[str]:
        return [str(self.executable_path)] + args

    def check_output(self, args: List[str] = [],
                     *popenargs: Any, **kwargs: Any) -> str:
        """ Run command with arguments and return its output as a byte string.

        See subprocess.check_output() for details.
        @param args: a list of program arguments
        @param popenargs: other positional arguments to pass along
        @param kwargs: keyword arguments to pass along
        @return the command's output
        """
        try:
            startupinfo = subprocess.STARTUPINFO()  # type: ignore[attr-defined]

            # Needed on Windows, fails elsewhere.
            # noinspection PyUnresolvedReferences
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW  # type: ignore[attr-defined]
            kwargs.setdefault('startupinfo', startupinfo)
        except AttributeError:
            pass
        kwargs.setdefault('universal_newlines', True)
        final_args = self.build_args(args)
        try:
            return subprocess.check_output(final_args, *popenargs, **kwargs)
        except OSError as ex:
            original_error = ex.strerror or 'Failed'
            ex.strerror = f'{original_error} for command {final_args}.'
            raise

    def create_process(self, args: List[str] = [],
                       *popenargs: Any, **kwargs: Any) -> subprocess.Popen:
        """ Execute a child program in a new process.

        See subprocess.Popen class for details.
        @param args: a list of program arguments
        @param popenargs: other positional arguments to pass along
        @param kwargs: keyword arguments to pass along
        @return the new Popen object
        """
        try:
            startupinfo = subprocess.STARTUPINFO()  # type: ignore[attr-defined]
            # Needed on Windows, fails elsewhere.
            # noinspection PyUnresolvedReferences
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW  # type: ignore[attr-defined]
            kwargs.setdefault('startupinfo', startupinfo)
        except AttributeError:
            pass
        kwargs.setdefault('universal_newlines', True)
        with open(os.devnull) as devnull:
            kwargs.setdefault('stdin', devnull)
            return subprocess.Popen(self.build_args(args or []), *popenargs, **kwargs)

    def log_call(self, args: List[str], format_string: str = '%s') -> None:
        """ Launch a subprocess, and log any output to the debug logger.

        Raise an exception if the return code is not zero. This assumes only a
        small amount of output, and holds it all in memory before logging it.
        Logged output includes both stdout and stderr.
        @param args: A list of arguments to pass to check_output().
        @param format_string: A template for the debug message that will have
        each line of output formatted with it.
        """

        output = self.check_output(args, stderr=subprocess.STDOUT)
        for line in output.splitlines():
            self.logger.debug(format_string, line)

    def yield_output(self, args: List[str], *popenargs: Any, **kwargs: Any
                     ) -> Iterator[str]:
        """ Launch a subprocess, and yield the lines of standard output.

        Raise an exception if the return code is not zero.
        Standard error is written to standard error and not returned, unless
        you specify stderr=subprocessing.STDOUT as a keyword argument.
        @param args: A list of arguments to pass to subprocess.Popen().
        @param popenargs: other positional arguments to pass along
        @param kwargs: keyword arguments to pass along
        """
        p = self.create_process(args,
                                stdout=subprocess.PIPE,
                                *popenargs,
                                **kwargs)
        assert p.stdout is not None
        for line in p.stdout:
            yield line
        p.wait()
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode,
                                                self.build_args(args))

    def redirect_call(self,
                      args: List[str],
                      outpath: Path,
                      format_string: str = '%s',
                      ignored: Optional[re.Pattern] = None) -> None:
        """ Launch a subprocess, and redirect the output to a file.

        Raise an exception if the return code is not zero.
        Standard error is logged to the warn logger.
        @param args: A list of arguments to pass to subprocess.Popen().
        @param outpath: a filepath that stdout should be redirected to. If you
        don't need to redirect the output, then just use subprocess.check_call().
        @param format_string: A template for the debug message that will have each
        line of standard error formatted with it.
        @param ignored: A regular expression pattern for stderr messages that
        should not be logged.
        """

        with open(outpath, 'w') as outfile:
            p = self.create_process(args, stdout=outfile, stderr=subprocess.PIPE)
            assert p.stderr is not None
            for line in p.stderr:
                if not ignored or not re.search(ignored, line):
                    self.logger.warn(format_string, line.rstrip())
            p.wait()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode,
                                                    self.build_args(args))


class CutAdapt(CommandWrapper):
    @property
    def executable_name(self) -> str:
        return 'cutadapt'

    @property
    def expected_version(self) -> Optional[str]:
        return None

    def get_version(self):
        stdout = self.check_output(['--version'], stderr=subprocess.STDOUT)
        return stdout.split('\n')[0].split()[-1]


class Bowtie2(CommandWrapper):
    @property
    def executable_name(self) -> str:
        return 'bowtie2-align-s'

    @property
    def expected_version(self) -> Optional[str]:
        return '2.2.8'

    def get_version(self):
        stdout = self.check_output(['--version'], stderr=subprocess.STDOUT)
        return stdout.split('\n')[0].split()[-1]


class Bowtie2Build(CommandWrapper):
    @property
    def executable_name(self) -> str:
        return 'bowtie2-build-s'

    @property
    def expected_version(self) -> Optional[str]:
        return '2.2.8'

    def get_version(self):
        stdout = self.check_output(['--version'], stderr=subprocess.STDOUT)
        return stdout.split('\n')[0].split()[-1]

    def build(self, ref_path, reffile_template):
        """ Build an index from a reference file.

        @param ref_path: path to a FASTA file with reference sequences. Must
            be small enough to use with a small index (4GB or less).
        @param reffile_template: file name template for the index files.
        """

        small_index_max_size = 4 * 1024**3 - 200  # From bowtie2-build wrapper
        assert os.stat(ref_path).st_size <= small_index_max_size
        lines = []
        try:
            for line in self.yield_output(['--wrapper',
                                           'micall-0',
                                           '--quiet',
                                           '-f',
                                           ref_path,
                                           reffile_template],
                                          stderr=subprocess.STDOUT):
                lines.append(line.rstrip())
        except Exception:
            for line in lines:
                self.logger.error(line)
            raise
        for line in lines:
            if line != 'Building a SMALL index\n':
                self.logger.debug(line)


class LineCounter:
    """ Run the wc command to count lines in a file.

    Fall back to pure Python if wc command is not available.

    Inspired by zed: https://gist.github.com/zed/0ac760859e614cd03652
    """
    def __init__(self):
        self.command: str = 'wc'

    def count(self, filename: Path, gzip: bool = False) -> int:
        if self.command:
            try:
                if gzip:
                    p = subprocess.Popen(['gunzip', '-c', filename], stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT)
                    wc_output = subprocess.check_output(['wc', '-l'], stdin=p.stdout,
                                                        stderr=subprocess.STDOUT)
                else:
                    wc_output = subprocess.check_output(['wc', '-l', filename],
                                                        stderr=subprocess.STDOUT)

                return int(wc_output.split()[0])
            except CalledProcessError:
                self.command = ''
        return self.buffered_count(filename)

    def buffered_count(self, filename: Path) -> int:
        with open(filename) as f:
            lines = 0
            buf_size = 1024 * 1024
            read_f = f.read  # loop optimization

            buf = read_f(buf_size)
            while buf:
                lines += buf.count('\n')
                buf = read_f(buf_size)

        return lines


class Blastn(CommandWrapper):
    BLAST_COLUMNS = ['qaccver',
                     'saccver',
                     'qlen',
                     'pident',
                     'score',
                     'qcovhsp',
                     'qstart',
                     'qend',
                     'sstart',
                     'send',
                     ]

    @property
    def executable_name(self) -> str:
        return 'blastn'

    @property
    def expected_version(self) -> Optional[str]:
        return None

    def get_version(self):
        stdout = self.check_output(['-version'], stderr=subprocess.STDOUT)
        return stdout.split('\n')[0].split()[-1]

    def genotype(self, database: Path, contigs_fasta: Path) -> str:
        blast_format = '10 ' + " ".join(Blastn.BLAST_COLUMNS)
        program = ["-outfmt", blast_format,
                   "-query", str(contigs_fasta),
                   "-db", str(database),
                   "-evalue", "0.0001",
                   "-max_target_seqs", "5000",
                   "-gapopen", "5",
                   "-gapextend", "2",
                   "-penalty", "-3",
                   "-reward", "1",
                   ]
        stdout = self.check_output(program)
        return stdout


class PackagedResource(ExternalResource):
    @property
    def identifier(self) -> str:
        return 'micall://' + str(self.relative_path)

    @contextlib.contextmanager
    def path(self) -> Iterator[Path]:
        with PackagedResource.root_dir() as root:
            yield root / self.relative_path

    @abstractproperty
    def relative_path(self) -> Path: ...

    @staticmethod
    @contextlib.contextmanager
    def root_dir() -> Iterator[Path]:
        """
        A context manager handling the path to packaged MiCall directory.

        The complexity of the function arises from the need to maintain
        compatibility with multiple python versions due to changes in APIs
        of the `importlib.resources` package.

        It first tries to fetch the resource using `resources.files`
        function introduced in Python 3.9. If it fails, it falls back on
        `resources.path`.  It further ensures that the obtained resource
        is returned as a Path instance regardless of it being a string,
        Path, or contextlib context-manager instance.

        Note: `resources.path` is set to be deprecated in future Python
        versions, hence the intended primary method is using
        `resources.files`.

        Yields: Path: A path-like object pointing to
                  the root directory of 'micall' package.
        """

        try:
            ret = resources.as_file(resources.files('micall'))  # type: ignore
        except AttributeError:
            ret = resources.path('micall')  # type: ignore

        if isinstance(ret, str):
            yield Path(ret)
        elif isinstance(ret, Path):
            yield ret
        else:
            with ret as path:
                yield path


class DefaultBlastDatabase(PackagedResource):
    @property
    def relative_path(self) -> Path:
        return Path(".") / "blast_db" / "refs.fasta"


class ProjectsFile(PackagedResource):
    @property
    def relative_path(self) -> Path:
        return Path(".") / "projects.json"


class ProjectsScoringFile(PackagedResource):
    @property
    def relative_path(self) -> Path:
        return Path(".") / "project_scoring.json"
