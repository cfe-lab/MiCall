import subprocess
import os
import sys
import re
import shutil


class AssetWrapper(object):
    """ Wraps a packaged asset, and finds its path. """
    def __init__(self, path):
        dir, fn = os.path.split(path)
        if len(dir) == 0:
            # <path> is probably execname
            self.path = shutil.which(path)
            if self.path is None:
                print("Failed to find executable {}".format(path))
                sys.exit()

        elif len(fn) == 0:
            # <path> is a directory, trailing separator
            print("Incomplete path {}".format(path))
            sys.exit()
        else:
            # user provided a relative or absolute path
            if not os.path.exists(path):
                print("Path {} does not exist".format(path))
                sys.exit()
            self.path = path



class CommandWrapper(AssetWrapper):
    """ Wraps an external tool, and builds the command lines for it. """
    def __init__(self, version, execname, logger=None, *args, **kwargs):
        super(CommandWrapper, self).__init__(path=execname, *args, **kwargs)
        self.version = version
        self.logger = logger

    def build_args(self, args):
        return [self.path] + args

    def check_output(self, args=[], *popenargs, **kwargs):
        """ Run command with arguments and return its output as a byte string.

        See subprocess.check_output() for details.
        @param args: a list of program arguments
        @param popenargs: other positional arguments to pass along
        @param kwargs: keyword arguments to pass along
        @return the command's output
        """
        try:
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            kwargs.setdefault('startupinfo', startupinfo)
        except:
            pass
        kwargs.setdefault('universal_newlines', True)
        kwargs.setdefault('stdin', sys.stdin)
        final_args = self.build_args(args)
        try:
            return subprocess.check_output(final_args, *popenargs, **kwargs)
        except OSError as ex:
            ex.strerror += ' for command {}'.format(final_args)
            raise

    def create_process(self, args=[], *popenargs, **kwargs):
        """ Execute a child program in a new process.

        See subprocess.Popen class for details.
        @param args: a list of program arguments
        @param popenargs: other positional arguments to pass along
        @param kwargs: keyword arguments to pass along
        @return the new Popen object
        """
        try:
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            kwargs.setdefault('startupinfo', startupinfo)
        except:
            pass
        kwargs.setdefault('universal_newlines', True)
        kwargs.setdefault('stdin', sys.stdin)
        return subprocess.Popen(self.build_args(args), *popenargs, **kwargs)

    def check_logger(self):
        """ Raise an exception if no logger is set for this command. """
        if self.logger is None:
            raise RuntimeError('logger not set for command {}'.format(self.path))

    def log_call(self, args, format_string='%s'):
        """ Launch a subprocess, and log any output to the debug logger.

        Raise an exception if the return code is not zero. This assumes only a
        small amount of output, and holds it all in memory before logging it.
        Logged output includes both stdout and stderr.
        @param args: A list of arguments to pass to check_output().
        @param logger: the logger to log to
        @param format_string: A template for the debug message that will have
        each line of output formatted with it.
        """
        self.check_logger()
        output = self.check_output(args, stderr=subprocess.STDOUT)
        for line in output.splitlines():
            self.logger.debug(format_string, line)

    def yield_output(self, args, *popenargs, **kwargs):
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
        for line in p.stdout:
            yield line
        p.wait()
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode,
                                                self.build_args(args))

    def redirect_call(self, args, outpath, format_string='%s', ignored=None):
        """ Launch a subprocess, and redirect the output to a file.

        Raise an exception if the return code is not zero.
        Standard error is logged to the warn logger.
        @param args: A list of arguments to pass to subprocess.Popen().
        @param outpath: a filename that stdout should be redirected to. If you
        don't need to redirect the output, then just use subprocess.check_call().
        @param format_string: A template for the debug message that will have each
        line of standard error formatted with it.
        @param ignored: A regular expression pattern for stderr messages that
        should not be logged.
        """
        self.check_logger()

        with open(outpath, 'w') as outfile:
            p = self.create_process(args, stdout=outfile, stderr=subprocess.PIPE)
            for line in p.stderr:
                if not ignored or not re.search(ignored, line):
                    self.logger.warn(format_string, line.rstrip())
            p.wait()
            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode,
                                                    self.build_args(args))

    def validate_version(self, version_found):
        if self.version != version_found:
            message = '{} version incompatibility: expected {}, found {}'.format(
                self.path,
                self.version,
                version_found)
            raise RuntimeError(message)


class Bowtie2(CommandWrapper):
    def __init__(self, version=None, execname='bowtie2', logger=None, *args, **kwargs):
        super(Bowtie2, self).__init__(version, execname, logger, *args, **kwargs)
        stdout = self.check_output(['--version'], stderr=subprocess.STDOUT)
        self.version = stdout.split('\n')[0].split()[-1]
        #self.validate_version(version_found)


class Bowtie2Build(CommandWrapper):
    def __init__(self,
                 version=None,
                 execname='bowtie2-build',
                 logger=None,
                 *args,
                 **kwargs):
        super(Bowtie2Build, self).__init__(version,
                                           execname,
                                           logger,
                                           *args,
                                           **kwargs)
        stdout = self.check_output(['--version'], stderr=subprocess.STDOUT)
        self.version = stdout.split('\n')[0].split()[-1]
        #self.validate_version(version_found)

    def build(self, ref_path, reffile_template):
        """ Build an index from a reference file.

        @param ref_path: path to a FASTA file with reference sequences. Must
            be small enough to use with a small index (4GB or less).
        @param reffile_template: file name template for the index files.
        """
        SMALL_INDEX_MAX_SIZE = 4 * 1024**3 - 200  # From bowtie2-build wrapper
        assert os.stat(ref_path).st_size <= SMALL_INDEX_MAX_SIZE
        self.check_logger()
        for line in self.yield_output(['--wrapper',
                                       'micall-0',
                                       '--quiet',
                                       '-f',
                                       ref_path,
                                       reffile_template],
                                      stderr=subprocess.STDOUT):
            if line != 'Building a SMALL index\n':
                self.logger.debug(line)


class LineCounter():
    """ Run the wc command to count lines in a file.

    Fall back to pure Python if wc command is not available.

    Inspired by zed: https://gist.github.com/zed/0ac760859e614cd03652
    """
    def __init__(self):
        self.command = 'wc'

    def count(self, filename, gzip=False):
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
            except:
                self.command = None
        return self.buffered_count(filename)

    def buffered_count(self, filename):
        with open(filename) as f:
            lines = 0
            buf_size = 1024 * 1024
            read_f = f.read  # loop optimization

            buf = read_f(buf_size)
            while buf:
                lines += buf.count('\n')
                buf = read_f(buf_size)

        return lines
