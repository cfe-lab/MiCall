import subprocess

class Samtools():
    def __init__(self, version, execname='samtools'):
        self.version = version
        self.execname = execname

        try:
            stdout = subprocess.check_output([self.execname, 'help'])
            version_line = filter(lambda x: x.startswith('Version:'), stdout.split('\n'))
            assert version_line, 'Failed to parse version number from samtools help text'

            local_version = version_line[0].split()[1]
            assert self.version == local_version, 'samtools version incompatibility %s != %s' % (
                self.version, local_version)
        except OSError:
            raise RuntimeError('samtools not found; check if it is installed and in $PATH\n')

    def run(self, cmd, options):
        """
        Generic wrapper function to call samtools.
        :param cmd: <command> argument; e.g., faidx
        :param options: a list of arguments expected by samtools
        :return:
        """
        subprocess.check_call([self.execname, cmd] + options)

