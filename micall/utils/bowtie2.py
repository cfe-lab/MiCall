"""
Wrapper script for bowtie2-align and bowtie2-build with version control
"""
import subprocess

class Bowtie2():
    def __init__(self, version, align_exec='bowtie2', build_exec='bowtie2-build'):
        self.version = version  # assume same version number for both executables
        self.align_exec = align_exec
        self.build_exec = build_exec

        # checks for bowtie2(-align)
        try:
            subprocess.check_output([self.align_exec, '-h'])
        except OSError:
            raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

        stdout = subprocess.check_output([self.align_exec, '--version'])
        local_version = stdout.split('\n')[0].split()[-1]
        assert self.version == local_version, 'bowtie2 version incompatibility %s != %s' % (
            self.version, local_version)

        # checks for bowtie2-build
        try:
            subprocess.check_output([self.build_exec, '-h'])
        except OSError:
            raise RuntimeError('bowtie2-build not found; check if it is installed and in $PATH\n')

        stdout = subprocess.check_output([self.build_exec, '--version'])
        local_version = stdout.split('\n')[0].split()[-1]
        assert self.version == local_version, 'bowtie2-build version incompatibility %s != %s' % (
            self.version, local_version)

    def build(fasta, refpath=None):
        """
        Construct bowtie2 indices (.bt2 files)
        :param version: Enforces bowtie2 version number (e.g., '2.2.1')
        :param fasta: Path to FASTA containing reference sequences
        :return:
        """
        if refpath is None:
            refpath = fasta
        subprocess.check_call(['bowtie2-build', '-q', '-f', fasta, refpath])

    def align_paired(refpath, fastq1, fastq2, nthreads, flags=('--quiet', '--no-unal', '--local')):
        """
        Call bowtie2-align on paired read data
        :param refpath: Path to bowtie2 index files (.bt2)
        :param fastq1: Files with #1 mates
        :param fastq2: Files with #2 mates
        :param flags: A tuple of bowtie2 flags such as --local
        :return:
        """

        # stream output from bowtie2
        bowtie_args = ['bowtie2', '-x', refpath, '-1', fastq1, '-2', fastq2, '-p', str(nthreads)] + list(flags)
        p2 = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
        for line in p2.stdout:
            yield line

        # exception handling


def main():
    bowtie2 = Bowtie2(version='2.1.3')
    iter = bowtie2.align_paired('gb-ref.fa', 'test1.fastq', 'test2.fastq')
    for line in iter:
        print line

if __name__ == '__main__':
    main()
