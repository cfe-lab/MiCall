import os
import sys

import Tkinter as tk
from ttk import Progressbar
import tkFileDialog

import shutil
import subprocess
from micall.utils.sample_sheet_parser import sample_sheet_parser
from micall.core.prelim_map import prelim_map
from micall.core.remap import remap
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts
from micall.g2p.pssm_lib import Pssm
from micall.g2p.sam_g2p import sam_g2p
from micall.utils.coverage_plots import coverage_plot

from tempfile import gettempdir
from micall.settings import pipeline_version
from multiprocessing import cpu_count
from micall.utils import collate

import re

fastq_re = re.compile('_L001_R[12]_001.fastq')

files_to_collate = (('amino_frequencies.csv', '.amino.csv'),
                    ('collated_conseqs.csv', '.conseq.csv'),
                    ('failed_align.csv', None),
                    ('failed_read.csv', None),
                    ('g2p.csv', None),
                    ('conseq_ins.csv', None),
                    ('coord_ins.csv', None),
                    ('nucleotide_frequencies.csv', '.nuc.csv'),
                    ('nuc_variants.csv', None),
                    ('collated_counts.csv', '.remap_counts.csv'))

def resource_path(target):
    return os.path.join('' if not hasattr(sys, '_MEIPASS') else sys._MEIPASS, target)


class Redirector(object):
    # see https://gist.github.com/K-DawG007/7729555
    def __init__(self, widget):
        self.widget = widget

    def write(self, msg):
        self.widget.config(state=tk.NORMAL)
        self.widget.insert(tk.END, msg, ('ERROR',))
        self.widget.see(tk.END)
        self.widget.config(state=tk.DISABLED)

class MiCall(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        self.__version__ = '0.3'
        self.pssm = Pssm(path_to_lookup=resource_path('micall/g2p/g2p_fpr.txt'),
                         path_to_matrix=resource_path('micall/g2p/g2p.matrix'))

        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.rundir = None  # path to MiSeq run folder containing data
        self.workdir = gettempdir()  # default to temp directory
        os.chdir(self.workdir)

        self.run_info = None
        self.target_files = []

        self.button_frame = tk.Frame(self)
        self.button_frame.pack(side='top')
        self.console_frame = tk.Frame(self)
        self.console_frame.pack(side='top', fill='both', expand=True)

        #self.button_setwd = tk.Button(
        #    self.button_frame, text="Set working directory", command=self.set_workdir
        #)
        #self.button_setwd.grid(row=0, column=0, sticky='W')

        self.button_load = tk.Button(
            self.button_frame, text="Select run folder", command=self.open_files
        )
        self.button_load.grid(row=0, column=0, sticky='W')

        self.button_run = tk.Button(
            self.button_frame, text="Run", command=self.process_files, state=tk.DISABLED
        )
        self.button_run.grid(row=0, column=1, sticky='W')

        self.button_quit = tk.Button(
            self.button_frame, text='Quit', command=self.quit
        )
        self.button_quit.grid(row=0, column=2, sticky='W')

        self.thread_label = tk.Label(self.button_frame, text='#threads')
        self.thread_label.grid(row=0, column=3, sticky='E')

        self.nthreads = tk.IntVar()
        self.nthreads.set(max(1, cpu_count()/2))
        self.thread_select = tk.OptionMenu(self.button_frame, self.nthreads, *[(i+1) for i in range(cpu_count())])
        self.thread_select.grid(row=0, column=4, sticky='E')

        self.progress_bar = Progressbar(self.button_frame, orient='horizontal', length=500, mode='determinate')
        self.progress_bar.grid(row=1, columnspan=5)

        self.console = tk.Text(
            self.console_frame, bg='black', fg='white', cursor='xterm'
        )
        self.console.pack(side='top', fill='both', expand=True)
        self.console.tag_configure('ERROR', foreground="red")

        # redirect stderr to Text widget
        sys.stderr = Redirector(self.console)

        self.write('Welcome to MiCall v%s (pipeline v%s)\n' % (self.__version__, pipeline_version))
        self.write('Default working directory %s\n' % (self.workdir))

    def set_workdir(self):
        """
        Set destination folder to write intermediate and output files.
        :return:
        """
        folder = tkFileDialog.askdirectory()
        if folder == '':
            # cancelled file dialog
            return
        self.workdir = folder
        self.write('Set destination folder to %s\n' % (self.workdir,))
        os.chdir(self.workdir)  # set working directory
        self.button_load.config(state=tk.ACTIVE)

    def write(self, msg):
        """
        Write to Text widget, scroll to bottom, and then disable editing again.
        :param msg:
        :return:
        """
        self.console.config(state=tk.NORMAL)
        self.console.insert(tk.END, msg)
        self.console.see(tk.END)
        self.console.config(state=tk.DISABLED)

    def open_files(self):
        """
        Transfer FASTQ files into working folder, uncompress.
        :return:
        """
        #self.button_setwd.config(state=tk.DISABLED)
        #self.button_load.config(state=tk.DISABLED)
        self.rundir = tkFileDialog.askdirectory()
        self.write('Selected folder %s\n' % (self.rundir,))

        # check for presence of FASTQ files and SampleSheet.csv
        fastq_files = []
        run_info = None
        for root, dirs, files in os.walk(self.rundir):
            for file in files:
                if file == 'SampleSheet.csv' and run_info is None:
                    try:
                        with open(os.path.join(root, file), 'rU') as handle:
                            run_info = sample_sheet_parser(handle)
                    except:
                        raise

                if len(fastq_re.findall(file)) > 0:
                    fastq_files.append(os.path.join(root, file))


        if not fastq_files:  # empty list
            self.write('Error, folder does not seem to contain any FASTQ(.gz) files!\n')
            self.rundir = None
            #self.button_setwd.config(state=tk.ACTIVE)
            self.button_load.config(state=tk.ACTIVE)
            return False

        if run_info is None:
            self.write('Warning, failed to locate run manifest (SampleSheet.csv).\n')

        self.write('Found %d FASTQ files.\n' % (len(fastq_files),))

        # transfer FASTQ.gz files to working folder
        self.target_files = []
        for src in fastq_files:
            filename = os.path.basename(src)
            prefix = filename.split('.')[0]
            dest = os.path.join(self.workdir, filename)

            if '_R1_001' in dest:
                self.target_files.append(dest.replace('.gz', ''))

            if os.path.exists(os.path.join(self.workdir, prefix+'.fastq')):
                # uncompressed file already present
                continue
            elif os.path.exists(os.path.join(self.workdir, prefix+'.fastq.gz')):
                subprocess.check_call(['gunzip', dest])
                continue
            else:
                # neither file type is present
                shutil.copy(src, dest)
                if dest.endswith('.gz'):
                    subprocess.check_call(['gunzip', dest])


        # remove duplicate entries
        self.target_files = list(set(self.target_files))

        self.write('Transferred %d sets of FASTQ files to working directory.\n' % len(self.target_files))
        self.button_run.config(state=tk.ACTIVE)

    def line_count(self, file):
        """
        Count number of records in a FASTQ file.  This is simply the number of
        lines divided by 4, and multiplied by 2 because these are paired FASTQs.
        :param file:
        :return:
        """
        p = subprocess.Popen(['wc', '-l', file], stdout=subprocess.PIPE)
        output = p.communicate()[0]
        return int(output.strip(' \n').split()[0]) / 2

    def callback(self, msg):
        if type(msg) is int:
            self.progress_bar['value'] = msg
        elif type(msg) is str:
            self.write(msg+'\n')

        self.parent.update_idletasks()


    def process_files(self):
        """
        Perform MiCall data processing on FASTQ files in working directory.
        :return:
        """
        working_files = []
        image_paths = []
        prefixes = []

        # look for FASTQ files
        if len(self.target_files) == 0:
            print 'ERROR: No files to process'
            return

        for fastq1 in self.target_files:
            fastq2 = fastq1.replace('_R1_001', '_R2_001')
            if not os.path.exists(fastq2):
                self.write('ERROR: Missing R2 file for', fastq1)
                continue

            prefix = os.path.basename(fastq1).replace('_L001_R1_001.fastq', '')
            prefixes.append(prefix)
            output_csv = fastq1.replace('_L001_R1_001.fastq', '.prelim.csv')
            working_files.append(output_csv)

            self.write('Processing sample %s\n... preliminary mapping\n' % (prefix,))
            nrecords = self.line_count(fastq1)
            self.progress_bar['value'] = 0
            self.progress_bar['maximum'] = nrecords
            self.parent.update_idletasks()  # flush buffer

            with open(output_csv, 'w') as handle:
                prelim_map(fastq1, fastq2, handle, self.workdir, nthreads=self.nthreads.get(), callback=self.callback)

            # prepare file handles for remap stage
            prelim_csv = open(output_csv, 'rU')
            remap_csv = open(os.path.join(self.workdir, prefix+'.remap.csv'), 'w')
            counts_csv = open(os.path.join(self.workdir, prefix+'.remap_counts.csv'), 'w')
            conseq_csv = open(os.path.join(self.workdir, prefix+'.remap_conseq.csv'), 'w')
            unmapped1 = open(os.path.join(self.workdir, prefix+'.unmapped1.fastq'), 'w')
            unmapped2 = open(os.path.join(self.workdir, prefix+'.unmapped2.fastq'), 'w')
            working_files += map(lambda x: x.name, [remap_csv, counts_csv, conseq_csv, unmapped1, unmapped2])

            self.write('... remapping\n')
            self.parent.update_idletasks()
            self.progress_bar['value'] = 0
            remap(fastq1, fastq2, prelim_csv, remap_csv, counts_csv, conseq_csv, unmapped1, unmapped2, self.workdir,
                  nthreads=self.nthreads.get(), callback=self.callback)

            # prepare file handles for conversion from SAM format to alignment
            remap_csv = open(os.path.join(self.workdir, prefix+'.remap.csv'), 'rU')
            aligned_csv = open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'w')
            insert_csv = open(os.path.join(self.workdir, prefix+'.insert.csv'), 'w')
            failed_csv = open(os.path.join(self.workdir, prefix+'.failed.csv'), 'w')
            working_files += map(lambda x: x.name, [aligned_csv, insert_csv, failed_csv])

            self.write('... converting into alignment\n')
            self.parent.update_idletasks()
            sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv, nthreads=self.nthreads.get())

            aligned_csv = open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'rU')
            nuc_csv = open(os.path.join(self.workdir, prefix+'.nuc.csv'), 'w')
            amino_csv = open(os.path.join(self.workdir, prefix+'.amino.csv'), 'w')
            coord_ins_csv = open(os.path.join(self.workdir, prefix+'.coord_ins.csv'), 'w')
            conseq_csv = open(os.path.join(self.workdir, prefix+'.conseq.csv'), 'w')
            failed_align_csv = open(os.path.join(self.workdir, prefix+'.failed_align.csv'), 'w')
            nuc_variants_csv = open(os.path.join(self.workdir, prefix+'.nuc_variants.csv'), 'w')

            self.write('... extracting statistics from alignments\n')
            self.parent.update_idletasks()
            aln2counts(aligned_csv, nuc_csv, amino_csv, coord_ins_csv, conseq_csv, failed_align_csv, nuc_variants_csv,
                       self.workdir)

            self.write('... generating coverage plots\n')
            self.parent.update_idletasks()
            amino_csv = open(os.path.join(self.workdir, prefix+'.amino.csv'), 'rU')
            image_paths += coverage_plot(amino_csv)

            self.write('... performing g2p scoring on samples covering HIV-1 V3\n')
            self.parent.update_idletasks()
            remap_csv = open(os.path.join(self.workdir, prefix+'.remap.csv'), 'rU')
            nuc_csv = open(os.path.join(self.workdir, prefix+'.nuc.csv'), 'rU')
            g2p_csv = open(os.path.join(self.workdir, prefix+'.g2p.csv'), 'w')
            sam_g2p(pssm=self.pssm, remap_csv=remap_csv, nuc_csv=nuc_csv, g2p_csv=g2p_csv)


        # prevent rerun until a new folder is loaded
        self.button_run.config(state=tk.DISABLED)

        savedir = tkFileDialog.askdirectory(title='Select folder to save results')
        # TODO: if user hits cancel by accident, prevent MiCall from deleting result files

        # collate results to results folder
        for target_file, extension in files_to_collate:
            if extension is None:
                extension = '.' + target_file
            collate.collate_named_files(src_dir=self.workdir, sample_list=prefixes, extension=extension,
                                        output_path=os.path.join(savedir, target_file))

        # copy coverage plots
        imagedir = os.path.join(savedir, 'coverage')
        if not os.path.exists(imagedir):
            os.mkdir(imagedir)
        for src in image_paths:
            dest = os.path.join(imagedir, os.path.basename(src))
            shutil.move(src, dest)

        # clean up working files
        for fn in working_files:
            os.remove(fn)

        # clean up raw data in working directory
        for fn in self.target_files:
            os.remove(fn)

        # reactivate buttons
        self.button_load.config(state=tk.ACTIVE)
        #self.button_setwd.config(state=tk.DISABLED)

        self.write('Run complete.\n')


root = tk.Tk()  # parent widget
root.wm_title('MiCall')

app = MiCall(root).pack(fill='both', expand=True)

root.mainloop()  # enter Tk event loop
root.destroy()
