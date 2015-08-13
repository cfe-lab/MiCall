import gzip
import json
from multiprocessing import cpu_count
import os
import re
import shutil
from tempfile import gettempdir

import Tkinter as tk
from ttk import Progressbar
import tkFileDialog

from micall.core.prelim_map import prelim_map
from micall.core.remap import remap
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts
from micall.g2p.pssm_lib import Pssm
from micall.g2p.sam_g2p import sam_g2p
from micall.settings import pipeline_version
from micall.utils import collate
from micall.utils.coverage_plots import coverage_plot
from micall.utils.externals import AssetWrapper, LineCounter
from micall.utils.sample_sheet_parser import sample_sheet_parser

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
    CONFIG_FILE = os.path.expanduser("~/.micall.config")
    def __init__(self, parent, *args, **kwargs):
        self.__version__ = '0.3'
        self.pssm = Pssm(path_to_lookup=AssetWrapper('micall/g2p/g2p_fpr.txt').path,
                         path_to_matrix=AssetWrapper('micall/g2p/g2p.matrix').path)

        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.rundir = None  # path to MiSeq run folder containing data
        self.workdir = gettempdir()  # default to temp directory
        os.chdir(self.workdir)
        
        self.line_counter = LineCounter()

        self.run_info = None
        self.target_files = []

        self.button_frame = tk.Frame(self)
        self.button_frame.pack(side='top')
        self.console_frame = tk.Frame(self)
        self.console_frame.pack(side='top', fill='both', expand=True)
        
        try:
            with open(MiCall.CONFIG_FILE, 'rU') as f:
                self.config = json.load(f)
        except:
            self.config = {}

        #self.button_setwd = tk.Button(
        #    self.button_frame, text="Set working directory", command=self.set_workdir
        #)
        #self.button_setwd.grid(row=0, column=0, sticky='W')

        self.button_run = tk.Button(
            self.button_frame, text="Run", command=self.process_files
        )
        self.button_run.grid(row=0, column=1, sticky='W')

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
        #sys.stderr = Redirector(self.console)

        self.write('Welcome to MiCall v%s (pipeline v%s)\n' % (self.__version__, pipeline_version))

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
        
    def write_config(self):
        try:
            with open(MiCall.CONFIG_FILE, 'w') as f:
                json.dump(self.config, f)
        except:
            pass # For now, we don't care if config fails

    def open_files(self):
        """
        Transfer FASTQ files into working folder, uncompress.
        :return:
        """
        fastq_files = []
        setting_name = 'run_path'
        run_path = self.config.get(setting_name, '')
        self.rundir = tkFileDialog.askdirectory(
            title='Choose a folder of FASTQ files:',
            initialdir=run_path)
        if not self.rundir:
            return fastq_files
        self.config[setting_name] = self.rundir
        self.write_config()
        self.write('Selected folder %s\n' % (self.rundir,))

        # check for presence of FASTQ files and SampleSheet.csv
        run_info = None
        for root, _dirs, files in os.walk(self.rundir):
            for file in files:
                if file == 'SampleSheet.csv' and run_info is None:
                    try:
                        with open(os.path.join(root, file), 'rU') as handle:
                            run_info = sample_sheet_parser(handle)
                    except:
                        raise

                if len(fastq_re.findall(file)) > 0:
                    fastq_files.append(os.path.join(root, file))


        if fastq_files:
            self.write('Found %d FASTQ files.\n' % (len(fastq_files),))
        else:  # empty list
            self.write('Error, folder does not seem to contain any FASTQ(.gz) files!\n')
            self.rundir = None

        if run_info is None:
            self.write('Warning, failed to locate run manifest (SampleSheet.csv).\n')

        return fastq_files

    def callback(self, msg):
        if type(msg) is int:
            self.progress_bar['value'] = msg
        elif type(msg) is str:
            self.write(msg+'\n')

        self.parent.update_idletasks()

    def make_tree(self, path):
        if not os.path.isdir(path):
            parent = os.path.dirname(path)
            self.make_tree(parent)
            os.mkdir(path)

    def process_files(self):
        """
        Perform MiCall data processing on FASTQ files in working directory.
        :return:
        """
        image_paths = []
        prefixes = []

        # look for FASTQ files
        fastq_files = self.open_files()
        if len(fastq_files) == 0:
            return

        setting_name = 'results_path'
        savedir = self.config.get(setting_name, '')
        savedir = tkFileDialog.askdirectory(title='Choose a folder to save results',
                                            initialdir=savedir)
        if not savedir:
            return
        self.config[setting_name] = savedir
        self.write_config()
        self.workdir = os.path.join(savedir, 'working')
        if os.path.exists(self.workdir):
            shutil.rmtree(self.workdir)
        self.make_tree(self.workdir)
        os.chdir(self.workdir)

        # transfer FASTQ.gz files to working folder
        self.target_files = []
        for src in fastq_files:
            if src.startswith(self.workdir):
                # Working file from previous run.
                continue
            filename = os.path.basename(src)
            prefix = filename.split('.')[0]
            dest = os.path.join(self.workdir, filename)

            if '_R1_001' in dest:
                self.target_files.append(dest.replace('.gz', ''))

            # neither file type is present
            if not dest.endswith('.gz'):
                shutil.copy(src, dest)
            else:
                dest = os.path.join(self.workdir, prefix+'.fastq')
                with gzip.open(src, 'rb') as zip_src, open(dest, 'w') as fastq_dest:
                    shutil.copyfileobj(zip_src, fastq_dest)


        # remove duplicate entries
        self.target_files = sorted(set(self.target_files))

        self.write('Transferred %d sets of FASTQ files to working directory.\n' % len(self.target_files))

        for fastq1 in self.target_files:
            fastq2 = fastq1.replace('_R1_001', '_R2_001')
            if not os.path.exists(fastq2):
                self.write('ERROR: Missing R2 file for', fastq1)
                continue

            prefix = os.path.basename(fastq1).replace('_L001_R1_001.fastq', '')
            prefixes.append(prefix)
            output_csv = fastq1.replace('_L001_R1_001.fastq', '.prelim.csv')

            self.write('Processing sample %s\n... preliminary mapping\n' % (prefix,))
            # four lines per read, two files
            nrecords = self.line_counter.count(fastq1) / 2
            self.progress_bar['value'] = 0
            self.progress_bar['maximum'] = nrecords
            self.parent.update_idletasks()  # flush buffer

            with open(output_csv, 'wb') as handle:
                prelim_map(fastq1,
                           fastq2,
                           handle,
                           nthreads=self.nthreads.get(),
                           callback=self.callback)

            # prepare file handles for remap stage
            with open(output_csv, 'rU') as prelim_csv, \
                 open(os.path.join(self.workdir, prefix+'.remap.csv'), 'wb') as remap_csv, \
                 open(os.path.join(self.workdir, prefix+'.remap_counts.csv'), 'wb') as counts_csv, \
                 open(os.path.join(self.workdir, prefix+'.remap_conseq.csv'), 'wb') as conseq_csv, \
                 open(os.path.join(self.workdir, prefix+'.unmapped1.fastq'), 'w') as unmapped1, \
                 open(os.path.join(self.workdir, prefix+'.unmapped2.fastq'), 'w') as unmapped2:

                self.write('... remapping\n')
                self.parent.update_idletasks()
                self.progress_bar['value'] = 0
                remap(fastq1,
                      fastq2,
                      prelim_csv,
                      remap_csv,
                      counts_csv,
                      conseq_csv,
                      unmapped1,
                      unmapped2,
                      self.workdir,
                      nthreads=self.nthreads.get(),
                      callback=self.callback)

            # prepare file handles for conversion from SAM format to alignment
            with open(os.path.join(self.workdir, prefix+'.remap.csv'), 'rU') as remap_csv, \
                 open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'wb') as aligned_csv, \
                 open(os.path.join(self.workdir, prefix+'.insert.csv'), 'wb') as insert_csv, \
                 open(os.path.join(self.workdir, prefix+'.failed.csv'), 'wb') as failed_csv:
                

                self.write('... converting into alignment\n')
                self.parent.update_idletasks()
                sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv, nthreads=self.nthreads.get())

            with open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'rU') as aligned_csv, \
                 open(os.path.join(self.workdir, prefix+'.nuc.csv'), 'wb') as nuc_csv, \
                 open(os.path.join(self.workdir, prefix+'.amino.csv'), 'wb') as amino_csv, \
                 open(os.path.join(self.workdir, prefix+'.coord_ins.csv'), 'wb') as coord_ins_csv, \
                 open(os.path.join(self.workdir, prefix+'.conseq.csv'), 'wb') as conseq_csv, \
                 open(os.path.join(self.workdir, prefix+'.failed_align.csv'), 'wb') as failed_align_csv, \
                 open(os.path.join(self.workdir, prefix+'.nuc_variants.csv'), 'wb') as nuc_variants_csv:

                self.write('... extracting statistics from alignments\n')
                self.parent.update_idletasks()
                aln2counts(aligned_csv,
                           nuc_csv,
                           amino_csv,
                           coord_ins_csv,
                           conseq_csv,
                           failed_align_csv,
                           nuc_variants_csv)

            self.write('... generating coverage plots\n')
            self.parent.update_idletasks()
            with open(os.path.join(self.workdir, prefix+'.amino.csv'), 'rU') as amino_csv:
                image_paths += coverage_plot(amino_csv)

            self.write('... performing g2p scoring on samples covering HIV-1 V3\n')
            self.parent.update_idletasks()
            with open(os.path.join(self.workdir, prefix+'.remap.csv'), 'rU') as remap_csv, \
                 open(os.path.join(self.workdir, prefix+'.nuc.csv'), 'rU') as nuc_csv, \
                 open(os.path.join(self.workdir, prefix+'.g2p.csv'), 'wb') as g2p_csv:
                
                sam_g2p(pssm=self.pssm, remap_csv=remap_csv, nuc_csv=nuc_csv, g2p_csv=g2p_csv)

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
        os.chdir(savedir)
        shutil.rmtree(self.workdir)

        self.write('Run complete.\n')


root = tk.Tk()  # parent widget
root.wm_title('MiCall')

app = MiCall(root).pack(fill='both', expand=True)

root.mainloop()  # enter Tk event loop
