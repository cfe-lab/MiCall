import os
import Tkinter as tk
import sys
import tkFileDialog
import shutil
import subprocess
from micall.utils.sample_sheet_parser import sample_sheet_parser
from micall.core.prelim_map import prelim_map
from micall.core.remap import remap
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts
from tempfile import gettempdir

class MiCall(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        self.__version__ = 0.1

        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.rundir = None
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
        self.button_load.grid(row=0, column=1, sticky='W')

        self.button_run = tk.Button(
            self.button_frame, text="Run", command=self.process_files, state=tk.DISABLED
        )
        self.button_run.grid(row=0, column=2, sticky='W')

        self.button_quit = tk.Button(
            self.button_frame, text='Quit', command=self.quit
        )
        self.button_quit.grid(row=0, column=3, sticky='W')

        self.console = tk.Text(
            self.console_frame, bg='black', fg='white', cursor='xterm'
        )
        self.console.pack(side='top', fill='both', expand=True)
        self.write('Welcome to MiCall v%1.1f!\n' % (self.__version__,))
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
        self.button_load.config(state=tk.DISABLED)

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
                if file.endswith('_R1_001.fastq.gz'):
                    fastq_files.append(os.path.join(root, file))

        if not fastq_files:  # empty list
            self.write('Error, folder does not seem to contain any FASTQ.gz files!\n')
            self.rundir = None
            #self.button_setwd.config(state=tk.ACTIVE)
            self.button_load.config(state=tk.ACTIVE)
            return False

        if run_info is None:
            self.write('Warning, failed to locate run manifest (SampleSheet.csv).\n')

        self.write('Found %d sets of FASTQ files.\n' % (len(fastq_files),))

        # transfer FASTQ.gz files to working folder
        self.target_files = []
        for src in fastq_files:
            filename = os.path.basename(src)
            dest = os.path.join(self.workdir, filename)
            self.target_files.append(dest.replace('.gz', ''))

            # don't overwrite file if already present
            if os.path.exists(os.path.join(self.workdir, filename.replace('.gz', ''))):
                self.write('Found %s in working folder, skipping transfer.\n' % (filename,))
                continue

            # copy file
            #self.write('Copying from %s to %s\n' % (src, dest))
            shutil.copy(src, dest)

            src2 = src.replace('_R1_001', '_R2_001')
            dest2 = os.path.join(self.workdir, os.path.basename(src2))
            shutil.copy(src2, dest2)

            # decompress files
            #self.write('Uncompressing files\n')
            p = subprocess.Popen(['gunzip', dest], stdout=subprocess.PIPE)
            p = subprocess.Popen(['gunzip', dest2], stdout=subprocess.PIPE)

        self.write('Transferred and uncompressed FASTQ files to working directory.\n')
        self.button_run.config(state=tk.ACTIVE)


    def process_files(self):
        """
        Perform MiCall data processing on FASTQ files in working directory.
        :return:
        """
        working_files = []
        output_files = []

        # look for FASTQ files
        if len(self.target_files) == 0:
            print 'ERROR: No files to process'
            return

        for fastq1 in self.target_files:
            fastq2 = fastq1.replace('_R1_001', '_R2_001')
            prefix = os.path.basename(fastq1).replace('_L001_R1_001.fastq', '')
            output_csv = fastq1.replace('_L001_R1_001.fastq', '.prelim.csv')
            working_files.append(output_csv)

            self.write('Preliminary map of %s\n' % (prefix,))
            self.parent.update_idletasks()  # flush buffer
            with open(output_csv, 'w') as handle:
                prelim_map(fastq1, fastq2, handle, self.workdir)

            # prepare file handles for remap stage
            prelim_csv = open(output_csv, 'rU')
            remap_csv = open(os.path.join(self.workdir, prefix+'.remap.csv'), 'w')
            counts_csv = open(os.path.join(self.workdir, prefix+'.remap_counts.csv'), 'w')
            conseq_csv = open(os.path.join(self.workdir, prefix+'.remap_conseq.csv'), 'w')
            unmapped1 = open(os.path.join(self.workdir, prefix+'.unmapped1.fastq'), 'w')
            unmapped2 = open(os.path.join(self.workdir, prefix+'.unmapped2.fastq'), 'w')
            working_files += map(lambda x: x.name, [remap_csv, counts_csv, conseq_csv, unmapped1, unmapped2])

            self.write('Iterative remap of %s\n' % (prefix,))
            self.parent.update_idletasks()
            remap(fastq1, fastq2, prelim_csv, remap_csv, counts_csv, conseq_csv, unmapped1, unmapped2, self.workdir)

            # prepare file handles for conversion from SAM format to alignment
            remap_csv = open(os.path.join(self.workdir, prefix+'.remap.csv'), 'rU')
            aligned_csv = open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'w')
            insert_csv = open(os.path.join(self.workdir, prefix+'.insert.csv'), 'w')
            failed_csv = open(os.path.join(self.workdir, prefix+'.failed.csv'), 'w')
            working_files += map(lambda x: x.name, [aligned_csv, insert_csv, failed_csv])

            self.write('Converting mapped reads into alignment\n')
            self.parent.update_idletasks()
            sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv)

            aligned_csv = open(os.path.join(self.workdir, prefix+'.aligned.csv'), 'rU')
            nuc_csv = open(os.path.join(self.workdir, prefix+'.nuc.csv'), 'w')
            amino_csv = open(os.path.join(self.workdir, prefix+'.amino.csv'), 'w')
            coord_ins_csv = open(os.path.join(self.workdir, prefix+'.coord_ins.csv'), 'w')
            conseq_csv = open(os.path.join(self.workdir, prefix+'.conseq.csv'), 'w')
            failed_align_csv = open(os.path.join(self.workdir, prefix+'.failed_align.csv'), 'w')
            nuc_variants_csv = open(os.path.join(self.workdir, prefix+'.nuc_variants.csv'), 'w')
            output_files += map(lambda x: x.name, [nuc_csv, amino_csv, coord_ins_csv, conseq_csv, failed_align_csv,
                                                    nuc_variants_csv])

            self.write('Extracting statistics from alignments\n')
            self.parent.update_idletasks()
            aln2counts(aligned_csv, nuc_csv, amino_csv, coord_ins_csv, conseq_csv, failed_align_csv, nuc_variants_csv,
                       self.workdir)

        # prevent rerun until a new folder is loaded
        self.button_run.config(state=tk.DISABLED)

        # collate results and specify save location, then reset buttons
        savedir = tkFileDialog.askdirectory(title='Select folder to save results')
        for src in output_files:
            dest = os.path.join(savedir, os.path.basename(src))
            shutil.move(src, dest)

        # clean up working files
        for fn in working_files:
            os.remove(fn)

        # clean up raw data
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
