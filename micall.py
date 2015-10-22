import gzip
import json
import multiprocessing
from multiprocessing import cpu_count
import os
import re
import shutil
import sys
from tempfile import gettempdir
import traceback

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

fastq_re = re.compile('(_L001_R[12]_001|.censored[12]).fastq')
sys.stdin = open(os.devnull, 'r') # Fixes a problem when Windows launches without a console.

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
        self.pssm = Pssm(path_to_lookup=AssetWrapper('micall/g2p/g2p_fpr.txt').path,
                         path_to_matrix=AssetWrapper('micall/g2p/g2p.matrix').path)

        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        parent.report_callback_exception = self.report_callback_exception

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

        self.nthreads = self.config.get('threads', None)
        if not self.nthreads:
            self.nthreads = int(round(cpu_count() * 0.5))
            self.config['threads'] = self.nthreads
            self.write_config()

        self.button_run = tk.Button(
            self.button_frame, text="Run", command=self.process_files
        )
        self.button_run.grid(row=0, column=1, sticky='W')

        self.progress_bar = Progressbar(self.button_frame, orient='horizontal', length=500, mode='determinate')
        self.progress_bar.grid(row=1, columnspan=5)

        scrollbar = tk.Scrollbar(self.console_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.console = tk.Text(self.console_frame,
                               bg='black',
                               fg='white',
                               yscrollcommand=scrollbar.set)
        self.console.pack(side=tk.LEFT, fill=tk.BOTH)
        self.console.tag_configure('ERROR', foreground="red")
        scrollbar.config(command=self.console.yview)

        # redirect stderr to Text widget
        #sys.stderr = Redirector(self.console)

        self.write('Welcome to MiCall v{}, running with {} threads.\n'.format(
            pipeline_version,
            self.nthreads))

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
                json.dump(self.config, f, indent=4)
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

        # check for presence of FASTQ files
        for root, _dirs, files in os.walk(self.rundir):
            for file in files:
                if len(fastq_re.findall(file)) > 0:
                    fastq_files.append(os.path.join(root, file))

        if fastq_files:
            self.write('Found %d FASTQ files.\n' % (len(fastq_files),))
        else:  # empty list
            self.write('Error, folder does not seem to contain any FASTQ(.gz) files!\n')
            self.rundir = None

        return fastq_files

    def callback(self, message=None, progress=None, max_progress=None):
        """ Report progress to the user.
        
        @param message: a string to write in the window, or None to not write
        @param progress: an integer value to update the progress bar
        @param max_progress: an integer value that represents completed progress
        """
        if max_progress:
            self.progress_bar['maximum'] = max_progress
        if progress is not None:
            self.progress_bar['value'] = progress
        if message is not None:
            self.write(message + '\n')

        self.parent.update()

    def make_tree(self, path):
        if not os.path.isdir(path):
            parent = os.path.dirname(path)
            self.make_tree(parent)
            os.mkdir(path)


    def process_sample(self, fastq1, progress, prefixes, image_paths, error_log):
        fastq2 = fastq1.replace('_R1_001', '_R2_001').replace('censored1',
                                                              'censored2')
        if not os.path.exists(fastq2):
            raise IOError('ERROR: Missing R2 file for {}'.format(fastq1))
        
        prefix = os.path.basename(fastq1).replace('_L001_R1_001.fastq',
                                                  '').replace('.censored1.fastq',
                                                              '')
        prefixes.append(prefix)
        output_csv = prefix + '.prelim.csv'
        self.write('Processing sample {} ({})\n'.format(prefix, progress))
        with open(output_csv, 'wb') as handle:
            prelim_map(fastq1,
                       fastq2,
                       handle,
                       nthreads=self.nthreads,
                       callback=self.callback,
                       stderr=error_log)
        
        # prepare file handles for remap stage
        with open(output_csv, 'rU') as prelim_csv, \
             open(os.path.join(self.workdir, prefix + '.remap.csv'), 'wb') as remap_csv, \
             open(os.path.join(self.workdir, prefix + '.remap_counts.csv'), 'wb') as counts_csv, \
             open(os.path.join(self.workdir, prefix + '.remap_conseq.csv'), 'wb') as conseq_csv, \
             open(os.path.join(self.workdir, prefix + '.unmapped1.fastq'), 'w') as unmapped1, \
             open(os.path.join(self.workdir, prefix + '.unmapped2.fastq'), 'w') as unmapped2:
            
            self.write('... remapping\n')
            self.parent.update()
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
                  nthreads=self.nthreads,
                  callback=self.callback,
                  stderr=error_log)
            
        # prepare file handles for conversion from SAM format to alignment
        with open(os.path.join(self.workdir, prefix + '.remap.csv'), 'rU') as remap_csv, \
             open(os.path.join(self.workdir, prefix + '.aligned.csv'), 'wb') as aligned_csv, \
             open(os.path.join(self.workdir, prefix + '.insert.csv'), 'wb') as insert_csv, \
             open(os.path.join(self.workdir, prefix + '.failed.csv'), 'wb') as failed_csv:
            
            self.write('... converting into alignment\n')
            self.parent.update()
            sam2aln(remap_csv,
                    aligned_csv,
                    insert_csv,
                    failed_csv,
                    nthreads=self.nthreads)
            
        with open(os.path.join(self.workdir, prefix + '.aligned.csv'), 'rU') as aligned_csv, \
             open(os.path.join(self.workdir, prefix + '.nuc.csv'), 'wb') as nuc_csv, \
             open(os.path.join(self.workdir, prefix + '.amino.csv'), 'wb') as amino_csv, \
             open(os.path.join(self.workdir, prefix + '.coord_ins.csv'), 'wb') as coord_ins_csv, \
             open(os.path.join(self.workdir, prefix + '.conseq.csv'), 'wb') as conseq_csv, \
             open(os.path.join(self.workdir, prefix + '.failed_align.csv'), 'wb') as failed_align_csv, \
             open(os.path.join(self.workdir, prefix + '.nuc_variants.csv'), 'wb') as nuc_variants_csv:
            
            self.parent.update()
            aln2counts(aligned_csv,
                       nuc_csv,
                       amino_csv,
                       coord_ins_csv,
                       conseq_csv,
                       failed_align_csv,
                       nuc_variants_csv,
                       callback=self.callback)
            
        self.write('... generating coverage plots\n')
        self.parent.update()
        with open(os.path.join(self.workdir, prefix + '.amino.csv'), 'rU') as amino_csv:
            image_paths += coverage_plot(amino_csv)
        self.write('... performing g2p scoring on samples covering HIV-1 V3\n')
        self.parent.update()
        with open(os.path.join(self.workdir, prefix + '.remap.csv'), 'rU') as remap_csv, \
             open(os.path.join(self.workdir, prefix + '.nuc.csv'), 'rU') as nuc_csv, \
             open(os.path.join(self.workdir, prefix + '.g2p.csv'), 'wb') as g2p_csv:
            
            sam_g2p(pssm=self.pssm,
                    remap_csv=remap_csv,
                    nuc_csv=nuc_csv,
                    g2p_csv=g2p_csv)
    
    def build_run_summary(self, sample_count, error_count):
        summary = '{} sample'.format(sample_count)
        if sample_count != 1:
            summary += 's'
        if error_count > 0:
            summary += ' with {} error'.format(error_count)
            if error_count > 1:
                summary += 's'
        return summary

    def process_files(self):
        """
        Perform MiCall data processing on FASTQ files in working directory.
        """
        # look for FASTQ files
        fastq_files = self.open_files()
        if len(fastq_files) == 0:
            return

        savedir_setting_name = 'results_path'
        savedir = self.config.get(savedir_setting_name, '')
        savedir = tkFileDialog.askdirectory(title='Choose a folder to save results',
                                            initialdir=savedir)
        if not savedir:
            return
        self.config[savedir_setting_name] = savedir
        if os.path.exists(savedir) and os.listdir(savedir):
            self.write('Run FAILED - results folder is not empty: {}\n'.format(
                savedir))
            return
        
        workdir_purged_setting_name = 'workdir_purged'
        is_workdir_purged = self.config.get(workdir_purged_setting_name)
        if is_workdir_purged is None:
            self.config[workdir_purged_setting_name] = is_workdir_purged = True
        self.write_config()
        
        self.workdir = os.path.join(savedir, 'working')
        if os.path.exists(self.workdir):
            shutil.rmtree(self.workdir)
        self.make_tree(self.workdir)
        os.chdir(self.workdir)

        # transfer FASTQ.gz files to working folder
        self.callback(progress=0, max_progress=len(fastq_files))
        self.target_files = []
        for file_index, src in enumerate(fastq_files):
            if src.startswith(self.workdir):
                # Working file from previous run.
                continue
            self.callback(progress=file_index)
            filename = os.path.basename(src)
            prefix = filename.split('.')[0]
            dest = os.path.join(self.workdir, filename)

            if '_R1_001' in dest or 'censored1' in dest:
                self.target_files.append(dest.replace('.gz', ''))

            # neither file type is present
            if not dest.endswith('.gz'):
                shutil.copy(src, dest)
            else:
                dest = os.path.join(self.workdir, prefix+'.fastq')
                with gzip.open(src, 'rb') as zip_src, open(dest, 'w') as fastq_dest:
                    shutil.copyfileobj(zip_src, fastq_dest)
        self.callback(progress=len(fastq_files))

        # remove duplicate entries
        self.target_files = sorted(set(self.target_files))

        self.write('Transferred %d sets of FASTQ files to working directory.\n' % len(self.target_files))

        image_paths = []
        prefixes = []
        error_count = 0
        sample_count = len(self.target_files)
        run_summary = self.build_run_summary(sample_count, error_count)
        with open(os.path.join(savedir, 'micall.log'), 'w') as error_log:
            for sample_index, fastq1 in enumerate(self.target_files):
                try:
                    progress = '{} of {}'.format(sample_index+1, run_summary)
                    self.process_sample(fastq1, progress, prefixes, image_paths, error_log)
                except:
                    msg = traceback.format_exc()
                    self.write(msg)
                    error_log.write(msg)
                    error_count += 1
                    run_summary = self.build_run_summary(sample_count, error_count)

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
        if is_workdir_purged:
            shutil.rmtree(self.workdir)

        self.write('Run complete: {}.\n'.format(run_summary))

    def report_callback_exception(self, exc_type, exc_value, exc_traceback):
        message = ''.join(traceback.format_exception(exc_type,
                                                     exc_value,
                                                     exc_traceback))
        self.write(message)

def main():
    root = tk.Tk()  # parent widget
    root.wm_title('MiCall')
    
    MiCall(root).pack(fill='both', expand=True)
    
    root.mainloop()  # enter Tk event loop

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
