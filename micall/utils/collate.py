import glob
import os

# TODO: This can probably be deleted. Wait until we revisit cross-contamination filter to be sure.
# def collate_frequencies (run_path, output_path, output_type):
#     """
#     Collate amino/nuc .freq files into a single summary file.
#     """
#
#     if output_type == "amino":
#         file_extension = "amino.csv"
#         header = "sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*"
#     elif output_type == "nuc":
#         file_extension = "nuc.csv"
#         header = "sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T"
#     else:
#         raise Exception("Incorrect output_type parameter")
#
#     # Summarize uncleaned freqs
#     with open(output_path, "w") as f_out:
#         f_out.write("{}\n".format(header))
#         for path in [x for x in glob.glob("{}/*.{}".format(run_path, file_extension))
#                      if "clean.{}".format(file_extension) not in x]:
#             prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
#             sample = prefix.split(".")[0]
#             with open(path,"r") as f:
#                 lines = f.readlines()
#
#             for line in lines:
#                 f_out.write("{},{}\n".format(sample, line.rstrip("\n")))
#
#     # Summarize clean freqs
#     file_extension = "clean.{}".format(file_extension)
#
#     output_path = output_path.replace("_frequencies.csv", "_cleaned_frequencies.csv")
#     with open(output_path, "w") as f_out:
#         f_out.write("{}\n".format(header))
#         for path in glob.glob("{}/*.{}".format(run_path,file_extension)):
#             prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
#             sample, region, q = prefix.split(".")[:3]
#             with open(path,"r") as f:
#                 lines = f.readlines()
#
#             for j, line in enumerate(lines):
#                 if j == 0:
#                     continue
#                 f_out.write("{},{},{},{}\n".format(sample, region, q, line.rstrip("\n")))


def collate_named_files(src_dir, sample_list, extension, output_path):
    """
    Collate files
    :param src_dir: directory to copy files from
    :param sample_list: sample prefix to match by wildcard
    :param extension: filetype to collate (e.g., '.nuc.csv')
    :param output_path: path to write results to
    :return:
    """
    with open(output_path, 'w') as fout:
        is_header_written = False
        for samplename in sample_list:
            srcfile = os.path.join(src_dir, samplename+extension)
            if not os.path.exists(srcfile):
                continue
            with open(srcfile, 'rU') as fin:
                for i, line in enumerate(fin):
                    if i == 0:
                        if not is_header_written:
                            fout.write('sample,')
                            fout.write(line)
                            is_header_written = True
                    else:
                        fout.write(samplename + ',')
                        fout.write(line)


def collate_labeled_files(pattern, output_path):
    """ Collate files that are already labeled on each row.

    Just remove redundant headers and concatenate all the files together.
    @param pattern: search pattern for all the files to collate
    @param param: output_path the path to write the results to
    """
    with open(output_path, 'w') as fout:
        is_header_written = False
        filenames = glob.glob(pattern)
        filenames.sort()
        for filename in filenames:
            basename = os.path.basename(filename)
            samplename = basename.split('.')[0]
            with open(filename, 'rU') as fin:
                for i, line in enumerate(fin):
                    if i == 0:
                        if not is_header_written:
                            fout.write('sample,')
                            fout.write(line)
                            is_header_written = True
                    else:
                        fout.write(samplename + ',')
                        fout.write(line)
