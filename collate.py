import glob, os

def collate_conseqs(run_path,output_path):
    """
    Collate .conseq.csv files into a single CSV summary file.
    """

    files = glob.glob("{}/*.conseq.csv".format(run_path))

    with open(output_path,"w") as f_out:
        f_out.write("sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence\n")

        for path in files:
            prefix = (os.path.basename(path)).rstrip(".conseq.csv")
            sample = prefix.split(".")[0]
            sname, snum = sample.split('_')
            with open(path,"r") as f:
                for line in f:
                    region, qcut, cut, conseq = line.rstrip().split(',')
                    f_out.write(','.join((sname,
                                          region,
                                          qcut,
                                          snum,
                                          cut,
                                          conseq)) + '\n')

def collate_frequencies (run_path, output_path, output_type):
    """
    Collate amino/nuc .freq files into a single summary file.
    """

    if output_type == "amino":
        file_extension = "amino.csv"
        header = "sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*"
    elif output_type == "nuc":
        file_extension = "nuc.csv"
        header = "sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T"
    else:
        raise Exception("Incorrect output_type parameter")

    # Summarize uncleaned freqs
    with open(output_path, "w") as f_out:
        f_out.write("{}\n".format(header))
        for path in [x for x in glob.glob("{}/*.{}".format(run_path,file_extension)) if "clean.{}".format(file_extension) not in x]:
            prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
            sample = prefix.split(".")[0]
            with open(path,"r") as f:
                lines = f.readlines()

            for line in lines:
                f_out.write("{},{}\n".format(sample, line.rstrip("\n")))

    # Summarize clean freqs
    file_extension = "clean.{}".format(file_extension)

    output_path = output_path.replace("_frequencies.csv", "_cleaned_frequencies.csv")
    with open(output_path, "w") as f_out:
        f_out.write("{}\n".format(header))
        for path in glob.glob("{}/*.{}".format(run_path,file_extension)):
            prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
            sample, region, q = prefix.split(".")[:3]
            with open(path,"r") as f:
                lines = f.readlines()

            for j, line in enumerate(lines):
                if j == 0:
                    continue
                f_out.write("{},{},{},{}\n".format(sample, region, q, line.rstrip("\n")))

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
            with open(filename, 'rU') as fin:
                for i, line in enumerate(fin):
                    if i == 0:
                        if not is_header_written:
                            fout.write(line)
                            is_header_written = True
                    else:
                        fout.write(line)
