from prop_x4 import convert_fasta
def collate_counts(run_path,output_path):
    """
    Collate .counts files into a single CSV summary file.
    """
    import glob,os

    with open(output_path, "w") as collated_file:
        for count_file in [f for f in glob.glob("{}/*.counts".format(run_path))]:
            prefix = (os.path.basename(count_file)).rstrip(".counts")
            with open(count_file,"r") as f_in:
                for _, line in enumerate(f_in.readlines()):
                    collated_file.write("{},{}".format(prefix,line))

def collate_conseqs(run_path,output_path):
    """
    Collate .conseq files into a single CSV summary file.
    """
    import glob,os

    files = [f for f in glob.glob("{}/*.conseq".format(run_path)) if 'pileup' not in f]

    with open(output_path,"w") as f_out:
        f_out.write("sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence\n")

        for path in files:
            prefix = (os.path.basename(path)).rstrip(".conseq")
            _, region, q = prefix.split(".")[:3]

            with open(path,"r") as f:
                fasta = convert_fasta(f.readlines())

            for header, sequence in fasta:
                fasta_sample, s_number, consensus_percentage = header.split("_")
                f_out.write("{},{},{},{},{},{}\n".format(fasta_sample,region,q,s_number,consensus_percentage,sequence))

def collate_frequencies (run_path, output_path, output_type):
    """
    Collate amino/nuc .freq files into a single summary file.
    """
    import glob,os

    if output_type == "amino":
        file_extension = "amino.freqs"
        header = "sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*"
    elif output_type == "nuc":
        file_extension = "nuc.freqs"
        header = "sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T"
    else:
        raise Exception("Incorrect output_type parameter")

    # Summarize uncleaned freqs
    with open(output_path, "w") as f_out:
        f_out.write("{}\n".format(header))
        for path in [x for x in glob.glob("{}/*.{}".format(run_path,file_extension)) if "clean.{}".format(file_extension) not in x]:
            prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
            sample, region, q = prefix.split(".")[:3]
            with open(path,"r") as f:
                lines = f.readlines()

            for j, line in enumerate(lines):
                if j == 0:
                    continue
                f_out.write("{},{},{},{}\n".format(sample, region, q, line.rstrip("\n")))

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
