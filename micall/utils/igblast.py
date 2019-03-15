"""wrapper for igBlast
"""
from subprocess import check_output
from tempfile import NamedTemporaryFile

default = {'igblast_path':\
    '/bin/igblastn',\
    'db_path':\
        '/data/ncbi-igblast-1.13.0/'}

def igblast(seq, config=default):
    result = None 
    with NamedTemporaryFile() as fasta_in:
        fasta_fd = open(fasta_in.name, 'w')
        print(">sample\n{}".format(seq), file=fasta_fd)
        print(fasta_in.name)
        fasta_fd.close()
        db = config['db_path']
        p = [config['igblast_path'], '-germline_db_V', db + 'bin/tcr_v',\
            '-germline_db_J', db + 'bin/tcr_j', '-germline_db_D', db + 'bin/tcr_d',\
            '-organism', 'human', '-query', fasta_in.name, '-auxiliary_data',\
            db + 'optional_file/human_gl.aux', '-ig_seqtype', 'TCR',\
            '-show_translation', '-outfmt', '3']
        print(' '.join(p))

        result = check_output(p) 
    return result

def parse_fmt3(result):
    r = result.split()
    print(len(r))
    print(r[10:20])
