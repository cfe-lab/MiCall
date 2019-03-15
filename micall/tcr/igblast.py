"""wrapper for igBlast
"""
import sys
from subprocess import check_output
from tempfile import NamedTemporaryFile
#import externals, os

# make sure $IGDATA is set to this
default_path = '/opt/igblast/'

# default_db = __file__ + "/micall/data"
default_db = '/opt/micall/micall/data/'


def igblast_seq(seq, db=default_db, path=default_path):
    with NamedTemporaryFile() as fasta_in:
        fasta_fd = open(fasta_in.name, 'w')
        print(">micall_tempsample\n{}".format(seq), file=fasta_fd)
        fasta_fd.close()
        return igblast(fasta_in.name, db=db, path=path)

def igblast(fasta, db=default_db, path=default_path):
    p = [path + "bin/igblastn", '-germline_db_V', db + 'tcr_v',\
        '-germline_db_J', db + 'tcr_j', '-germline_db_D', db + 'tcr_d',\
        '-organism', 'human', '-query', fasta, '-auxiliary_data',\
         path + 'optional_file/human_gl.aux', '-ig_seqtype', 'TCR',\
        '-show_translation', '-outfmt', '3']
    return check_output(p) 

def parse_fmt3(result):
    queries = []
    for query in result.decode().split('Query= ')[1:]:
        sample = query.split()[0]
        scores = []
        for line in query.splitlines():
            line = line.split()
            if "Domain" in line:
                break
            if len(line) == 3:
                if line[0] == 'Lambda':
                    break
                scores.append((line[0], float(line[1]), float(line[2])))
        queries.append((sample, scores))
    return queries

if __name__ == "__main__":
    print(parse_fmt3(igblast_seq(sys.argv[1])))
