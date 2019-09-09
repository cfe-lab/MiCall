from sierralocal.main import scorefile
from sierralocal.hivdb import HIVdb
from sierralocal.jsonwriter import JSONWriter

from csv import DictReader, DictWriter
import sys
import tempfile
import argparse


def score_conseq(handle):
    reader = DictReader(handle)
    if 'consensus-percent-cutoff' not in reader.fieldnames:
        print("Error: input CSV does not appear to be a MiCall conseq output file.")
        sys.exit()

    algorithm = HIVdb()
    _ = JSONWriter(algorithm)

    # write sequences to temporary FASTA file
    tf = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for row in reader:
        tf.write('>{}\n{}\n'.format(row['consensus-percent-cutoff'], row['sequence']))

    sequence_headers, sequence_scores, ordered_mutation_list, file_genes, sequence_lengths, \
    file_trims, subtypes = scorefile(tf.name, algorithm)

    for i in range(len(sequence_headers)):
        header = sequence_headers[i]

        scores = sequence_scores[i]
        mutations = ordered_mutation_list[i]
        genes = file_genes[i]

        score_out = {}
        mutout = {}
        for j, gene in enumerate(genes):
            mutlist = ','.join('{}{}{}'.format(wt, pos, mut) for pos, mut, wt in mutations[j])
            mutout.update({gene[0]: mutlist})

            d = dict([(k, v[0]) for k, v in scores[j].items()])
            score_out.update({gene[0]: d})

        yield header, score_out, mutout


def main():
    parser = argparse.ArgumentParser(
        description="Uses sierra-local to generate HIV drug resistance predictions from "
                    "consensus sequences produced by MiCall-Lite, using the HIVdb algorithm.")

    parser.add_argument('conseq', type=argparse.FileType('r'),
                        help='<input> *.conseq.csv file from MiCall run.')
    parser.add_argument('out', type=argparse.FileType('w'),
                        help='<output> path to write tab-separated output.')

    args = parser.parse_args()

    results = score_conseq(args.conseq)

    # prepare output file based on first row of results
    header, scores, mutlist = next(results)

    fieldnames = ['header']
    genes = mutlist.keys()
    for gene in genes:
        for drug in scores[gene].keys():
            fieldnames.append(drug)
        fieldnames.append(gene+'.mutations')

    writer = DictWriter(args.out, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    # handle first entry
    row = {'header': header}
    for gene in genes:
        row.update(scores[gene])
        row.update({gene + '.mutations': mutlist[gene]})
    writer.writerow(row)

    # iterate through rest of results
    for header, scores, mutlist in results:
        row = {'header': header}
        for gene in genes:
            row.update(scores[gene])
            row.update({gene+'.mutations': mutlist[gene]})

        writer.writerow(row)


if __name__ == '__main__':
    main()
