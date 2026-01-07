#!/usr/bin/env python3

import csv
import os
import subprocess as sp
import typing
from argparse import ArgumentParser
from collections import namedtuple
from pathlib import Path

from Bio import SeqIO
from xopen import xopen

import pyaragorn


class Gene(typing.NamedTuple):
    id: str
    kind: str
    start: int
    stop: int
    strand: int


def aragorn_predict(genome: Path, txt_file: Path, aragorn_bin: Path, linear: bool, transl_table: int = 11, ps: float = 100.0):
    """
    Predict CDS with aragorn
    :param genome: Path to the genome (fasta)
    :param train_file: Path to the training file
    :param txt_file: Path to text output file
    :param aragorn_bin: Path to the aragorn binary
    :param closed: Closed ends.
    :param transl_table: Translation table
    :param ps: Scoring threshold percentage (default 100.0)
    """
    cmd = [
        aragorn_bin,
        '-o', str(txt_file),
        '-wa',
        '-e',
        '-l' if linear else '-c',
        f'-ps{ps}',
        f'-gc{transl_table}',
        str(genome),
    ]

    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'aragorn error! error code: {proc.returncode}')


def parse_txt_output(txt_path: Path, pyaragorn_processed: bool = False) -> set[Gene]:
    """
    Parse the aragorn output.
    :param gff_path: Path to the GFF output file
    :param pyaragorn_processed: pyaragorn output
    :return: Set of namedtuples of CDS coordinates
    """
    cdss: set[Gene] = set()
    contig_id = None
    with txt_path.open() as fh:
        for line in fh:
            if line[0] == '>':
                contig_id = line[1:].split()[0]
                fh.readline()
            else:
                #print(line)
                (n, kind, coords, energy, x, y) = line.strip().split()
                if coords.startswith("c"):
                    strand = -1
                    coords = coords[1:]
                else:
                    strand = 1
                start, stop = map(int, coords[1:-1].split(","))
                gene: namedtuple = Gene(contig_id,
                                        kind,
                                        int(start),
                                        int(stop),
                                        strand)
                cdss.add(gene)
    return cdss



def pyaragorn_predict(rna_finder: pyaragorn.RNAFinder, genome: Path) -> set[Gene]:
    """
    Predict all genes using pyaragorn.
    :param orf_finder: Trained pyaragorn OrfFinder instance
    :param genome: Path to the contig (fasta)
    :return:
    """
    cdss: set[Gene] = set()
    with xopen(genome, mode='r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            for prediction in rna_finder.find_rna(bytes(record.seq)):
                kind = prediction.type
                if prediction.type == "tRNA":
                    if prediction.amino_acid == "???":
                        kind = f"{kind}-?({'|'.join(prediction.amino_acids)})"
                    else:
                        kind = f"{kind}-{prediction.amino_acid}"
                elif prediction.type == "tmRNA":
                    mask = [["", "p"], ["*", "p*"]]
                    kind = "tmRNA{}".format(mask[prediction.permuted][prediction.energy < 100.0])
                gene: namedtuple = Gene(record.id,
                                        kind,
                                        prediction.begin,
                                        prediction.end,
                                        prediction.strand)
                cdss.add(gene)
    return cdss


def compare_mismatches(aragorn_genes: set[Gene], pyaragorn_genes: set[Gene]) -> list:
    """
    Compare aragorn/pyaragorn mismatches to find related hits.
    :param aragorn_genes: Mismatched aragorn genes
    :param pyaragorn_genes: Mismatched pyaragorn genes
    :return: Related hits
    """
    ordered_pairs: list = []
    empty_gene = Gene(None, None, None, None, None)

    longer, shorter = (aragorn_genes, pyaragorn_genes) if len(aragorn_genes) >= len(pyaragorn_genes) else (pyaragorn_genes, aragorn_genes)
    aragorn_more_hits = True if len(aragorn_genes) >= len(pyaragorn_genes) else False

    for main_cds in longer:
        for candidate_cds in shorter:
            if main_cds.id == candidate_cds.id and main_cds.strand == candidate_cds.strand:  # matching contig AND strand
                if main_cds.stop == candidate_cds.stop or main_cds.start == candidate_cds.start:  # matching stop or start
                    mismatch = {'aragorn': main_cds,
                                'pyaragorn': candidate_cds,
                                'strand': main_cds.strand} if aragorn_more_hits else {'aragorn': candidate_cds,
                                                                                       'pyaragorn': main_cds,
                                                                                       'strand': main_cds.strand}
                    ordered_pairs.append(mismatch)
                    break
        else:
            mismatch = {'aragorn': main_cds,
                        'pyaragorn': empty_gene,
                        'strand': main_cds.strand} if aragorn_more_hits else {'aragorn': empty_gene,
                                                                               'pyaragorn': main_cds,
                                                                               'strand': main_cds.strand}
            ordered_pairs.append(mismatch)

    return ordered_pairs


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Compare CDS predictions of aragorn to the predictions of pyaragorn and save the differences in a TSV file.')
    parser.add_argument('--genome', '-g', type=Path, nargs='+', action='append', help='Input genomes (/some/path/*.fasta)')
    parser.add_argument('--aragorn', '-a', type=Path, default=Path('aragorn'),
                        help='Path to a newly compiled aragorn binary.')
    parser.add_argument('--linear', '-l', action='store_true', help='Closed ends. Do not allow genes to run off edges.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output path (default="./comparison")')
    parser.add_argument('--ps', '-p', type=float, default=100.0, help='Scoring threshold percentage (default=100.0)')
    args = parser.parse_args()
    return args.genome, args.aragorn, args.linear, args.output, args.ps


def main():
    genomes, aragorn_bin, linear, out_path, ps = parse_arguments()
    out_path.resolve()
    print(f'Genomes linear={linear} ps={ps}')

    tsv_file = open(out_path.joinpath('mismatches.tsv'), mode='w')
    tsv_writer = csv.writer(tsv_file, delimiter='\t', lineterminator='\n')
    header: list[str] = ['genome', 'strand', 'aragorn_kind', 'aragorn_start', 'aragorn_stop', 'pyaragorn_kind', 'pyaragorn_start', 'pyaragorn_stop']
    tsv_writer.writerow(header)

    for genome in genomes[0]:
        genome.resolve()
        prefix: str = str(genome).split('/')[-1].split('.')[0]

        tmp_path = out_path.joinpath(f'tmp/{prefix}')
        try:
            os.makedirs(tmp_path)
        except FileExistsError:
            pass

        # Predict aragorn
        aragorn_txt_file = tmp_path.joinpath(f'{prefix}.aragorn.txt')
        aragorn_predict(genome=genome,
                         txt_file=aragorn_txt_file,
                         aragorn_bin=aragorn_bin,
                         linear=linear,
                         ps=ps)
        aragorn_genes: set[dict] = parse_txt_output(txt_path=aragorn_txt_file)
        # Predict pyaragorn
        rna_finder = pyaragorn.RNAFinder(11, linear=linear, threshold_scale=ps/100.0)
        pyaragorn_genes: set[dict] = pyaragorn_predict(rna_finder=rna_finder,
                                                       genome=genome)

        print(f'Hits genome={prefix}: aragorn={len(aragorn_genes)}, pyaragorn={len(pyaragorn_genes)}, equal={aragorn_genes == pyaragorn_genes}')

        if aragorn_genes != pyaragorn_genes:
            difference_aragorn = aragorn_genes.difference(pyaragorn_genes)
            difference_pyaragorn = pyaragorn_genes.difference(aragorn_genes)
            ordered_pairs = compare_mismatches(difference_aragorn, difference_pyaragorn)

            for mismatch in ordered_pairs:
                # genome, strand, aragorn_start, aragorn_stop, aragorn_edge, pyaragorn_start, pyaragorn_stop, pyaragorn_edge
                line = [prefix,
                        mismatch['strand'],
                        mismatch['aragorn'].kind,
                        mismatch['aragorn'].start,
                        mismatch['aragorn'].stop,
                        mismatch['pyaragorn'].kind,
                        mismatch['pyaragorn'].start,
                        mismatch['pyaragorn'].stop]

                tsv_writer.writerow(line)

    tsv_file.close()


if __name__ == '__main__':
    main()
