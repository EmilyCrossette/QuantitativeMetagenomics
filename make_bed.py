import sys

import pandas as pd
import click

@click.command()
@click.argument('fasta_name')
@click.argument('fasta')
@click.option('--output_fp', '-o', default = 'genome.bed')

def make_bed(fasta_name, fasta, output_fp):
    """Obtain the start and stop positions from a multi-fasta downloaded from NCBI

    Args:
        fasta_name (str): name of chromosome (must aggree with genome-wide fasta header)
        fasta (str): /path/and/name/to.fa[.gz]
            (file may be gziped)
        output_fp (str): output file path (optional)

    Yields:
        BED file containing
            header (str): header sequence of fasta entry (excludes '>')
            start (int): position of the gene's start in genome
            stop (int): position of the gene's end in genome
            gene_name (str): name of gene

    This script was derived from a fasta sequence parser developed by M. Sherman.
    """

    if '.gz' in fasta:
        file_object = gzip.open(fasta, 'rt')
    else:
        file_object = open(fasta, 'rt')

    beds = []
    for line in file_object:
        if line.startswith('>'):
            gene_name = line[1:].split()[0]
            header = line[1:]
            if "complement" in header:
                # complement = "Y"
                gene_start = header.partition("complement(")[2].partition("..")[0]
                gene_end = header.partition("..")[2].partition(")]")[0]
            else:
                # complement = "N"
                gene_start = header.partition("location=")[2].partition("..")[0]
                gene_end = header.partition("..")[2].partition("] ")[0]

            beds.append({'header': fasta_name,
                         'start': gene_start,
                         'end': gene_end,
                         'gene_name': gene_name})

    output = pd.DataFrame(beds)
    output.to_csv(output_fp, sep='\t', index=None, header=None)


if __name__ == '__main__':
    make_bed()
