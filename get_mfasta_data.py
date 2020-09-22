import pandas as pd
import screed
import re

import click
import os

from tqdm import tqdm


@click.command()
@click.argument('multifasta_file')
@click.option('--output_fp', '-o', default='gene_info.tsv')

def get_mfasta_data(multifasta_file, output_fp):
    """Produce a table with the target ID, length and GC content

    Args:
        multifasta
        tsv_name (str): optional file name including ".tsv"

    Returns:
        TSV file:
            target (str): header sequence of fasta entry (excludes '>')
            Gene_len (int): length of sequence
            GC-content (float): GC-content of sequence

    Psuedocode:

    Capture header and append to table
    Count number of letters in line below header
    Count number of of GC content

    """


    with screed.open(multifasta_file) as seqfile:
        output = []
        for seq in seqfile:
            output.append({
                'target':      seq.name,
                'gc_content':  len(re.findall('[GCgc]', seq.sequence)) / len(seq.sequence),
                'gene_length': len(seq.sequence),
            })

    df = pd.DataFrame(output)

    df.to_csv(output_fp, sep='\t', index=None)

if __name__ == '__main__':
    get_mfasta_data()
